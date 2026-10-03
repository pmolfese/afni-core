// PUBLIC DOMAIN NOTICE
//
// This file is part of afni-core, which was written by employees of the United
// States Government (National Institutes of Health) as part of their official
// duties. It is a "United States Government Work" (17 U.S.C. 105) and is in the
// public domain; outside the US, rights are waived under CC0 1.0. See LICENSE.
//
// ---------------------------------------------------------------------------
// WHAT THIS FILE IS
//
// Pure operations on sets of surface nodes (ROIs): grow and shrink a region by
// rings of neighbors or by distance along the surface, find its boundary, split
// it into connected pieces and clean up small ones, find the shortest path between
// two nodes, close an open path, and fill the area a closed path encloses.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `roi.rs` defines `NodeSet` (a sorted set of node indices) and `Roi`; these
//   functions take and return `NodeSet`s so they compose: `dilate(erode(x))`.
// * `topology.rs` supplies the connectivity (`SurfaceTopology`); `mesh.rs` adds the
//   coordinates (`SurfaceMesh`) needed when a distance is involved.
// * `roi_edit.rs` builds undoable edit commands out of these operations. Nothing
//   here knows about a mouse, a window or a file.
//
// CONVENTIONS
//
// * "Neighbors" are nodes joined by a mesh edge. A "ring" is one edge step.
// * Distances are shortest paths ALONG mesh edges (Dijkstra), in the mesh's
//   coordinate units; that is the surface distance SUMA draws and measures with.
// * Every function is deterministic: ties are broken by the smaller node index.
// * A node set may contain nodes in several disconnected pieces; the operations
//   work on all of them.
// ---------------------------------------------------------------------------

//! Operations on node sets: grow, shrink, boundary, components, paths, fill.

use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};

use crate::error::{Error, Result};
use crate::mesh::SurfaceMesh;
use crate::numeric::ensure_finite;
use crate::roi::NodeSet;
use crate::topology::SurfaceTopology;

// ---------------------------------------------------------------------------
// Boundary, grow and shrink by rings
// ---------------------------------------------------------------------------

/// Check every member is a node of the surface (and return the node count).
fn check_set(topology: &SurfaceTopology, set: &NodeSet) -> Result<usize> {
    let n = topology.node_count();
    set.validate_within(n)?;
    Ok(n)
}

/// The members that have at least one neighbor outside the set (the inner edge of
/// the region). A member at the rim of an open surface is not a boundary node for
/// that reason alone.
pub fn boundary_nodes(topology: &SurfaceTopology, set: &NodeSet) -> Result<NodeSet> {
    check_set(topology, set)?;
    let mut out = Vec::new();
    for node in set.iter() {
        // `neighbors` cannot fail: the node was checked above.
        if topology
            .neighbors(node)?
            .iter()
            .any(|&nb| !set.contains(nb))
        {
            out.push(node);
        }
    }
    Ok(NodeSet::new(out))
}

/// Grow a region by `rings` steps: every node within `rings` mesh edges of a
/// member joins. `rings == 0` returns the set unchanged.
pub fn dilate(topology: &SurfaceTopology, set: &NodeSet, rings: usize) -> Result<NodeSet> {
    let n = check_set(topology, set)?;
    let mut inside = set.to_mask(n)?;
    // A frontier-by-frontier search: only the newest nodes need their neighbors
    // looked at in the next round.
    let mut frontier: Vec<u32> = set.iter().collect();
    for _ in 0..rings {
        let mut next = Vec::new();
        for &node in &frontier {
            for &nb in topology.neighbors(node)? {
                if !inside[nb as usize] {
                    inside[nb as usize] = true;
                    next.push(nb);
                }
            }
        }
        if next.is_empty() {
            break;
        }
        frontier = next;
    }
    Ok(NodeSet::from_mask(&inside))
}

/// Shrink a region by `rings` steps: members within `rings` edges of a node that is
/// NOT in the set are removed. Equivalent to removing the boundary `rings` times.
pub fn erode(topology: &SurfaceTopology, set: &NodeSet, rings: usize) -> Result<NodeSet> {
    let n = check_set(topology, set)?;
    let mut inside = set.to_mask(n)?;
    // Start from the current boundary; each round peels it off and the nodes that
    // become exposed are the next boundary.
    let mut frontier: Vec<u32> = boundary_nodes(topology, set)?.iter().collect();
    for _ in 0..rings {
        if frontier.is_empty() {
            break;
        }
        for &node in &frontier {
            inside[node as usize] = false;
        }
        let mut next = Vec::new();
        let mut queued = vec![false; n];
        for &node in &frontier {
            for &nb in topology.neighbors(node)? {
                if inside[nb as usize] && !queued[nb as usize] {
                    queued[nb as usize] = true;
                    next.push(nb);
                }
            }
        }
        frontier = next;
    }
    Ok(NodeSet::from_mask(&inside))
}

// ---------------------------------------------------------------------------
// Grow and shrink by distance along the surface
// ---------------------------------------------------------------------------

/// One entry of the priority queue: smallest distance first, ties by node index.
#[derive(PartialEq)]
struct Reach {
    distance: f64,
    node: u32,
}

impl Eq for Reach {}

impl Ord for Reach {
    fn cmp(&self, other: &Self) -> Ordering {
        // `BinaryHeap` pops the LARGEST item, so reverse the comparison.
        other
            .distance
            .total_cmp(&self.distance)
            .then_with(|| other.node.cmp(&self.node))
    }
}

impl PartialOrd for Reach {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Length of the mesh edge between two nodes.
fn edge_length(mesh: &SurfaceMesh, a: u32, b: u32) -> f64 {
    let (p, q) = (mesh.vertices()[a as usize], mesh.vertices()[b as usize]);
    let d: [f64; 3] = std::array::from_fn(|k| f64::from(p[k]) - f64::from(q[k]));
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// Shortest surface distance from the nearest of `sources` to every node it can
/// reach within `limit` (`None` = anywhere). `distance[i]` is `INFINITY` for a node
/// not reached.
fn distances_from(
    mesh: &SurfaceMesh,
    sources: impl IntoIterator<Item = u32>,
    limit: Option<f64>,
) -> Result<Vec<f64>> {
    let topology = mesh.topology();
    let n = topology.node_count();
    let mut dist = vec![f64::INFINITY; n];
    let mut heap = BinaryHeap::new();
    for s in sources {
        dist[s as usize] = 0.0;
        heap.push(Reach {
            distance: 0.0,
            node: s,
        });
    }
    while let Some(Reach { distance, node }) = heap.pop() {
        if distance > dist[node as usize] {
            continue; // an older, longer entry for this node
        }
        for &nb in topology.neighbors(node)? {
            let d = distance + edge_length(mesh, node, nb);
            if limit.is_some_and(|l| d > l) {
                continue;
            }
            if d < dist[nb as usize] {
                dist[nb as usize] = d;
                heap.push(Reach {
                    distance: d,
                    node: nb,
                });
            }
        }
    }
    Ok(dist)
}

fn check_distance(distance: f64) -> Result<()> {
    ensure_finite("distance", distance)?;
    if distance < 0.0 {
        return Err(Error::InvalidParameter {
            name: "distance".into(),
            reason: format!("{distance} is negative"),
        });
    }
    Ok(())
}

/// Grow a region to every node within `distance` (shortest path along edges, in
/// the mesh's units) of a member.
pub fn dilate_by_distance(mesh: &SurfaceMesh, set: &NodeSet, distance: f64) -> Result<NodeSet> {
    check_distance(distance)?;
    check_set(mesh.topology(), set)?;
    let dist = distances_from(mesh, set.iter(), Some(distance))?;
    Ok(NodeSet::from_iter(
        dist.iter()
            .enumerate()
            .filter(|(_, d)| d.is_finite())
            .map(|(i, _)| i as u32),
    ))
}

/// Shrink a region: keep only members that are MORE than `distance` (along the
/// surface) from every node outside the set.
pub fn erode_by_distance(mesh: &SurfaceMesh, set: &NodeSet, distance: f64) -> Result<NodeSet> {
    check_distance(distance)?;
    let n = check_set(mesh.topology(), set)?;
    let inside = set.to_mask(n)?;
    let outside = (0..n as u32).filter(|&i| !inside[i as usize]);
    let dist = distances_from(mesh, outside, Some(distance))?;
    // A member the search reached is within `distance` of an outsider: remove it.
    Ok(NodeSet::from_iter(
        set.iter().filter(|&node| !dist[node as usize].is_finite()),
    ))
}

// ---------------------------------------------------------------------------
// Connected components
// ---------------------------------------------------------------------------

/// Split a set into its connected pieces (nodes joined by mesh edges that are both
/// in the set). Largest piece first; equal sizes ordered by their smallest node.
pub fn connected_components(topology: &SurfaceTopology, set: &NodeSet) -> Result<Vec<NodeSet>> {
    let n = check_set(topology, set)?;
    let inside = set.to_mask(n)?;
    let mut seen = vec![false; n];
    let mut pieces: Vec<NodeSet> = Vec::new();
    for start in set.iter() {
        if seen[start as usize] {
            continue;
        }
        seen[start as usize] = true;
        let mut queue = VecDeque::from([start]);
        let mut members = vec![start];
        while let Some(node) = queue.pop_front() {
            for &nb in topology.neighbors(node)? {
                if inside[nb as usize] && !seen[nb as usize] {
                    seen[nb as usize] = true;
                    members.push(nb);
                    queue.push_back(nb);
                }
            }
        }
        pieces.push(NodeSet::new(members));
    }
    // Pieces were found in order of their smallest node, so a stable sort by size
    // keeps that order among equals.
    pieces.sort_by_key(|p| std::cmp::Reverse(p.len()));
    Ok(pieces)
}

/// Keep only the largest connected piece (ties: the one with the smallest node).
/// An empty set stays empty.
pub fn keep_largest_component(topology: &SurfaceTopology, set: &NodeSet) -> Result<NodeSet> {
    Ok(connected_components(topology, set)?
        .into_iter()
        .next()
        .unwrap_or_default())
}

/// Drop connected pieces with fewer than `min_nodes` nodes.
pub fn remove_small_components(
    topology: &SurfaceTopology,
    set: &NodeSet,
    min_nodes: usize,
) -> Result<NodeSet> {
    Ok(connected_components(topology, set)?
        .into_iter()
        .filter(|p| p.len() >= min_nodes)
        .fold(NodeSet::empty(), |acc, p| acc.union(&p)))
}

// ---------------------------------------------------------------------------
// Paths
// ---------------------------------------------------------------------------

/// The shortest path along mesh edges from `from` to `to`, both included, in
/// order. `from == to` gives the single node. Errors if the nodes are in different
/// connected pieces of the surface. Equal-length alternatives are resolved toward
/// the smaller node index, so the same call always returns the same path.
pub fn shortest_path(mesh: &SurfaceMesh, from: u32, to: u32) -> Result<Vec<u32>> {
    let topology = mesh.topology();
    let n = topology.node_count();
    for node in [from, to] {
        if node as usize >= n {
            return Err(Error::IndexOutOfRange {
                index: i64::from(node),
                len: n,
            });
        }
    }
    let mut dist = vec![f64::INFINITY; n];
    let mut parent: Vec<Option<u32>> = vec![None; n];
    let mut heap = BinaryHeap::new();
    dist[from as usize] = 0.0;
    heap.push(Reach {
        distance: 0.0,
        node: from,
    });
    while let Some(Reach { distance, node }) = heap.pop() {
        if node == to {
            break;
        }
        if distance > dist[node as usize] {
            continue;
        }
        for &nb in topology.neighbors(node)? {
            let d = distance + edge_length(mesh, node, nb);
            // Strictly shorter, or equal with a smaller parent: deterministic ties.
            let better = d < dist[nb as usize]
                || (d == dist[nb as usize] && parent[nb as usize].is_some_and(|p| node < p));
            if better {
                dist[nb as usize] = d;
                parent[nb as usize] = Some(node);
                heap.push(Reach {
                    distance: d,
                    node: nb,
                });
            }
        }
    }
    if !dist[to as usize].is_finite() {
        return Err(Error::InvalidParameter {
            name: "path".into(),
            reason: format!("node {to} cannot be reached from node {from}"),
        });
    }
    let mut path = vec![to];
    while let Some(p) = parent[*path.last().unwrap_or(&to) as usize] {
        path.push(p);
        if p == from {
            break;
        }
    }
    path.reverse();
    Ok(path)
}

/// Length of a path along the mesh edges, in the mesh's units. Consecutive nodes
/// must be neighbors.
pub fn path_length(mesh: &SurfaceMesh, path: &[u32]) -> Result<f64> {
    let report = check_path(mesh.topology(), path)?;
    if let Some(i) = report.first_bad_link {
        return Err(Error::InvalidParameter {
            name: "path".into(),
            reason: format!("nodes {} and {} are not neighbors", path[i], path[i + 1]),
        });
    }
    Ok(path.windows(2).map(|w| edge_length(mesh, w[0], w[1])).sum())
}

/// What [`check_path`] found out about a list of nodes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PathReport {
    /// Index `i` of the first pair `(path[i], path[i + 1])` that is neither the
    /// same node nor joined by an edge, or `None` if every step is a real step.
    pub first_bad_link: Option<usize>,
    /// `true` when the path ends where it began (at least a triangle: 4 entries).
    pub is_closed: bool,
    /// `true` if some node appears more than once other than as the closing node.
    pub has_repeats: bool,
}

impl PathReport {
    /// Every step joins neighbors.
    pub fn is_connected(&self) -> bool {
        self.first_bad_link.is_none()
    }
}

/// Examine a list of nodes as a path: is each step along an edge, does it close on
/// itself, does it revisit nodes? Errors only for a node outside the surface.
pub fn check_path(topology: &SurfaceTopology, path: &[u32]) -> Result<PathReport> {
    let n = topology.node_count();
    for &node in path {
        if node as usize >= n {
            return Err(Error::IndexOutOfRange {
                index: i64::from(node),
                len: n,
            });
        }
    }
    let mut first_bad_link = None;
    for (i, w) in path.windows(2).enumerate() {
        let joined = w[0] == w[1] || topology.neighbors(w[0])?.binary_search(&w[1]).is_ok();
        if !joined {
            first_bad_link = Some(i);
            break;
        }
    }
    let is_closed = path.len() >= 4 && path.first() == path.last();
    // Repeats among the nodes, not counting the closing duplicate.
    let body = if is_closed {
        &path[..path.len() - 1]
    } else {
        path
    };
    let mut sorted = body.to_vec();
    sorted.sort_unstable();
    let has_repeats = sorted.windows(2).any(|w| w[0] == w[1]);
    Ok(PathReport {
        first_bad_link,
        is_closed,
        has_repeats,
    })
}

/// The segment that closes an open path: the shortest path from its last node back
/// to its first, both included (so it starts where the path ends and ends where the
/// path began). This is SUMA's "join ends". Errors for an empty path.
pub fn join_ends(mesh: &SurfaceMesh, path: &[u32]) -> Result<Vec<u32>> {
    match (path.first(), path.last()) {
        (Some(&first), Some(&last)) => shortest_path(mesh, last, first),
        _ => Err(Error::Empty("path".into())),
    }
}

// ---------------------------------------------------------------------------
// Fill
// ---------------------------------------------------------------------------

/// The result of [`fill_enclosed`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fill {
    /// The boundary plus every node reached from the seed without crossing it.
    pub nodes: NodeSet,
    /// `true` if the fill reached a node on the rim of an open surface, meaning the
    /// boundary did not really enclose the seed's side (the fill leaked).
    pub touches_surface_rim: bool,
}

/// Fill the area around `seed` bounded by `boundary`: start at the seed and spread
/// through neighbors that are not in the boundary. The result is the boundary plus
/// the region reached, SUMA's "fill area".
///
/// On a closed surface any loop splits it into two sides, so which side is filled
/// depends on the seed. On an open surface a boundary that does not close lets the
/// fill escape to the rim; `touches_surface_rim` reports that. The seed must not be
/// on the boundary.
pub fn fill_enclosed(topology: &SurfaceTopology, boundary: &NodeSet, seed: u32) -> Result<Fill> {
    let n = check_set(topology, boundary)?;
    if seed as usize >= n {
        return Err(Error::IndexOutOfRange {
            index: i64::from(seed),
            len: n,
        });
    }
    if boundary.contains(seed) {
        return Err(Error::InvalidParameter {
            name: "seed".into(),
            reason: format!("node {seed} is on the boundary; pick a node on one side of it"),
        });
    }
    // Which nodes lie on the rim of the surface (an edge used by one triangle).
    let mut rim = vec![false; n];
    for edge in topology.edges() {
        if edge.faces.len() == 1 {
            rim[edge.nodes[0] as usize] = true;
            rim[edge.nodes[1] as usize] = true;
        }
    }
    let blocked = boundary.to_mask(n)?;
    let mut reached = vec![false; n];
    reached[seed as usize] = true;
    let mut queue = VecDeque::from([seed]);
    let mut touches_surface_rim = false;
    while let Some(node) = queue.pop_front() {
        touches_surface_rim |= rim[node as usize];
        for &nb in topology.neighbors(node)? {
            if !blocked[nb as usize] && !reached[nb as usize] {
                reached[nb as usize] = true;
                queue.push_back(nb);
            }
        }
    }
    Ok(Fill {
        nodes: NodeSet::from_mask(&reached).union(boundary),
        touches_surface_rim,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A (w+1) x (h+1) grid of unit squares in the z = 0 plane, each split into two
    /// triangles. Node (i, j) is `j * (w + 1) + i`.
    fn grid(w: u32, h: u32) -> SurfaceMesh {
        let stride = w + 1;
        let mut vertices = Vec::new();
        for j in 0..=h {
            for i in 0..=w {
                vertices.push([i as f32, j as f32, 0.0]);
            }
        }
        let mut faces = Vec::new();
        for j in 0..h {
            for i in 0..w {
                let a = j * stride + i;
                let (b, c, d) = (a + 1, a + stride, a + stride + 1);
                faces.push([a, b, d]);
                faces.push([a, d, c]);
            }
        }
        SurfaceMesh::from_triangles(vertices, faces).unwrap()
    }

    fn id(w: u32, i: u32, j: u32) -> u32 {
        j * (w + 1) + i
    }

    /// A solid block of grid nodes `[i0, i1] x [j0, j1]`.
    fn block(w: u32, i0: u32, i1: u32, j0: u32, j1: u32) -> NodeSet {
        let mut v = Vec::new();
        for j in j0..=j1 {
            for i in i0..=i1 {
                v.push(id(w, i, j));
            }
        }
        NodeSet::new(v)
    }

    #[test]
    fn boundary_of_a_block_is_its_rim() {
        let m = grid(8, 8);
        let b = block(8, 2, 5, 2, 5); // 4 x 4 nodes
        let edge = boundary_nodes(m.topology(), &b).unwrap();
        // A 4 x 4 block has 12 rim nodes; the inner 2 x 2 are not boundary.
        assert_eq!(edge.len(), 12);
        assert!(!edge.contains(id(8, 3, 3)) && edge.contains(id(8, 2, 2)));
        assert!(boundary_nodes(m.topology(), &NodeSet::new([999])).is_err());
    }

    #[test]
    fn dilate_and_erode_by_rings() {
        let m = grid(10, 10);
        let t = m.topology();
        let single = NodeSet::new([id(10, 5, 5)]);
        // On this grid a node touches 6 neighbors (4 axis + one diagonal each side).
        assert_eq!(dilate(t, &single, 0).unwrap(), single);
        assert_eq!(dilate(t, &single, 1).unwrap().len(), 7);
        let big = dilate(t, &single, 2).unwrap();
        assert!(big.len() > 7);
        // Eroding what was dilated by the same number of rings recovers the seed.
        assert_eq!(
            erode(t, &dilate(t, &single, 2).unwrap(), 2).unwrap(),
            single
        );
        assert!(erode(t, &single, 1).unwrap().is_empty());
        // Many rings stop at the whole surface.
        assert_eq!(dilate(t, &single, 50).unwrap().len(), t.node_count());
        // Erosion by 0 is a no-op; erosion removes exactly the boundary once.
        let b = block(10, 2, 7, 2, 7);
        assert_eq!(erode(t, &b, 0).unwrap(), b);
        assert_eq!(erode(t, &b, 1).unwrap(), block(10, 3, 6, 3, 6));
    }

    #[test]
    fn distance_based_grow_and_shrink() {
        let m = grid(10, 10);
        let single = NodeSet::new([id(10, 5, 5)]);
        // Edge lengths are 1 (axis) and sqrt(2) (diagonal): 1.0 reaches 4 axis
        // neighbors, 1.5 adds the 2 diagonal ones, 0.5 reaches nothing new.
        assert_eq!(dilate_by_distance(&m, &single, 0.5).unwrap(), single);
        assert_eq!(dilate_by_distance(&m, &single, 1.0).unwrap().len(), 5);
        assert_eq!(dilate_by_distance(&m, &single, 1.5).unwrap().len(), 7);
        assert!(dilate_by_distance(&m, &single, -1.0).is_err());
        assert!(dilate_by_distance(&m, &single, f64::NAN).is_err());
        // Shrinking a 6 x 6 block by distance 1.0 removes nodes within 1.0 of an
        // outsider: the rim.
        let b = block(10, 2, 7, 2, 7);
        assert_eq!(
            erode_by_distance(&m, &b, 1.0).unwrap(),
            block(10, 3, 6, 3, 6)
        );
    }

    #[test]
    fn components_are_ordered_and_cleanup_works() {
        let m = grid(10, 10);
        let t = m.topology();
        let a = block(10, 0, 2, 0, 2); // 9 nodes
        let b = block(10, 6, 7, 6, 7); // 4 nodes
        let c = NodeSet::new([id(10, 9, 0)]); // 1 node
        let all = a.union(&b).union(&c);
        let parts = connected_components(t, &all).unwrap();
        assert_eq!(parts, vec![a.clone(), b.clone(), c.clone()]);
        assert_eq!(keep_largest_component(t, &all).unwrap(), a);
        assert_eq!(remove_small_components(t, &all, 4).unwrap(), a.union(&b));
        assert_eq!(
            remove_small_components(t, &all, 10).unwrap(),
            NodeSet::empty()
        );
        assert_eq!(
            keep_largest_component(t, &NodeSet::empty()).unwrap(),
            NodeSet::empty()
        );
    }

    #[test]
    fn shortest_paths_are_deterministic_and_checked() {
        let m = grid(6, 6);
        let from = id(6, 0, 0);
        let to = id(6, 4, 0);
        let p = shortest_path(&m, from, to).unwrap();
        // Straight along the bottom edge: 5 nodes, length 4.
        assert_eq!(p, vec![0, 1, 2, 3, 4]);
        assert!((path_length(&m, &p).unwrap() - 4.0).abs() < 1e-9);
        // Diagonal steps are shorter than the L-shaped route.
        let diag = shortest_path(&m, id(6, 0, 0), id(6, 3, 3)).unwrap();
        assert_eq!(diag.len(), 4);
        assert!((path_length(&m, &diag).unwrap() - 3.0 * 2.0_f64.sqrt()).abs() < 1e-5);
        assert_eq!(shortest_path(&m, 5, 5).unwrap(), vec![5]);
        assert_eq!(shortest_path(&m, from, to).unwrap(), p, "repeatable");
        assert!(shortest_path(&m, 0, 999).is_err());
    }

    #[test]
    fn unreachable_nodes_are_an_error() {
        // Two separate triangles.
        let m = SurfaceMesh::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [5.0, 0.0, 0.0],
                [6.0, 0.0, 0.0],
                [5.0, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [3, 4, 5]],
        )
        .unwrap();
        assert!(shortest_path(&m, 0, 4).is_err());
    }

    #[test]
    fn path_checks_and_joining_the_ends() {
        let m = grid(6, 6);
        let t = m.topology();
        let open = vec![id(6, 1, 1), id(6, 2, 1), id(6, 3, 1), id(6, 3, 2)];
        let r = check_path(t, &open).unwrap();
        assert!(r.is_connected() && !r.is_closed && !r.has_repeats);
        // A gap is found and located.
        let broken = vec![0, 1, 3, 4];
        assert_eq!(check_path(t, &broken).unwrap().first_bad_link, Some(1));
        assert!(path_length(&m, &broken).is_err());
        // Joining the ends adds the way back; the result is closed.
        let back = join_ends(&m, &open).unwrap();
        assert_eq!((back[0], *back.last().unwrap()), (open[3], open[0]));
        let mut closed = open.clone();
        closed.extend_from_slice(&back[1..]);
        let r = check_path(t, &closed).unwrap();
        assert!(r.is_connected() && r.is_closed);
        // The shortest way back here cuts the corner through a node the path already
        // used ((2,1) is a diagonal neighbor of (3,2)), and `has_repeats` says so.
        assert!(r.has_repeats);
        assert!(join_ends(&m, &[]).is_err());
        // A revisited node is flagged.
        assert!(check_path(t, &[1, 2, 1]).unwrap().has_repeats);
        assert!(check_path(t, &[0, 999]).is_err());
    }

    #[test]
    fn fill_stops_at_a_closed_loop_and_reports_leaks() {
        let m = grid(8, 8);
        let t = m.topology();
        // The rim of a block is a closed loop around its inner 2 x 2 nodes.
        let ring = boundary_nodes(t, &block(8, 2, 5, 2, 5)).unwrap();
        let fill = fill_enclosed(t, &ring, id(8, 3, 3)).unwrap();
        assert_eq!(fill.nodes, block(8, 2, 5, 2, 5));
        assert!(!fill.touches_surface_rim);
        // Seeding the other side floods the rest of the surface (and reaches its rim).
        let outside = fill_enclosed(t, &ring, id(8, 0, 0)).unwrap();
        assert!(outside.touches_surface_rim);
        assert_eq!(outside.nodes.len(), t.node_count() - 4);
        // A seed on the boundary is rejected.
        assert!(fill_enclosed(t, &ring, id(8, 2, 2)).is_err());
        // An open line cannot enclose anything on an open surface: the fill leaks.
        let line = block(8, 1, 6, 4, 4);
        assert!(
            fill_enclosed(t, &line, id(8, 3, 3))
                .unwrap()
                .touches_surface_rim
        );
    }
}
