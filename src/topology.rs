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
// Triangle-mesh TOPOLOGY: which nodes are joined to which, with no coordinates.
// A `SurfaceTopology` is a validated list of triangles over `node_count` nodes
// together with everything derived from connectivity alone: node neighbors, the
// faces around each node, edges and the faces on each edge, face neighbors,
// connected components, ring neighborhoods, and a diagnostic report (boundary,
// non-manifold, inconsistent winding, ...).
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `mesh.rs` adds coordinates (`SurfaceMesh`) to compute normals, areas and
//   distances; `cluster.rs` uses both to group suprathreshold nodes.
// * `domain.rs` (Phase 1) says how many nodes a surface has and may carry its
//   id; the topology's own `TopologyId` identifies the connectivity, so two
//   surfaces with the same node count are not mistaken for each other.
// * `afni-io` supplies triangles (`Surface::faces`); an adapter converts them.
//
// DESIGN DECISIONS
//
// * Only out-of-range node indices and an empty node set are ERRORS. Real
//   meshes have holes, repeated vertices and the occasional bow-tie, and a viewer
//   must still be able to open them, so those are reported in `TopologyReport`
//   instead (call `TopologyReport::is_clean` or `require_clean` to be strict).
// * Everything is deterministic: neighbor lists and edge lists are sorted, so the
//   same triangles always give the same answers in the same order.
// ---------------------------------------------------------------------------

//! Validated triangle-mesh topology and connectivity queries.

use std::collections::{BTreeMap, VecDeque};

use crate::error::{Error, Result};
use crate::numeric::usize_to_u32;

/// A stable identity for a mesh's connectivity: a 64-bit FNV-1a hash of the node
/// count and the triangle list, in order. Equal triangles over the same number of
/// nodes always give the same id, on every platform; any change to a triangle or
/// to the order of triangles changes it. (It is an identity for caching and for
/// "is this the same surface?", not a cryptographic hash.)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TopologyId(u64);

impl TopologyId {
    /// The raw 64-bit value.
    pub fn value(self) -> u64 {
        self.0
    }
}

/// One undirected edge and the triangles that contain it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Edge {
    /// The two end nodes, smaller index first.
    pub nodes: [u32; 2],
    /// Indices of the triangles that contain this edge, ascending. One entry means
    /// a boundary edge, two a normal interior edge, more a non-manifold edge.
    pub faces: Vec<u32>,
}

/// Connected components of a mesh's nodes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Components {
    /// Component number of each node. Numbers are assigned in order of the lowest
    /// node index in each component, so they are deterministic.
    pub labels: Vec<u32>,
    /// Number of nodes in each component, indexed by component number.
    pub sizes: Vec<usize>,
}

impl Components {
    /// How many components there are.
    pub fn count(&self) -> usize {
        self.sizes.len()
    }
}

/// Everything unusual about a mesh's connectivity.
///
/// A well-formed closed surface has all lists empty, one component, and an Euler
/// characteristic of 2 (a sphere-like surface).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologyReport {
    /// Nodes that belong to no triangle.
    pub isolated_nodes: Vec<u32>,
    /// Triangles that name the same node twice (zero area by construction).
    pub degenerate_faces: Vec<u32>,
    /// Triangles that repeat an earlier triangle's node set (any order).
    pub duplicate_faces: Vec<u32>,
    /// Edges with exactly one triangle: the surface boundary. Indices into
    /// [`SurfaceTopology::edges`].
    pub boundary_edges: Vec<usize>,
    /// Edges shared by more than two triangles.
    pub non_manifold_edges: Vec<usize>,
    /// Interior edges whose two triangles traverse the edge in the SAME direction,
    /// which means the two triangles are wound inconsistently.
    pub inconsistent_winding_edges: Vec<usize>,
    /// Nodes whose triangles do not form a single fan (a "bow-tie" vertex where two
    /// parts of the surface touch at a point).
    pub non_manifold_nodes: Vec<u32>,
    /// Nodes that touch a boundary edge.
    pub boundary_nodes: Vec<u32>,
    /// `V - E + F` counting only nodes that belong to a triangle.
    pub euler_characteristic: i64,
    /// Number of connected pieces (isolated nodes count as pieces of their own).
    pub component_count: usize,
}

impl TopologyReport {
    /// True when nothing unusual was found: no isolated nodes, degenerate or
    /// duplicate triangles, boundary, non-manifold or inconsistently wound edges,
    /// or bow-tie nodes. (A single closed sphere-like surface satisfies this;
    /// an open patch does not, because of its boundary.)
    pub fn is_clean(&self) -> bool {
        self.isolated_nodes.is_empty()
            && self.degenerate_faces.is_empty()
            && self.duplicate_faces.is_empty()
            && self.boundary_edges.is_empty()
            && self.non_manifold_edges.is_empty()
            && self.inconsistent_winding_edges.is_empty()
            && self.non_manifold_nodes.is_empty()
    }

    /// True when the surface has no boundary and no non-manifold edges.
    pub fn is_closed_manifold(&self) -> bool {
        self.boundary_edges.is_empty() && self.non_manifold_edges.is_empty()
    }
}

/// A validated triangle mesh's connectivity.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceTopology {
    node_count: usize,
    faces: Vec<[u32; 3]>,
    edges: Vec<Edge>,
    /// Sorted neighbor list per node.
    neighbors: Vec<Vec<u32>>,
    /// Ascending list of triangles around each node.
    node_faces: Vec<Vec<u32>>,
    id: TopologyId,
}

/// FNV-1a over a stream of `u32`s (little-endian bytes).
fn fnv1a(values: impl Iterator<Item = u32>) -> u64 {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for v in values {
        for byte in v.to_le_bytes() {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    hash
}

impl SurfaceTopology {
    /// Build the topology of `faces` over `node_count` nodes.
    ///
    /// Fails if `node_count` is zero, if it does not fit a `u32` index, or if any
    /// triangle names a node `>= node_count`. Everything else is reported by
    /// [`report`](Self::report).
    pub fn new(node_count: usize, faces: Vec<[u32; 3]>) -> Result<Self> {
        if node_count == 0 {
            return Err(Error::Empty("surface nodes".into()));
        }
        usize_to_u32(node_count)?; // node indices are u32
        usize_to_u32(faces.len())?;
        for face in &faces {
            for &n in face {
                if n as usize >= node_count {
                    return Err(Error::IndexOutOfRange {
                        index: i64::from(n),
                        len: node_count,
                    });
                }
            }
        }
        // Edges, keyed by sorted end points so the order is deterministic.
        let mut edge_faces: BTreeMap<(u32, u32), Vec<u32>> = BTreeMap::new();
        let mut node_faces: Vec<Vec<u32>> = vec![Vec::new(); node_count];
        for (fi, face) in faces.iter().enumerate() {
            let fi = fi as u32;
            for k in 0..3 {
                let (a, b) = (face[k], face[(k + 1) % 3]);
                if a != b {
                    edge_faces.entry((a.min(b), a.max(b))).or_default().push(fi);
                }
                // A triangle that repeats a node is listed once per node.
                if node_faces[a as usize].last() != Some(&fi) {
                    node_faces[a as usize].push(fi);
                }
            }
        }
        let mut neighbors: Vec<Vec<u32>> = vec![Vec::new(); node_count];
        let mut edges = Vec::with_capacity(edge_faces.len());
        for ((a, b), mut fs) in edge_faces {
            fs.dedup(); // a degenerate triangle can name an edge twice
            neighbors[a as usize].push(b);
            neighbors[b as usize].push(a);
            edges.push(Edge {
                nodes: [a, b],
                faces: fs,
            });
        }
        for list in &mut neighbors {
            list.sort_unstable();
        }
        let id = TopologyId(fnv1a(
            std::iter::once(node_count as u32).chain(faces.iter().flatten().copied()),
        ));
        Ok(Self {
            node_count,
            faces,
            edges,
            neighbors,
            node_faces,
            id,
        })
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.node_count
    }

    /// Number of triangles.
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// The triangles, as given.
    pub fn faces(&self) -> &[[u32; 3]] {
        &self.faces
    }

    /// The undirected edges, sorted by `(smaller node, larger node)`.
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    /// The stable connectivity identity.
    pub fn id(&self) -> TopologyId {
        self.id
    }

    /// The neighbors of `node` (joined to it by an edge), ascending.
    pub fn neighbors(&self, node: u32) -> Result<&[u32]> {
        self.neighbors
            .get(node as usize)
            .map(Vec::as_slice)
            .ok_or(Error::IndexOutOfRange {
                index: i64::from(node),
                len: self.node_count,
            })
    }

    /// All neighbor lists, indexed by node. For algorithms that walk the whole mesh.
    pub fn all_neighbors(&self) -> &[Vec<u32>] {
        &self.neighbors
    }

    /// The triangles that contain `node`, ascending.
    pub fn faces_of_node(&self, node: u32) -> Result<&[u32]> {
        self.node_faces
            .get(node as usize)
            .map(Vec::as_slice)
            .ok_or(Error::IndexOutOfRange {
                index: i64::from(node),
                len: self.node_count,
            })
    }

    /// The triangles that share an edge with triangle `face` (ascending, without
    /// itself). Up to three for a manifold surface.
    pub fn face_neighbors(&self, face: u32) -> Result<Vec<u32>> {
        let f = self
            .faces
            .get(face as usize)
            .ok_or(Error::IndexOutOfRange {
                index: i64::from(face),
                len: self.faces.len(),
            })?;
        let mut out = Vec::new();
        for k in 0..3 {
            let (a, b) = (f[k], f[(k + 1) % 3]);
            if a == b {
                continue;
            }
            if let Ok(i) = self
                .edges
                .binary_search_by_key(&(a.min(b), a.max(b)), |e| (e.nodes[0], e.nodes[1]))
            {
                out.extend(self.edges[i].faces.iter().copied().filter(|&g| g != face));
            }
        }
        out.sort_unstable();
        out.dedup();
        Ok(out)
    }

    /// Connected components of the node graph (isolated nodes are components of
    /// one).
    pub fn connected_components(&self) -> Components {
        let mut labels = vec![u32::MAX; self.node_count];
        let mut sizes = Vec::new();
        for start in 0..self.node_count {
            if labels[start] != u32::MAX {
                continue;
            }
            let id = sizes.len() as u32;
            let mut size = 0;
            let mut queue = VecDeque::from([start as u32]);
            labels[start] = id;
            while let Some(n) = queue.pop_front() {
                size += 1;
                for &m in &self.neighbors[n as usize] {
                    if labels[m as usize] == u32::MAX {
                        labels[m as usize] = id;
                        queue.push_back(m);
                    }
                }
            }
            sizes.push(size);
        }
        Components { labels, sizes }
    }

    /// Breadth-first layers around `seed`: layer 0 is `[seed]`, layer `k` holds the
    /// nodes whose shortest edge path from `seed` has exactly `k` edges, each layer
    /// ascending. Stops after `max_rings` layers beyond the seed (so the result has
    /// at most `max_rings + 1` layers) or when the component is exhausted.
    pub fn ring_layers(&self, seed: u32, max_rings: usize) -> Result<Vec<Vec<u32>>> {
        self.neighbors(seed)?; // range check
        let mut seen = vec![false; self.node_count];
        seen[seed as usize] = true;
        let mut layers = vec![vec![seed]];
        for _ in 0..max_rings {
            let mut next = Vec::new();
            for &n in layers.last().expect("at least the seed layer") {
                for &m in &self.neighbors[n as usize] {
                    if !seen[m as usize] {
                        seen[m as usize] = true;
                        next.push(m);
                    }
                }
            }
            if next.is_empty() {
                break;
            }
            next.sort_unstable();
            layers.push(next);
        }
        Ok(layers)
    }

    /// Every node within `rings` edges of `seed` (including `seed`), ascending.
    pub fn within_rings(&self, seed: u32, rings: usize) -> Result<Vec<u32>> {
        let mut all: Vec<u32> = self
            .ring_layers(seed, rings)?
            .into_iter()
            .flatten()
            .collect();
        all.sort_unstable();
        Ok(all)
    }

    /// Diagnose the connectivity. Cost is linear in the mesh size.
    pub fn report(&self) -> TopologyReport {
        let isolated_nodes: Vec<u32> = (0..self.node_count as u32)
            .filter(|&n| self.node_faces[n as usize].is_empty())
            .collect();

        let degenerate_faces: Vec<u32> = self
            .faces
            .iter()
            .enumerate()
            .filter(|(_, f)| f[0] == f[1] || f[1] == f[2] || f[0] == f[2])
            .map(|(i, _)| i as u32)
            .collect();

        // A triangle repeats an earlier one if its sorted node triple was seen.
        let mut seen: BTreeMap<[u32; 3], u32> = BTreeMap::new();
        let mut duplicate_faces = Vec::new();
        for (i, f) in self.faces.iter().enumerate() {
            let mut key = *f;
            key.sort_unstable();
            if seen.insert(key, i as u32).is_some() {
                duplicate_faces.push(i as u32);
            }
        }

        let mut boundary_edges = Vec::new();
        let mut non_manifold_edges = Vec::new();
        let mut inconsistent_winding_edges = Vec::new();
        let mut boundary_flag = vec![false; self.node_count];
        for (ei, edge) in self.edges.iter().enumerate() {
            match edge.faces.len() {
                1 => {
                    boundary_edges.push(ei);
                    boundary_flag[edge.nodes[0] as usize] = true;
                    boundary_flag[edge.nodes[1] as usize] = true;
                }
                2 => {
                    // Consistent winding traverses a shared edge in OPPOSITE directions.
                    let dir = |face: u32| -> Option<bool> {
                        let f = self.faces[face as usize];
                        (0..3).find_map(|k| {
                            let (a, b) = (f[k], f[(k + 1) % 3]);
                            (a.min(b) == edge.nodes[0] && a.max(b) == edge.nodes[1])
                                .then_some(a < b)
                        })
                    };
                    if dir(edge.faces[0]) == dir(edge.faces[1]) {
                        inconsistent_winding_edges.push(ei);
                    }
                }
                n if n > 2 => non_manifold_edges.push(ei),
                _ => {}
            }
        }
        let boundary_nodes: Vec<u32> = (0..self.node_count as u32)
            .filter(|&n| boundary_flag[n as usize])
            .collect();

        // A node is a bow-tie if the triangles around it fall into more than one
        // group when two triangles are grouped by sharing an edge THROUGH that node.
        let mut non_manifold_nodes = Vec::new();
        for n in 0..self.node_count as u32 {
            let around = &self.node_faces[n as usize];
            if around.len() < 2 {
                continue;
            }
            let mut group: Vec<usize> = (0..around.len()).collect();
            let find = |g: &mut Vec<usize>, mut x: usize| {
                while g[x] != x {
                    g[x] = g[g[x]];
                    x = g[x];
                }
                x
            };
            for (i, &fi) in around.iter().enumerate() {
                for (j, &fj) in around.iter().enumerate().skip(i + 1) {
                    // Do these two triangles share an edge that contains `n`?
                    let (a, b) = (self.faces[fi as usize], self.faces[fj as usize]);
                    let shared = a.iter().filter(|x| **x != n && b.contains(x)).count() >= 1;
                    if shared {
                        let (ri, rj) = (find(&mut group, i), find(&mut group, j));
                        group[ri] = rj;
                    }
                }
            }
            let roots: std::collections::BTreeSet<usize> =
                (0..around.len()).map(|i| find(&mut group, i)).collect();
            if roots.len() > 1 {
                non_manifold_nodes.push(n);
            }
        }

        let used_nodes = self.node_count - isolated_nodes.len();
        let euler_characteristic =
            used_nodes as i64 - self.edges.len() as i64 + self.faces.len() as i64;
        TopologyReport {
            isolated_nodes,
            degenerate_faces,
            duplicate_faces,
            boundary_edges,
            non_manifold_edges,
            inconsistent_winding_edges,
            non_manifold_nodes,
            boundary_nodes,
            euler_characteristic,
            component_count: self.connected_components().count(),
        }
    }

    /// Fail unless [`report`](Self::report) is [clean](TopologyReport::is_clean).
    /// The error names the first problem found.
    pub fn require_clean(&self) -> Result<()> {
        let r = self.report();
        let problem = if !r.isolated_nodes.is_empty() {
            Some(format!("{} isolated node(s)", r.isolated_nodes.len()))
        } else if !r.degenerate_faces.is_empty() {
            Some(format!(
                "{} degenerate triangle(s)",
                r.degenerate_faces.len()
            ))
        } else if !r.duplicate_faces.is_empty() {
            Some(format!("{} duplicate triangle(s)", r.duplicate_faces.len()))
        } else if !r.non_manifold_edges.is_empty() {
            Some(format!(
                "{} non-manifold edge(s)",
                r.non_manifold_edges.len()
            ))
        } else if !r.non_manifold_nodes.is_empty() {
            Some(format!(
                "{} non-manifold node(s)",
                r.non_manifold_nodes.len()
            ))
        } else if !r.inconsistent_winding_edges.is_empty() {
            Some(format!(
                "{} inconsistently wound edge(s)",
                r.inconsistent_winding_edges.len()
            ))
        } else if !r.boundary_edges.is_empty() {
            Some(format!("{} boundary edge(s)", r.boundary_edges.len()))
        } else {
            None
        };
        match problem {
            None => Ok(()),
            Some(p) => Err(Error::InvalidParameter {
                name: "surface topology".into(),
                reason: p,
            }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tetrahedron: 4 nodes, 4 triangles, closed, wound outward (consistently).
    fn tetra() -> SurfaceTopology {
        SurfaceTopology::new(4, vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]]).unwrap()
    }

    /// A flat square split into two triangles: an open patch.
    fn square() -> SurfaceTopology {
        SurfaceTopology::new(4, vec![[0, 1, 2], [0, 2, 3]]).unwrap()
    }

    #[test]
    fn construction_validates_indices_and_counts() {
        assert!(SurfaceTopology::new(0, vec![]).is_err());
        assert!(matches!(
            SurfaceTopology::new(3, vec![[0, 1, 3]]),
            Err(Error::IndexOutOfRange { index: 3, len: 3 })
        ));
        assert!(
            SurfaceTopology::new(3, vec![]).is_ok(),
            "no triangles is allowed (all isolated)"
        );
    }

    #[test]
    fn tetrahedron_is_a_clean_closed_sphere() {
        let t = tetra();
        assert_eq!((t.node_count(), t.face_count(), t.edges().len()), (4, 4, 6));
        let r = t.report();
        assert!(r.is_clean() && r.is_closed_manifold(), "{r:?}");
        assert_eq!(r.euler_characteristic, 2);
        assert_eq!(r.component_count, 1);
        assert!(t.require_clean().is_ok());
        for n in 0..4 {
            assert_eq!(t.neighbors(n).unwrap().len(), 3);
            assert_eq!(t.faces_of_node(n).unwrap().len(), 3);
        }
        for f in 0..4 {
            assert_eq!(t.face_neighbors(f).unwrap().len(), 3);
        }
    }

    #[test]
    fn open_patch_reports_its_boundary() {
        let t = square();
        let r = t.report();
        assert_eq!(r.boundary_edges.len(), 4);
        assert_eq!(r.boundary_nodes, vec![0, 1, 2, 3]);
        assert!(!r.is_closed_manifold() && !r.is_clean());
        assert_eq!(r.euler_characteristic, 1, "a disk");
        assert!(t
            .require_clean()
            .unwrap_err()
            .to_string()
            .contains("boundary"));
        assert_eq!(t.face_neighbors(0).unwrap(), vec![1]);
        // Edge 0-2 is the shared diagonal; its two triangles are consistently wound.
        assert!(r.inconsistent_winding_edges.is_empty());
    }

    #[test]
    fn neighbor_and_edge_lists_are_sorted_and_deterministic() {
        let t = SurfaceTopology::new(5, vec![[4, 2, 0], [0, 2, 1], [3, 4, 0]]).unwrap();
        for n in 0..5 {
            let list = t.neighbors(n).unwrap();
            assert!(list.windows(2).all(|w| w[0] < w[1]), "{list:?}");
        }
        let nodes: Vec<[u32; 2]> = t.edges().iter().map(|e| e.nodes).collect();
        let mut sorted = nodes.clone();
        sorted.sort();
        assert_eq!(nodes, sorted);
        assert!(t.edges().iter().all(|e| e.nodes[0] < e.nodes[1]));
        assert!(t.neighbors(5).is_err());
    }

    #[test]
    fn diagnostics_find_every_kind_of_defect() {
        // Node 5 is isolated; triangle 1 is degenerate; triangle 2 duplicates 0.
        let t = SurfaceTopology::new(6, vec![[0, 1, 2], [3, 3, 4], [2, 1, 0]]).unwrap();
        let r = t.report();
        assert_eq!(r.isolated_nodes, vec![5]);
        assert_eq!(r.degenerate_faces, vec![1]);
        assert_eq!(r.duplicate_faces, vec![2]);
        // Triangles 0 and 2 traverse their shared edges in opposite directions, so
        // that is consistent; but each of those edges now has TWO faces (fine).
        assert!(r.inconsistent_winding_edges.is_empty());
        assert!(t.require_clean().is_err());
    }

    #[test]
    fn inconsistent_winding_and_non_manifold_edges_are_detected() {
        // Two triangles sharing edge 0-1 but wound the same way along it.
        let flipped = SurfaceTopology::new(4, vec![[0, 1, 2], [0, 1, 3]]).unwrap();
        assert_eq!(flipped.report().inconsistent_winding_edges.len(), 1);
        // Three triangles on one edge.
        let fin = SurfaceTopology::new(5, vec![[0, 1, 2], [0, 1, 3], [0, 1, 4]]).unwrap();
        let r = fin.report();
        assert_eq!(r.non_manifold_edges.len(), 1);
        assert!(!r.is_closed_manifold());
    }

    #[test]
    fn bow_tie_nodes_are_found() {
        // Two triangles that touch only at node 2.
        let t = SurfaceTopology::new(5, vec![[0, 1, 2], [2, 3, 4]]).unwrap();
        let r = t.report();
        assert_eq!(r.non_manifold_nodes, vec![2]);
        assert_eq!(
            r.component_count, 1,
            "they are still connected through node 2"
        );
        // A fan around a node is not a bow-tie.
        let fan = SurfaceTopology::new(5, vec![[0, 1, 2], [0, 2, 3], [0, 3, 4]]).unwrap();
        assert!(fan.report().non_manifold_nodes.is_empty());
    }

    #[test]
    fn components_are_numbered_by_lowest_node() {
        let t = SurfaceTopology::new(7, vec![[0, 1, 2], [4, 5, 6]]).unwrap();
        let c = t.connected_components();
        assert_eq!(c.count(), 3); // {0,1,2}, {3}, {4,5,6}
        assert_eq!(c.labels, vec![0, 0, 0, 1, 2, 2, 2]);
        assert_eq!(c.sizes, vec![3, 1, 3]);
    }

    #[test]
    fn ring_layers_are_breadth_first_and_bounded() {
        // A path-like strip: 0-1-2-3-4 along triangles.
        let t = SurfaceTopology::new(5, vec![[0, 1, 2], [1, 2, 3], [2, 3, 4]]).unwrap();
        let layers = t.ring_layers(0, 10).unwrap();
        assert_eq!(layers, vec![vec![0], vec![1, 2], vec![3, 4]]);
        assert_eq!(t.ring_layers(0, 1).unwrap(), vec![vec![0], vec![1, 2]]);
        assert_eq!(t.ring_layers(0, 0).unwrap(), vec![vec![0]]);
        assert_eq!(t.within_rings(4, 1).unwrap(), vec![2, 3, 4]);
        assert!(t.ring_layers(9, 1).is_err());
    }

    #[test]
    fn identity_depends_on_connectivity_only() {
        let a = tetra();
        let b = tetra();
        assert_eq!(a.id(), b.id());
        let c = SurfaceTopology::new(4, vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 2, 3]]).unwrap();
        assert_ne!(a.id(), c.id(), "one triangle's winding differs");
        let d = SurfaceTopology::new(5, a.faces().to_vec()).unwrap();
        assert_ne!(
            a.id(),
            d.id(),
            "same triangles over more nodes is a different surface"
        );
        // Stable value: pinned so an accidental change of the hash is noticed.
        assert_eq!(tetra().id().value(), 0x5f91a5095f5a861_u64);
    }
}
