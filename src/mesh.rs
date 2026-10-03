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
// `SurfaceMesh`: a triangle mesh with coordinates, i.e. a `SurfaceTopology` plus
// one `[x, y, z]` per node. It supplies the geometry that needs coordinates:
// triangle normals, node normals, triangle and node areas, bounds, edge lengths,
// enclosed volume, and distance searches along the mesh's edges (the "graph
// distance" AFNI/SUMA use as an approximation of geodesic distance).
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `topology.rs` provides the connectivity this builds on.
// * `cluster.rs` calls `NeighborhoodSearcher` and `node_areas` to grow and
//   measure clusters.
// * `afni-io` supplies vertices and triangles; coordinates are `f32` (as stored in
//   files and uploaded to GPUs) but ALL arithmetic here is `f64`, per the crate's
//   numeric convention.
//
// CONVENTIONS (checked against SUMA's `SurfaceMetrics`/`SurfMeasures`, see
// tests/mesh_conformance.rs)
//
// * Triangle normal = (v1 - v0) x (v2 - v0), normalized: counter-clockwise
//   triangles, seen from outside, point outward.
// * NODE normal = the normalized SUM OF UNIT triangle normals around the node
//   (each adjacent triangle counts once, regardless of its size or corner angle).
//   Area-weighted and angle-weighted normals are different and were checked NOT to
//   be SUMA's choice.
// * NODE area = one third of the total area of the triangles that contain the node.
// * Degenerate (zero-area) triangles have a zero normal and contribute zero area.
// ---------------------------------------------------------------------------

//! Mesh geometry (normals, areas, bounds, volume) and graph-distance searches.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::error::{Error, Result};
use crate::numeric::ensure_finite;
use crate::topology::SurfaceTopology;

/// An axis-aligned bounding box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds {
    /// Smallest coordinate on each axis.
    pub min: [f64; 3],
    /// Largest coordinate on each axis.
    pub max: [f64; 3],
}

impl Bounds {
    /// The midpoint of the box.
    pub fn center(&self) -> [f64; 3] {
        [0, 1, 2].map(|k| 0.5 * (self.min[k] + self.max[k]))
    }

    /// The box's size on each axis.
    pub fn extent(&self) -> [f64; 3] {
        [0, 1, 2].map(|k| self.max[k] - self.min[k])
    }
}

/// Which way a closed surface's triangle winding faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    /// Triangle normals point away from the enclosed volume (positive volume).
    Outward,
    /// Triangle normals point into the enclosed volume (negative volume).
    Inward,
    /// The mesh is open or inconsistently wound, so "outward" is not defined.
    Undefined,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn widen(v: [f32; 3]) -> [f64; 3] {
    [f64::from(v[0]), f64::from(v[1]), f64::from(v[2])]
}

/// A normalized copy of `v`, or the zero vector if `v` has no length.
fn unit_or_zero(v: [f64; 3]) -> [f64; 3] {
    let n = length(v);
    if n > 0.0 && n.is_finite() {
        [v[0] / n, v[1] / n, v[2] / n]
    } else {
        [0.0; 3]
    }
}

/// A triangle mesh: connectivity plus finite coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceMesh {
    topology: SurfaceTopology,
    vertices: Vec<[f32; 3]>,
    /// Edge length to each neighbor, parallel to `topology.all_neighbors()`.
    neighbor_lengths: Vec<Vec<f64>>,
}

impl SurfaceMesh {
    /// Combine `topology` with one coordinate per node. Every coordinate must be
    /// finite and there must be exactly `topology.node_count()` of them.
    pub fn new(topology: SurfaceTopology, vertices: Vec<[f32; 3]>) -> Result<Self> {
        if vertices.len() != topology.node_count() {
            return Err(Error::LengthMismatch {
                what: "vertex coordinates".into(),
                expected: topology.node_count(),
                found: vertices.len(),
            });
        }
        for (i, v) in vertices.iter().enumerate() {
            for &c in v {
                ensure_finite(&format!("coordinate of node {i}"), f64::from(c))?;
            }
        }
        let neighbor_lengths = topology
            .all_neighbors()
            .iter()
            .enumerate()
            .map(|(n, list)| {
                list.iter()
                    .map(|&m| length(sub(widen(vertices[m as usize]), widen(vertices[n]))))
                    .collect()
            })
            .collect();
        Ok(Self {
            topology,
            vertices,
            neighbor_lengths,
        })
    }

    /// Build the topology from `faces` and attach `vertices`.
    pub fn from_triangles(vertices: Vec<[f32; 3]>, faces: Vec<[u32; 3]>) -> Result<Self> {
        let topology = SurfaceTopology::new(vertices.len(), faces)?;
        Self::new(topology, vertices)
    }

    /// The connectivity.
    pub fn topology(&self) -> &SurfaceTopology {
        &self.topology
    }

    /// The node coordinates.
    pub fn vertices(&self) -> &[[f32; 3]] {
        &self.vertices
    }

    /// The (unnormalized) normal of a triangle: its cross product, whose length is
    /// twice the triangle's area.
    fn raw_face_normal(&self, f: &[u32; 3]) -> [f64; 3] {
        let (a, b, c) = (
            widen(self.vertices[f[0] as usize]),
            widen(self.vertices[f[1] as usize]),
            widen(self.vertices[f[2] as usize]),
        );
        cross(sub(b, a), sub(c, a))
    }

    /// Unit normal of every triangle (the zero vector for a zero-area triangle).
    pub fn face_normals(&self) -> Vec<[f64; 3]> {
        self.topology
            .faces()
            .iter()
            .map(|f| unit_or_zero(self.raw_face_normal(f)))
            .collect()
    }

    /// Area of every triangle.
    pub fn face_areas(&self) -> Vec<f64> {
        self.topology
            .faces()
            .iter()
            .map(|f| 0.5 * length(self.raw_face_normal(f)))
            .collect()
    }

    /// Total surface area.
    pub fn total_area(&self) -> f64 {
        self.face_areas().iter().sum()
    }

    /// The area associated with each node: a third of the area of the triangles
    /// around it. These sum to [`total_area`](Self::total_area) (for every node
    /// that belongs to a triangle).
    pub fn node_areas(&self) -> Vec<f64> {
        let face_areas = self.face_areas();
        let mut areas = vec![0.0; self.vertices.len()];
        for (f, &area) in self.topology.faces().iter().zip(&face_areas) {
            // A triangle naming a node twice counts that node for each corner it
            // occupies; such triangles have zero area anyway.
            for &n in f {
                areas[n as usize] += area / 3.0;
            }
        }
        areas
    }

    /// Unit normal at each node: the normalized sum of the UNIT normals of its
    /// triangles. A node in no triangle (or whose triangle normals cancel exactly)
    /// gets the zero vector.
    pub fn vertex_normals(&self) -> Vec<[f64; 3]> {
        let face_normals = self.face_normals();
        let mut sums = vec![[0.0_f64; 3]; self.vertices.len()];
        for (f, n) in self.topology.faces().iter().zip(&face_normals) {
            // A triangle contributes once to each distinct node it contains.
            let mut seen = [u32::MAX; 3];
            for (k, &node) in f.iter().enumerate() {
                if seen[..k].contains(&node) {
                    continue;
                }
                seen[k] = node;
                for axis in 0..3 {
                    sums[node as usize][axis] += n[axis];
                }
            }
        }
        sums.into_iter().map(unit_or_zero).collect()
    }

    /// The box around all nodes.
    pub fn bounds(&self) -> Bounds {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for v in &self.vertices {
            let w = widen(*v);
            for k in 0..3 {
                min[k] = min[k].min(w[k]);
                max[k] = max[k].max(w[k]);
            }
        }
        Bounds { min, max }
    }

    /// The length of every edge, in the order of
    /// [`SurfaceTopology::edges`].
    pub fn edge_lengths(&self) -> Vec<f64> {
        self.topology
            .edges()
            .iter()
            .map(|e| {
                length(sub(
                    widen(self.vertices[e.nodes[0] as usize]),
                    widen(self.vertices[e.nodes[1] as usize]),
                ))
            })
            .collect()
    }

    /// The enclosed volume from the divergence theorem, positive if the triangles
    /// are wound outward. Meaningful only for a closed, consistently wound surface;
    /// see [`orientation`](Self::orientation).
    pub fn signed_volume(&self) -> f64 {
        self.topology
            .faces()
            .iter()
            .map(|f| {
                let (a, b, c) = (
                    widen(self.vertices[f[0] as usize]),
                    widen(self.vertices[f[1] as usize]),
                    widen(self.vertices[f[2] as usize]),
                );
                dot(a, cross(b, c)) / 6.0
            })
            .sum()
    }

    /// Whether a closed surface's winding faces outward or inward.
    pub fn orientation(&self) -> Orientation {
        let r = self.topology.report();
        if !r.is_closed_manifold()
            || !r.inconsistent_winding_edges.is_empty()
            || self.topology.face_count() == 0
        {
            return Orientation::Undefined;
        }
        match self.signed_volume() {
            v if v > 0.0 => Orientation::Outward,
            v if v < 0.0 => Orientation::Inward,
            _ => Orientation::Undefined,
        }
    }

    /// A reusable searcher for nodes within a distance of a seed along mesh edges.
    pub fn searcher(&self) -> NeighborhoodSearcher<'_> {
        NeighborhoodSearcher::new(self)
    }
}

/// A heap entry ordered so the smallest distance pops first.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Frontier {
    distance: f64,
    node: u32,
}

impl Eq for Frontier {}

impl Ord for Frontier {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reversed so `BinaryHeap` (a max-heap) yields the NEAREST node first; ties
        // pop the lower node number first for determinism.
        other
            .distance
            .total_cmp(&self.distance)
            .then_with(|| other.node.cmp(&self.node))
    }
}

impl PartialOrd for Frontier {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Finds the nodes within a given graph distance of a seed, reusing its buffers
/// between calls so that thousands of searches over a large mesh stay cheap (a
/// fresh `O(nodes)` allocation per search would make clustering quadratic).
#[derive(Debug)]
pub struct NeighborhoodSearcher<'a> {
    mesh: &'a SurfaceMesh,
    distance: Vec<f64>,
    touched: Vec<u32>,
    heap: BinaryHeap<Frontier>,
}

impl<'a> NeighborhoodSearcher<'a> {
    fn new(mesh: &'a SurfaceMesh) -> Self {
        Self {
            mesh,
            distance: vec![f64::INFINITY; mesh.vertices.len()],
            touched: Vec::new(),
            heap: BinaryHeap::new(),
        }
    }

    /// Every node whose SHORTEST PATH along mesh edges from `seed` is at most
    /// `radius`, with that distance, in increasing order of distance (ties by node
    /// number). Includes `seed` itself at distance 0.
    ///
    /// This is a true shortest-path (Dijkstra) distance over the edge graph, the
    /// usual reading of "graph distance". SUMA's own routine is a layered
    /// approximation of it; see the roadmap discovery log for how they differ.
    pub fn within_distance(&mut self, seed: u32, radius: f64) -> Result<Vec<(u32, f64)>> {
        if seed as usize >= self.distance.len() {
            return Err(Error::IndexOutOfRange {
                index: i64::from(seed),
                len: self.distance.len(),
            });
        }
        ensure_finite("search radius", radius)?;
        let mut found = Vec::new();
        if radius < 0.0 {
            return Ok(found);
        }
        self.distance[seed as usize] = 0.0;
        self.touched.push(seed);
        self.heap.push(Frontier {
            distance: 0.0,
            node: seed,
        });
        while let Some(Frontier { distance, node }) = self.heap.pop() {
            if distance > self.distance[node as usize] {
                continue; // a shorter path to this node was already settled
            }
            found.push((node, distance));
            let neighbors = &self.mesh.topology.all_neighbors()[node as usize];
            let lengths = &self.mesh.neighbor_lengths[node as usize];
            for (&next, &edge) in neighbors.iter().zip(lengths) {
                let candidate = distance + edge;
                if candidate <= radius && candidate < self.distance[next as usize] {
                    if self.distance[next as usize].is_infinite() {
                        self.touched.push(next);
                    }
                    self.distance[next as usize] = candidate;
                    self.heap.push(Frontier {
                        distance: candidate,
                        node: next,
                    });
                }
            }
        }
        // Reset only what this search touched.
        for &n in &self.touched {
            self.distance[n as usize] = f64::INFINITY;
        }
        self.touched.clear();
        Ok(found)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A regular tetrahedron with edge length 2*sqrt(2), wound outward.
    fn tetra() -> SurfaceMesh {
        SurfaceMesh::from_triangles(
            vec![
                [1.0, 1.0, 1.0],
                [1.0, -1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
            ],
            vec![[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]],
        )
        .unwrap()
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn construction_checks_vertices() {
        let t = SurfaceTopology::new(3, vec![[0, 1, 2]]).unwrap();
        assert!(SurfaceMesh::new(t.clone(), vec![[0.0; 3]; 2]).is_err());
        assert!(SurfaceMesh::new(
            t.clone(),
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, f32::NAN, 0.0]]
        )
        .is_err());
        assert!(
            SurfaceMesh::new(t, vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]).is_ok()
        );
    }

    #[test]
    fn areas_normals_and_lengths_of_a_known_solid() {
        let m = tetra();
        let edge = 2.0 * 2.0_f64.sqrt();
        let tri = 3.0_f64.sqrt() / 4.0 * edge * edge;
        for a in m.face_areas() {
            assert!(close(a, tri), "{a}");
        }
        assert!(close(m.total_area(), 4.0 * tri));
        assert!(m.edge_lengths().iter().all(|&l| close(l, edge)));
        // Node areas: each node touches three triangles, a third of each.
        for a in m.node_areas() {
            assert!(close(a, tri));
        }
        assert!(close(m.node_areas().iter().sum::<f64>(), m.total_area()));
        // Outward normals: point away from the centroid (the origin here).
        for (f, n) in m.topology().faces().iter().zip(m.face_normals()) {
            let c = f.iter().fold([0.0; 3], |s, &i| {
                let v = widen(m.vertices()[i as usize]);
                [s[0] + v[0] / 3.0, s[1] + v[1] / 3.0, s[2] + v[2] / 3.0]
            });
            assert!(dot(n, c) > 0.0, "normal {n:?} points inward at {c:?}");
            assert!(close(length(n), 1.0));
        }
        for (i, n) in m.vertex_normals().iter().enumerate() {
            assert!(close(length(*n), 1.0));
            assert!(dot(*n, widen(m.vertices()[i])) > 0.0);
        }
    }

    #[test]
    fn volume_and_orientation() {
        let m = tetra();
        // Regular tetrahedron with edge a has volume a^3 / (6 sqrt 2); a = 2 sqrt 2 gives 8/3.
        assert!(close(m.signed_volume(), 8.0 / 3.0), "{}", m.signed_volume());
        assert_eq!(m.orientation(), Orientation::Outward);
        let flipped = SurfaceMesh::from_triangles(
            m.vertices().to_vec(),
            m.topology()
                .faces()
                .iter()
                .map(|f| [f[0], f[2], f[1]])
                .collect(),
        )
        .unwrap();
        assert!(close(flipped.signed_volume(), -8.0 / 3.0));
        assert_eq!(flipped.orientation(), Orientation::Inward);
        // An open patch has no orientation.
        let patch = SurfaceMesh::from_triangles(
            vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap();
        assert_eq!(patch.orientation(), Orientation::Undefined);
    }

    #[test]
    fn degenerate_triangles_have_zero_area_and_zero_normal() {
        let m = SurfaceMesh::from_triangles(
            vec![[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2], [0, 1, 3]],
        )
        .unwrap();
        assert_eq!(m.face_areas()[0], 0.0, "collinear points");
        assert_eq!(m.face_normals()[0], [0.0; 3]);
        assert!(close(m.face_areas()[1], 0.5));
        // Node 2 belongs only to the degenerate triangle: zero normal, zero area.
        assert_eq!(m.vertex_normals()[2], [0.0; 3]);
        assert_eq!(m.node_areas()[2], 0.0);
    }

    #[test]
    fn vertex_normals_are_unweighted_sums_of_unit_face_normals() {
        // Node 0 sits where one large and one small triangle meet at different
        // angles; the unweighted sum differs from the area-weighted one.
        let m = SurfaceMesh::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [0.0, 10.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, -1.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 4, 3]],
        )
        .unwrap();
        let fn_ = m.face_normals();
        let expect = unit_or_zero([
            fn_[0][0] + fn_[1][0],
            fn_[0][1] + fn_[1][1],
            fn_[0][2] + fn_[1][2],
        ]);
        let got = m.vertex_normals()[0];
        for k in 0..3 {
            assert!(close(got[k], expect[k]), "{got:?} vs {expect:?}");
        }
        // Area-weighted would be dominated by the big triangle (normal +z).
        assert!(got[2] < 0.99);
    }

    #[test]
    fn bounds() {
        let b = tetra().bounds();
        assert_eq!(b.min, [-1.0; 3]);
        assert_eq!(b.max, [1.0; 3]);
        assert_eq!(b.center(), [0.0; 3]);
        assert_eq!(b.extent(), [2.0; 3]);
    }

    /// A strip of unit squares as triangles: nodes (i, 0) and (i, 1) for i in 0..=n.
    fn strip(n: u32) -> SurfaceMesh {
        let mut v = Vec::new();
        for i in 0..=n {
            v.push([i as f32, 0.0, 0.0]);
            v.push([i as f32, 1.0, 0.0]);
        }
        let mut f = Vec::new();
        for i in 0..n {
            let (a, b, c, d) = (2 * i, 2 * i + 1, 2 * i + 2, 2 * i + 3);
            f.push([a, c, b]);
            f.push([b, c, d]);
        }
        SurfaceMesh::from_triangles(v, f).unwrap()
    }

    #[test]
    fn distance_search_is_shortest_path_and_bounded() {
        let m = strip(4);
        let mut s = m.searcher();
        let found = s.within_distance(0, 1.0).unwrap();
        // From node 0 at (0,0): node 1 (0,1) at 1.0, node 2 (1,0) at 1.0; the diagonal
        // to node 3 (1,1) is sqrt(2) > 1.
        let nodes: Vec<u32> = found.iter().map(|d| d.0).collect();
        assert_eq!(nodes, vec![0, 1, 2]);
        assert_eq!(found[0], (0, 0.0));
        let wide = s.within_distance(0, 2.5).unwrap();
        // Distances are non-decreasing and shortest-path (node 4 at (2,0) is 2.0 away).
        assert!(wide.windows(2).all(|w| w[0].1 <= w[1].1));
        assert!(close(wide.iter().find(|d| d.0 == 4).unwrap().1, 2.0));
        // The searcher is reusable: the same query gives the same answer.
        assert_eq!(s.within_distance(0, 2.5).unwrap(), wide);
        // A zero radius is just the seed; a negative radius is empty.
        assert_eq!(s.within_distance(3, 0.0).unwrap(), vec![(3, 0.0)]);
        assert!(s.within_distance(3, -1.0).unwrap().is_empty());
        assert!(s.within_distance(99, 1.0).is_err());
        assert!(s.within_distance(0, f64::NAN).is_err());
    }
}
