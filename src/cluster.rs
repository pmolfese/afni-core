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
// Surface clustering: having decided which nodes are "active" (they passed a
// threshold), group the active nodes into connected clusters and keep only the
// clusters big enough to believe. This is the surface counterpart of
// `3dClusterize` and a port of SUMA's `SurfClust`, with sumaru's refinements.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * Uses `mesh.rs` (adjacency, node areas, coordinates, graph-distance search).
// * `overlay.rs` consumes the result: `ClusterLabels::survivor_mask` is exactly the
//   per-sample list its `cluster_survivors` input wants.
// * `threshold.rs` decides which nodes are active (`Threshold::passes`); this file
//   does not threshold, it only clusters.
//
// WHAT "CONNECTED" MEANS (the `Connectivity` setting; SurfClust's `-rmm`)
//
//   EdgeRings(k)       two nodes are joined if at most k mesh edges separate them
//                      (k = 1 is plain edge adjacency; SurfClust `-rmm -k`).
//   GraphDistance(r)   joined if the shortest path along mesh edges is at most r
//                      millimetres (SurfClust `-rmm r`, r > 0).
//
// The path may pass through nodes that are NOT active, which is how a wider setting
// bridges small gaps.
//
// FAITHFULNESS TO SurfClust, AND WHERE THIS DELIBERATELY DIFFERS
//
// * Clusters grow from seeds taken from the HIGHEST node index down, and ties in
//   cluster size keep that discovery order. That is SurfClust's order, so ranks
//   match its table.
// * SurfClust treats a node whose VALUE is exactly 0 as inactive, whatever the
//   threshold said. Here that is an option (`exclude_zero_values`), off by default.
// * SurfClust merges touching positive and negative regions into one cluster.
//   `Tails::Merged` does the same; `Tails::Separate` (sumaru's refinement, and what
//   `3dClusterize -bisided` does for volumes) keeps them apart.
// * SurfClust's millimetre search is a layered approximation of graph distance;
//   this crate uses true shortest paths (see the roadmap discovery log).
// * SurfClust's "central node" columns are not computed.
// ---------------------------------------------------------------------------

//! Connected-cluster labeling on a surface mesh.

use std::collections::VecDeque;

use crate::error::{Error, Result};
use crate::mesh::{NeighborhoodSearcher, SurfaceMesh};
use crate::numeric::ensure_finite;

/// How two nodes come to be in the same cluster.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Connectivity {
    /// Joined when at most this many mesh edges apart (at least 1).
    EdgeRings(u32),
    /// Joined when the shortest path along edges is at most this long, in the
    /// mesh's coordinate units (millimetres for a surface in mm). Must be finite
    /// and non-negative; 0 joins nothing, so every node is its own cluster.
    GraphDistance(f64),
}

/// How the two signs of a two-sided threshold are treated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Tails {
    /// Cluster on active status alone: positive and negative regions that touch
    /// become ONE cluster. This is SurfClust's behavior.
    #[default]
    Merged,
    /// Never join nodes whose sign differs (a positive blob touching a negative
    /// blob stays two clusters), like `3dClusterize -bisided`. The sign comes from
    /// the tail values (a node with value exactly 0 counts as positive).
    Separate,
}

/// How surviving clusters are ordered (and so numbered, rank 1 first).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ClusterSort {
    /// Largest total area first (SurfClust's default).
    #[default]
    Area,
    /// Most nodes first.
    Nodes,
    /// The order clusters were found: seeds from the highest node index down.
    Discovery,
}

/// Settings for [`label_clusters`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClusterParams {
    /// How nodes are joined.
    pub connectivity: Connectivity,
    /// How signs are treated.
    pub tails: Tails,
    /// Keep a cluster only if its area is at least this, in squared mesh units.
    pub min_area: Option<f64>,
    /// Keep a cluster only if it has at least this many nodes. When both limits are
    /// given a cluster must satisfy both (SurfClust drops a cluster failing either).
    pub min_nodes: Option<usize>,
    /// Order of the surviving clusters.
    pub sort: ClusterSort,
    /// Treat a node whose value is exactly zero as inactive (SurfClust does).
    pub exclude_zero_values: bool,
}

impl Default for ClusterParams {
    /// Edge adjacency, merged tails, no size limit, area order, zeros allowed.
    fn default() -> Self {
        Self {
            connectivity: Connectivity::EdgeRings(1),
            tails: Tails::Merged,
            min_area: None,
            min_nodes: None,
            sort: ClusterSort::Area,
            exclude_zero_values: false,
        }
    }
}

impl ClusterParams {
    fn validate(&self) -> Result<()> {
        match self.connectivity {
            Connectivity::EdgeRings(0) => {
                return Err(Error::InvalidParameter {
                    name: "edge rings".into(),
                    reason: "must be at least 1".into(),
                })
            }
            Connectivity::GraphDistance(r) => {
                ensure_finite("cluster radius", r)?;
                if r < 0.0 {
                    return Err(Error::InvalidParameter {
                        name: "cluster radius".into(),
                        reason: format!("{r} is negative"),
                    });
                }
            }
            _ => {}
        }
        if let Some(a) = self.min_area {
            ensure_finite("minimum cluster area", a)?;
        }
        Ok(())
    }
}

/// What a cluster contains. Mirrors the columns of `SurfClust`'s table that can be
/// computed from the data alone.
#[derive(Debug, Clone, PartialEq)]
pub struct ClusterSummary {
    /// Rank, starting at 1 (the value stored in [`ClusterLabels::labels`]).
    pub label: u32,
    /// The node the cluster was grown from.
    pub seed_node: u32,
    /// Number of nodes.
    pub node_count: usize,
    /// Total area: the sum of the member nodes' areas.
    pub area: f64,
    /// Mean value.
    pub mean: f64,
    /// Mean absolute value.
    pub mean_abs: f64,
    /// Smallest value and its node (the lowest node number among ties).
    pub min: (f64, u32),
    /// Largest value and its node.
    pub max: (f64, u32),
    /// Smallest absolute value and its node.
    pub min_abs: (f64, u32),
    /// Largest absolute value and its node.
    pub max_abs: (f64, u32),
    /// The node with the largest absolute value and that (signed) value: the
    /// cluster's peak. The same node as `max_abs`.
    pub peak: (u32, f64),
    /// Sample variance of the values (divides by `n - 1`; 0 for one node).
    pub variance: f64,
    /// Standard error of the mean, `sqrt(variance / n)`.
    pub std_error: f64,
    /// Value-weighted mean position: `sum(value * xyz) / sum(value)`. `NaN` if the
    /// values sum to exactly zero (the weights cancel).
    pub center_of_mass: [f64; 3],
    /// Absolute-value-weighted mean position: `sum(|value| * xyz) / sum(|value|)`.
    /// Unlike `center_of_mass` it stays well defined when positive and negative
    /// values cancel (merged tails); `NaN` only if every value is exactly zero.
    /// SurfClust has no such column, so this is an addition for viewers.
    pub center_of_mass_abs: [f64; 3],
    /// Plain mean position of the member nodes.
    pub centroid: [f64; 3],
}

/// The result of clustering.
#[derive(Debug, Clone, PartialEq)]
pub struct ClusterLabels {
    /// One entry per node: the cluster rank (1-based), or 0 for a node that is not in
    /// a surviving cluster.
    pub labels: Vec<u32>,
    /// The surviving clusters, rank 1 first.
    pub clusters: Vec<ClusterSummary>,
}

impl ClusterLabels {
    /// An empty result for a mesh of `node_count` nodes.
    pub fn empty(node_count: usize) -> Self {
        Self {
            labels: vec![0; node_count],
            clusters: Vec::new(),
        }
    }

    /// Nodes of cluster `label`, ascending (empty for 0 or an unknown label).
    pub fn nodes_for(&self, label: u32) -> Vec<u32> {
        if label == 0 {
            return Vec::new();
        }
        self.labels
            .iter()
            .enumerate()
            .filter(|(_, &l)| l == label)
            .map(|(n, _)| n as u32)
            .collect()
    }

    /// `true` for each node that belongs to a surviving cluster: the list an
    /// overlay wants as `cluster_survivors`.
    pub fn survivor_mask(&self) -> Vec<bool> {
        self.labels.iter().map(|&l| l != 0).collect()
    }
}

/// What to cluster.
#[derive(Debug, Clone, Copy)]
pub struct ClusterInput<'a> {
    /// The surface.
    pub mesh: &'a SurfaceMesh,
    /// `true` for nodes that passed the threshold (one per node).
    pub active: &'a [bool],
    /// The value at each node (one per node). Active nodes must be finite.
    pub values: &'a [f64],
    /// Values whose sign separates tails under [`Tails::Separate`]; `None` uses
    /// `values`. (Use the threshold column's values when the threshold column is not
    /// the intensity column.)
    pub tail_values: Option<&'a [f64]>,
}

/// For every node, the other nodes within the connectivity radius. Build it once
/// with [`ClusterNeighborhoods::build`] when clustering the same mesh repeatedly
/// (for example while a threshold slider moves) and pass it to
/// [`label_clusters_cached`]: the searches are then done once instead of once per
/// active node per update.
#[derive(Debug, Clone, PartialEq)]
pub struct ClusterNeighborhoods {
    connectivity: Connectivity,
    node_count: usize,
    offsets: Vec<usize>,
    nodes: Vec<u32>,
}

impl ClusterNeighborhoods {
    /// Compute the neighborhood of every node (excluding the node itself). Cost
    /// grows with the radius; memory is the total number of (node, neighbor) pairs.
    pub fn build(mesh: &SurfaceMesh, connectivity: Connectivity) -> Result<Self> {
        let params = ClusterParams {
            connectivity,
            ..Default::default()
        };
        params.validate()?;
        let n = mesh.vertices().len();
        let mut searcher = mesh.searcher();
        let mut offsets = Vec::with_capacity(n + 1);
        let mut nodes = Vec::new();
        offsets.push(0);
        for node in 0..n as u32 {
            nodes.extend(neighbors_of(mesh, &mut searcher, connectivity, node)?);
            offsets.push(nodes.len());
        }
        Ok(Self {
            connectivity,
            node_count: n,
            offsets,
            nodes,
        })
    }

    /// The neighbors of `node`.
    fn of(&self, node: u32) -> &[u32] {
        &self.nodes[self.offsets[node as usize]..self.offsets[node as usize + 1]]
    }

    /// Total number of stored (node, neighbor) pairs, a measure of memory use.
    pub fn pair_count(&self) -> usize {
        self.nodes.len()
    }
}

/// The nodes joined to `node` by `connectivity`, excluding `node`.
fn neighbors_of(
    mesh: &SurfaceMesh,
    searcher: &mut NeighborhoodSearcher<'_>,
    connectivity: Connectivity,
    node: u32,
) -> Result<Vec<u32>> {
    Ok(match connectivity {
        Connectivity::EdgeRings(1) => mesh.topology().neighbors(node)?.to_vec(),
        Connectivity::EdgeRings(k) => {
            let mut v: Vec<u32> = mesh.topology().within_rings(node, k as usize)?;
            v.retain(|&m| m != node);
            v
        }
        Connectivity::GraphDistance(r) => searcher
            .within_distance(node, r)?
            .into_iter()
            .map(|(m, _)| m)
            .filter(|&m| m != node)
            .collect(),
    })
}

/// Cluster the active nodes of `input.mesh`, computing neighborhoods as needed.
pub fn label_clusters(input: &ClusterInput<'_>, params: &ClusterParams) -> Result<ClusterLabels> {
    grow(input, params, None)
}

/// Like [`label_clusters`], reusing precomputed `neighborhoods`, which must have been
/// built for the same mesh size and the same `params.connectivity`.
pub fn label_clusters_cached(
    input: &ClusterInput<'_>,
    params: &ClusterParams,
    neighborhoods: &ClusterNeighborhoods,
) -> Result<ClusterLabels> {
    if neighborhoods.node_count != input.mesh.vertices().len()
        || neighborhoods.connectivity != params.connectivity
    {
        return Err(Error::InvalidParameter {
            name: "neighborhoods".into(),
            reason: "they were built for a different mesh or connectivity".into(),
        });
    }
    grow(input, params, Some(neighborhoods))
}

fn grow(
    input: &ClusterInput<'_>,
    params: &ClusterParams,
    cache: Option<&ClusterNeighborhoods>,
) -> Result<ClusterLabels> {
    params.validate()?;
    let n = input.mesh.vertices().len();
    let check = |what: &str, len: usize| {
        if len == n {
            Ok(())
        } else {
            Err(Error::LengthMismatch {
                what: what.into(),
                expected: n,
                found: len,
            })
        }
    };
    check("active flags", input.active.len())?;
    check("values", input.values.len())?;
    let tails_src = input.tail_values.unwrap_or(input.values);
    check("tail values", tails_src.len())?;

    // The nodes that may join a cluster.
    let mut eligible = vec![false; n];
    for (i, flag) in eligible.iter_mut().enumerate() {
        if !input.active[i] {
            continue;
        }
        ensure_finite(&format!("value of active node {i}"), input.values[i])?;
        *flag = !(params.exclude_zero_values && input.values[i] == 0.0);
    }
    let sign = |i: usize| tails_src[i] < 0.0;

    let node_areas = input.mesh.node_areas();
    let mut searcher = input.mesh.searcher();
    let mut assigned = vec![false; n];
    let mut found: Vec<(Vec<u32>, u32)> = Vec::new();

    // SurfClust starts from the highest-numbered unassigned active node.
    for seed in (0..n).rev() {
        if !eligible[seed] || assigned[seed] {
            continue;
        }
        let seed_sign = sign(seed);
        let mut members = vec![seed as u32];
        let mut queue = VecDeque::from([seed as u32]);
        assigned[seed] = true;
        while let Some(node) = queue.pop_front() {
            let owned;
            let reach: &[u32] = match cache {
                Some(c) => c.of(node),
                None => {
                    owned = neighbors_of(input.mesh, &mut searcher, params.connectivity, node)?;
                    &owned
                }
            };
            for &m in reach {
                let mi = m as usize;
                if !eligible[mi] || assigned[mi] {
                    continue;
                }
                if params.tails == Tails::Separate && sign(mi) != seed_sign {
                    continue;
                }
                assigned[mi] = true;
                members.push(m);
                queue.push_back(m);
            }
        }
        members.sort_unstable();
        found.push((members, seed as u32));
    }

    // Summaries, size filter, ordering.
    let mut kept: Vec<ClusterSummary> = Vec::new();
    let mut member_lists: Vec<Vec<u32>> = Vec::new();
    for (members, seed) in found {
        let summary = summarize(&members, seed, input, &node_areas);
        let big_enough = params.min_area.map_or(true, |a| summary.area >= a)
            && params.min_nodes.map_or(true, |k| summary.node_count >= k);
        if big_enough {
            kept.push(summary);
            member_lists.push(members);
        }
    }
    // Stable sorts, so ties keep discovery order (SurfClust's tie behavior).
    let mut order: Vec<usize> = (0..kept.len()).collect();
    match params.sort {
        ClusterSort::Area => order.sort_by(|&a, &b| kept[b].area.total_cmp(&kept[a].area)),
        ClusterSort::Nodes => order.sort_by(|&a, &b| kept[b].node_count.cmp(&kept[a].node_count)),
        ClusterSort::Discovery => {}
    }
    let mut labels = vec![0_u32; n];
    let mut clusters = Vec::with_capacity(kept.len());
    for (rank, &idx) in order.iter().enumerate() {
        let label = rank as u32 + 1;
        for &m in &member_lists[idx] {
            labels[m as usize] = label;
        }
        let mut s = kept[idx].clone();
        s.label = label;
        clusters.push(s);
    }
    Ok(ClusterLabels { labels, clusters })
}

fn summarize(
    members: &[u32],
    seed: u32,
    input: &ClusterInput<'_>,
    node_areas: &[f64],
) -> ClusterSummary {
    let n = members.len() as f64;
    let v = |m: u32| input.values[m as usize];
    let first = members[0];
    let (mut min, mut max) = ((v(first), first), (v(first), first));
    let (mut min_abs, mut max_abs) = ((v(first).abs(), first), (v(first).abs(), first));
    let (mut area, mut sum, mut sum_abs) = (0.0, 0.0, 0.0);
    let (mut weighted, mut centroid) = ([0.0_f64; 3], [0.0_f64; 3]);
    for &m in members {
        let val = v(m);
        area += node_areas[m as usize];
        sum += val;
        sum_abs += val.abs();
        // Ties keep the lowest node number because members are ascending.
        if val < min.0 {
            min = (val, m);
        }
        if val > max.0 {
            max = (val, m);
        }
        if val.abs() < min_abs.0 {
            min_abs = (val.abs(), m);
        }
        if val.abs() > max_abs.0 {
            max_abs = (val.abs(), m);
        }
        let xyz = input.mesh.vertices()[m as usize];
        for k in 0..3 {
            let c = f64::from(xyz[k]);
            weighted[k] += val * c;
            weighted_abs[k] += val.abs() * c;
            centroid[k] += c;
        }
    }
    let mean = sum / n;
    let variance = if members.len() > 1 {
        members
            .iter()
            .map(|&m| (v(m) - mean) * (v(m) - mean))
            .sum::<f64>()
            / (n - 1.0)
    } else {
        0.0
    };
    ClusterSummary {
        label: 0,
        seed_node: seed,
        node_count: members.len(),
        area,
        mean,
        mean_abs: sum_abs / n,
        min,
        max,
        min_abs,
        max_abs,
        peak: (max_abs.1, v(max_abs.1)),
        variance,
        std_error: (variance / n).sqrt(),
        center_of_mass: weighted.map(|w| if sum == 0.0 { f64::NAN } else { w / sum }),
        center_of_mass_abs: weighted_abs.map(|w| if sum_abs == 0.0 { f64::NAN } else { w / sum_abs }),
        centroid: centroid.map(|c| c / n),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A (w+1) x (h+1) grid of unit squares as triangles, in the z = 0 plane. Node
    /// (i, j) has index `j * (w + 1) + i`.
    fn grid(w: u32, h: u32) -> SurfaceMesh {
        let mut v = Vec::new();
        for j in 0..=h {
            for i in 0..=w {
                v.push([i as f32, j as f32, 0.0]);
            }
        }
        let at = |i: u32, j: u32| j * (w + 1) + i;
        let mut f = Vec::new();
        for j in 0..h {
            for i in 0..w {
                f.push([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
                f.push([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
            }
        }
        SurfaceMesh::from_triangles(v, f).unwrap()
    }

    fn run(
        mesh: &SurfaceMesh,
        active: &[bool],
        values: &[f64],
        params: &ClusterParams,
    ) -> ClusterLabels {
        label_clusters(
            &ClusterInput {
                mesh,
                active,
                values,
                tail_values: None,
            },
            params,
        )
        .unwrap()
    }

    fn flags(n: usize, on: &[usize]) -> Vec<bool> {
        let mut a = vec![false; n];
        for &i in on {
            a[i] = true;
        }
        a
    }

    #[test]
    fn separate_blobs_become_separate_clusters_ranked_by_area() {
        let m = grid(6, 1); // 14 nodes: row 0 = 0..=6, row 1 = 7..=13
        let n = 14;
        // Blob A: nodes 0, 1 (area from a corner node pair). Blob B: nodes 4, 5, 6, 11, 12, 13.
        let active = flags(n, &[0, 1, 4, 5, 6, 11, 12, 13]);
        let values = vec![1.0; n];
        let r = run(&m, &active, &values, &ClusterParams::default());
        assert_eq!(r.clusters.len(), 2);
        assert_eq!(r.clusters[0].node_count, 6, "the bigger blob is rank 1");
        assert_eq!(r.clusters[1].node_count, 2);
        assert_eq!(r.nodes_for(1), vec![4, 5, 6, 11, 12, 13]);
        assert_eq!(r.nodes_for(2), vec![0, 1]);
        assert_eq!(r.labels[2], 0, "inactive nodes are never labeled");
        assert_eq!(
            r.survivor_mask(),
            (0..n).map(|i| r.labels[i] != 0).collect::<Vec<_>>()
        );
        assert_eq!(r.nodes_for(0), Vec::<u32>::new());
    }

    #[test]
    fn rings_and_radius_bridge_gaps() {
        let m = grid(6, 1);
        let n = 14;
        let active = flags(n, &[0, 1, 4, 5]); // a one-node gap (2, 3) in row 0 between 1 and 4: two nodes
        let values = vec![1.0; n];
        let one = run(&m, &active, &values, &ClusterParams::default());
        assert_eq!(one.clusters.len(), 2);
        // Three edges reach from node 1 to node 4 along row 0; two rings do not.
        let two = ClusterParams {
            connectivity: Connectivity::EdgeRings(2),
            ..Default::default()
        };
        assert_eq!(run(&m, &active, &values, &two).clusters.len(), 2);
        let three = ClusterParams {
            connectivity: Connectivity::EdgeRings(3),
            ..Default::default()
        };
        assert_eq!(run(&m, &active, &values, &three).clusters.len(), 1);
        // The same in millimetres: nodes 1 and 4 are 3.0 apart along a row.
        let near = ClusterParams {
            connectivity: Connectivity::GraphDistance(2.9),
            ..Default::default()
        };
        let far = ClusterParams {
            connectivity: Connectivity::GraphDistance(3.0),
            ..Default::default()
        };
        assert_eq!(run(&m, &active, &values, &near).clusters.len(), 2);
        assert_eq!(run(&m, &active, &values, &far).clusters.len(), 1);
        // A zero radius joins nothing.
        let zero = ClusterParams {
            connectivity: Connectivity::GraphDistance(0.0),
            ..Default::default()
        };
        assert_eq!(run(&m, &active, &values, &zero).clusters.len(), 4);
    }

    #[test]
    fn tails_merge_or_separate() {
        let m = grid(3, 1); // nodes 0..=3 row 0, 4..=7 row 1
        let n = 8;
        let active = flags(n, &[1, 2, 5, 6]);
        let mut values = vec![0.0; n];
        values[1] = 3.0;
        values[5] = 3.0;
        values[2] = -3.0;
        values[6] = -3.0;
        let merged = run(&m, &active, &values, &ClusterParams::default());
        assert_eq!(
            merged.clusters.len(),
            1,
            "SurfClust joins touching opposite signs"
        );
        let separate = run(
            &m,
            &active,
            &values,
            &ClusterParams {
                tails: Tails::Separate,
                ..Default::default()
            },
        );
        assert_eq!(separate.clusters.len(), 2);
        assert_eq!(separate.clusters[0].node_count, 2);
        // The sign can come from a different column than the values.
        let tail = vec![1.0; n]; // all one sign: separate behaves like merged
        let r = label_clusters(
            &ClusterInput {
                mesh: &m,
                active: &active,
                values: &values,
                tail_values: Some(&tail),
            },
            &ClusterParams {
                tails: Tails::Separate,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(r.clusters.len(), 1);
    }

    #[test]
    fn size_limits_apply_to_area_and_nodes_together() {
        let m = grid(6, 1);
        let n = 14;
        let active = flags(n, &[0, 1, 4, 5, 6, 11, 12, 13]);
        let values = vec![1.0; n];
        let by_nodes = ClusterParams {
            min_nodes: Some(3),
            ..Default::default()
        };
        assert_eq!(run(&m, &active, &values, &by_nodes).clusters.len(), 1);
        // The big blob has area 3.0 (nodes 4,5,6,11,12,13: 1/3 per adjacent triangle).
        let big = run(&m, &active, &values, &ClusterParams::default()).clusters[0].area;
        let by_area = ClusterParams {
            min_area: Some(big + 0.001),
            ..Default::default()
        };
        assert_eq!(run(&m, &active, &values, &by_area).clusters.len(), 0);
        let exactly = ClusterParams {
            min_area: Some(big),
            ..Default::default()
        };
        assert_eq!(
            run(&m, &active, &values, &exactly).clusters.len(),
            1,
            "the limit is inclusive"
        );
        // Both limits: both must hold.
        let both = ClusterParams {
            min_area: Some(0.1),
            min_nodes: Some(100),
            ..Default::default()
        };
        assert!(run(&m, &active, &values, &both).clusters.is_empty());
    }

    #[test]
    fn ties_keep_discovery_order_from_the_highest_seed() {
        let m = grid(6, 1);
        let n = 14;
        // The two opposite corner nodes are mirror images: identical area.
        let active = flags(n, &[0, 13]);
        let values = vec![1.0; n];
        let r = run(&m, &active, &values, &ClusterParams::default());
        assert_eq!(r.clusters.len(), 2);
        assert!((r.clusters[0].area - r.clusters[1].area).abs() < 1e-12);
        // Seeds go from the highest index down, so node 13 is found first and keeps rank 1.
        assert_eq!((r.clusters[0].seed_node, r.clusters[1].seed_node), (13, 0));
        let nodes = run(
            &m,
            &active,
            &values,
            &ClusterParams {
                sort: ClusterSort::Nodes,
                ..Default::default()
            },
        );
        assert_eq!(
            nodes.clusters[0].seed_node, 13,
            "a tie in node count keeps the same order"
        );
        let discovery = run(
            &m,
            &active,
            &values,
            &ClusterParams {
                sort: ClusterSort::Discovery,
                ..Default::default()
            },
        );
        assert_eq!(discovery.clusters[0].seed_node, 13);
    }

    #[test]
    fn sort_modes_order_by_the_requested_measure() {
        // A line of 3 nodes on a fat patch vs a 2-node blob in a region of bigger triangles.
        let mut v = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ];
        v.extend([[10.0, 0.0, 0.0], [20.0, 0.0, 0.0], [10.0, 10.0, 0.0]]);
        let m = SurfaceMesh::from_triangles(v, vec![[0, 1, 2], [1, 3, 2], [4, 5, 6]]).unwrap();
        let active = flags(7, &[0, 1, 2, 4, 5]);
        let values = vec![1.0; 7];
        // Cluster A (nodes 0,1,2): 3 nodes. Cluster B (nodes 4,5): 2 nodes but a huge triangle.
        let by_area = run(
            &m,
            &active,
            &values,
            &ClusterParams {
                sort: ClusterSort::Area,
                ..Default::default()
            },
        );
        let by_nodes = run(
            &m,
            &active,
            &values,
            &ClusterParams {
                sort: ClusterSort::Nodes,
                ..Default::default()
            },
        );
        assert_eq!(by_area.clusters[0].node_count, 2);
        assert_eq!(by_nodes.clusters[0].node_count, 3);
    }

    #[test]
    fn summaries_report_the_statistics() {
        let m = grid(3, 1);
        let n = 8;
        let active = flags(n, &[1, 2, 5]);
        let mut values = vec![0.0; n];
        values[1] = 4.0;
        values[2] = -2.0;
        values[5] = 1.0;
        let r = run(&m, &active, &values, &ClusterParams::default());
        let c = &r.clusters[0];
        assert_eq!(c.node_count, 3);
        assert!((c.mean - 1.0).abs() < 1e-12 && (c.mean_abs - 7.0 / 3.0).abs() < 1e-12);
        assert_eq!((c.min, c.max), ((-2.0, 2), (4.0, 1)));
        assert_eq!((c.min_abs, c.max_abs), ((1.0, 5), (4.0, 1)));
        assert_eq!(c.peak, (1, 4.0));
        // Sample variance of {4, -2, 1}: mean 1, squares 9 + 9 + 0 = 18, / 2 = 9.
        assert!(
            (c.variance - 9.0).abs() < 1e-12
                && (c.std_error - (9.0_f64 / 3.0).sqrt()).abs() < 1e-12
        );
        // Centroid: mean of (1,0), (2,0), (1,1).
        assert!(
            (c.centroid[0] - 4.0 / 3.0).abs() < 1e-12 && (c.centroid[1] - 1.0 / 3.0).abs() < 1e-12
        );
        // Center of mass: sum(v * xyz) / sum(v) = (4*1 + -2*2 + 1*1, 0 + 0 + 1*1) / 3.
        assert!(
            (c.center_of_mass[0] - 1.0 / 3.0).abs() < 1e-12
                && (c.center_of_mass[1] - 1.0 / 3.0).abs() < 1e-12
        );
        // A single node has zero variance, and cancelling weights give NaN.
        let single = run(&m, &flags(n, &[1]), &values, &ClusterParams::default());
        assert_eq!(single.clusters[0].variance, 0.0);
        let cancel = {
            let mut v = vec![0.0; n];
            v[1] = 2.0;
            v[2] = -2.0;
            run(&m, &flags(n, &[1, 2]), &v, &ClusterParams::default())
        };
        assert!(cancel.clusters[0].center_of_mass[0].is_nan());
        // The absolute-weighted center survives the cancellation: equal weights at
        // nodes 1 and 2 (x = 1 and 2) put it halfway, at x = 1.5.
        assert!((cancel.clusters[0].center_of_mass_abs[0] - 1.5).abs() < 1e-12);
        // For {4, -2, 1} at x = 1, 2, 1 the weights are 4, 2, 1: (4 + 4 + 1) / 7.
        assert!((c.center_of_mass_abs[0] - 9.0 / 7.0).abs() < 1e-12);
    }

    #[test]
    fn zero_values_can_be_excluded_like_surfclust() {
        let m = grid(3, 1);
        let n = 8;
        let active = flags(n, &[1, 2, 5, 6]);
        let mut values = vec![1.0; n];
        values[2] = 0.0;
        let keep = run(&m, &active, &values, &ClusterParams::default());
        assert_eq!(keep.clusters[0].node_count, 4);
        let drop = run(
            &m,
            &active,
            &values,
            &ClusterParams {
                exclude_zero_values: true,
                ..Default::default()
            },
        );
        assert_eq!(drop.clusters.iter().map(|c| c.node_count).sum::<usize>(), 3);
        assert_eq!(drop.labels[2], 0);
    }

    #[test]
    fn bad_inputs_are_errors() {
        let m = grid(2, 1);
        let n = 6;
        let ok = flags(n, &[1]);
        let vals = vec![1.0; n];
        let p = ClusterParams::default();
        assert!(label_clusters(
            &ClusterInput {
                mesh: &m,
                active: &ok[..3],
                values: &vals,
                tail_values: None
            },
            &p
        )
        .is_err());
        assert!(label_clusters(
            &ClusterInput {
                mesh: &m,
                active: &ok,
                values: &vals[..3],
                tail_values: None
            },
            &p
        )
        .is_err());
        assert!(label_clusters(
            &ClusterInput {
                mesh: &m,
                active: &ok,
                values: &vals,
                tail_values: Some(&vals[..2])
            },
            &p
        )
        .is_err());
        let mut nan = vals.clone();
        nan[1] = f64::NAN;
        assert!(label_clusters(
            &ClusterInput {
                mesh: &m,
                active: &ok,
                values: &nan,
                tail_values: None
            },
            &p
        )
        .is_err());
        // A NaN at an INACTIVE node is irrelevant.
        nan[1] = 1.0;
        nan[4] = f64::NAN;
        assert!(label_clusters(
            &ClusterInput {
                mesh: &m,
                active: &ok,
                values: &nan,
                tail_values: None
            },
            &p
        )
        .is_ok());
        for bad in [
            Connectivity::EdgeRings(0),
            Connectivity::GraphDistance(-1.0),
            Connectivity::GraphDistance(f64::NAN),
        ] {
            let bp = ClusterParams {
                connectivity: bad,
                ..Default::default()
            };
            assert!(
                label_clusters(
                    &ClusterInput {
                        mesh: &m,
                        active: &ok,
                        values: &vals,
                        tail_values: None
                    },
                    &bp
                )
                .is_err(),
                "{bad:?}"
            );
        }
    }

    #[test]
    fn cached_neighborhoods_give_identical_results() {
        let m = grid(8, 3);
        let n = 36;
        let active = flags(n, &[0, 1, 2, 9, 10, 20, 21, 22, 30, 35]);
        let values: Vec<f64> = (0..n).map(|i| i as f64 - 10.0).collect();
        for conn in [
            Connectivity::EdgeRings(1),
            Connectivity::EdgeRings(2),
            Connectivity::GraphDistance(1.5),
            Connectivity::GraphDistance(2.2),
        ] {
            let params = ClusterParams {
                connectivity: conn,
                ..Default::default()
            };
            let direct = run(&m, &active, &values, &params);
            let cache = ClusterNeighborhoods::build(&m, conn).unwrap();
            let cached = label_clusters_cached(
                &ClusterInput {
                    mesh: &m,
                    active: &active,
                    values: &values,
                    tail_values: None,
                },
                &params,
                &cache,
            )
            .unwrap();
            assert_eq!(direct, cached, "{conn:?}");
            assert!(cache.pair_count() > 0);
        }
        // A cache for a different connectivity or mesh is refused.
        let cache = ClusterNeighborhoods::build(&m, Connectivity::EdgeRings(2)).unwrap();
        let wrong = ClusterParams {
            connectivity: Connectivity::EdgeRings(3),
            ..Default::default()
        };
        let input = ClusterInput {
            mesh: &m,
            active: &active,
            values: &values,
            tail_values: None,
        };
        assert!(label_clusters_cached(&input, &wrong, &cache).is_err());
        let other = grid(2, 2);
        assert!(ClusterNeighborhoods::build(&other, Connectivity::EdgeRings(2)).is_ok());
    }

    #[test]
    fn no_active_nodes_gives_no_clusters() {
        let m = grid(2, 2);
        let r = run(&m, &[false; 9], &[0.0; 9], &ClusterParams::default());
        assert!(r.clusters.is_empty() && r.labels.iter().all(|&l| l == 0));
        assert_eq!(r, ClusterLabels::empty(9));
    }

    #[test]
    fn clusters_feed_overlays_end_to_end() {
        use crate::color::ContinuousColorMap;
        use crate::overlay::{
            evaluate_rows, FailedThreshold, OverlayColors, OverlayInputs, OverlaySpec,
            RangeSelection,
        };
        use crate::threshold::Threshold;
        // Threshold, cluster (keep blobs of at least 3 nodes), then hide the rest.
        let m = grid(6, 1);
        let n = 14;
        let values: Vec<f64> = vec![
            5.0, 5.0, 0.0, 0.0, 5.0, 5.0, 5.0, 0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0,
        ];
        let threshold = Threshold::AbsoluteAbove(1.0);
        let active: Vec<bool> = values.iter().map(|&v| threshold.passes(v)).collect();
        let labels = run(
            &m,
            &active,
            &values,
            &ClusterParams {
                min_nodes: Some(3),
                ..Default::default()
            },
        );
        assert_eq!(
            labels.clusters.len(),
            1,
            "the 2-node blob (nodes 0, 1) is too small"
        );
        let mut spec = OverlaySpec::new(OverlayColors::Continuous(ContinuousColorMap::grayscale()));
        spec.intensity_range =
            RangeSelection::Manual(crate::column::ColumnRange::new(0.0, 5.0).unwrap());
        spec.threshold = threshold;
        spec.failed = FailedThreshold::Hide;
        let survivors = labels.survivor_mask();
        let r = evaluate_rows(
            &spec,
            &OverlayInputs {
                intensity: &values,
                cluster_survivors: Some(&survivors),
                ..Default::default()
            },
        )
        .unwrap();
        // Nodes 0 and 1 passed the threshold but their cluster was too small.
        assert!(!r.passed[0] && !r.passed[1] && r.colors[0].a == 0.0);
        assert!(r.passed[4] && r.passed[13]);
        assert_eq!(r.diagnostics.rejected_by_cluster, 2);
        assert_eq!(n, survivors.len());
    }
}
