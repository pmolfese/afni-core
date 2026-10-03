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
// Clustering of voxels in a 3D volume: threshold a statistic, join the voxels
// that survive when they touch (faces, faces+edges, or faces+edges+corners),
// drop clusters that are too small, rank the rest, and report each cluster's
// size, peak, centroid and bounding box. This is what AFNI's `3dClusterize`
// (ptaylor/3dClusterize.c, built on `NIH_find_clusters` in edt_clust.c) does.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `cluster.rs` does the same job on a SURFACE mesh. The two are kept apart on
//   purpose: a volume's neighbors come from grid arithmetic, a mesh's from its
//   triangles, and their size limits (voxel count and volume vs node count and
//   area) differ. They share vocabulary (labels, ranks, a survivor mask).
// * `domain.rs` supplies `VolumeDomain`: the grid size, the voxel order
//   (`i + nx * (j + ny * k)`) and the optional voxel-to-world matrix.
// * A caller (`afni-io`'s `Volume`, a viewer, a tool) reads the file and passes
//   plain slices of numbers here; nothing in this file touches a file.
//
// WHAT MATCHES 3dClusterize, AND WHAT IS DELIBERATELY DIFFERENT
//
// Verified against the live program (tests/volume_cluster_conformance.rs):
// * Neighbors: NN1 = 6 face neighbors, NN2 = 18 (faces + edges), NN3 = 26
//   (faces + edges + corners). Connectivity uses voxel INDEX offsets only; the
//   voxel size and the affine never change which voxels touch, even for an
//   oblique grid.
// * Thresholds are inclusive at their ends: right tail `v >= t`, left tail
//   `v <= t`, two-sided `v <= left || v >= right`, within-range `lo <= v <= hi`.
// * A voxel whose DATA value is exactly zero is never in a cluster (3dClusterize
//   clusters the non-zero voxels of the thresholded data). `exclude_zero_data`
//   keeps that default.
// * Clusters are ranked by voxel count, largest first, ties broken by discovery
//   order (the lowest voxel index first). 3dClusterize uses a stable bubble sort,
//   which gives exactly this.
// * With `Tails::Separate` ("bisided"), positive-tail and negative-tail voxels are
//   clustered independently, so a cluster never mixes them.
// * Mean is signed; SEM is `sqrt(s^2 / n)` with the sample variance; the peak is
//   the voxel with the largest absolute value (its signed value is reported), the
//   first one met when a tie occurs; the center of mass weights each voxel by its
//   absolute value.
// * Reported coordinates are the voxel centers pushed through the affine you
//   give. AFNI uses the CARDINAL (non-oblique) grid for its table, so pass the
//   cardinal matrix for byte-for-byte reports and the real (oblique) matrix for
//   true anatomical positions.
//
// Different on purpose (see the roadmap discovery log):
// * `3dClusterize -clust_vol V` is read by AFNI as "at least V VOXELS" (the volume
//   is negated and truncated like `-clust_nvox`). Here `min_volume` really is a
//   volume in cubic world units.
// * AFNI's report column labelled "Volume" prints the voxel COUNT. Here
//   `voxel_count` and `volume` are separate fields.
// * AFNI skips sorting when there are 3333 or more clusters. Here clusters are
//   always ranked.
// * AFNI compares thresholds as 32-bit floats. This API takes `f64`; to reproduce
//   AFNI on 32-bit data exactly, pass thresholds as `t as f32 as f64`.
// ---------------------------------------------------------------------------

//! Voxel clustering over a [`VolumeDomain`] (the volume counterpart of
//! [`crate::cluster`]).

use crate::domain::VolumeDomain;
use crate::error::{Error, Result};
use crate::numeric::ensure_finite;

// ---------------------------------------------------------------------------
// Connectivity
// ---------------------------------------------------------------------------

/// Which neighboring voxels count as touching. AFNI's "NN" levels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VoxelConnectivity {
    /// NN1: the 6 voxels sharing a face.
    Faces,
    /// NN2: the 18 voxels sharing a face or an edge.
    FacesEdges,
    /// NN3: all 26 surrounding voxels (face, edge or corner).
    FacesEdgesCorners,
}

impl VoxelConnectivity {
    /// The level from AFNI's `-NN` number (1, 2 or 3).
    pub fn from_nn(nn: u8) -> Result<Self> {
        match nn {
            1 => Ok(Self::Faces),
            2 => Ok(Self::FacesEdges),
            3 => Ok(Self::FacesEdgesCorners),
            _ => Err(Error::InvalidParameter {
                name: "NN".into(),
                reason: format!("{nn} is not 1, 2 or 3"),
            }),
        }
    }

    /// AFNI's `-NN` number.
    pub fn nn(self) -> u8 {
        match self {
            Self::Faces => 1,
            Self::FacesEdges => 2,
            Self::FacesEdgesCorners => 3,
        }
    }

    /// The offsets `(di, dj, dk)` to the touching voxels, in AFNI's order (`k`
    /// slowest, then `j`, then `i`; the voxel itself left out). The order matters
    /// only for which voxel is reported as the peak when two voxels tie.
    pub fn offsets(self) -> Vec<[i64; 3]> {
        // The largest squared distance (in voxel units) a neighbor may have:
        // 1 for a face, 2 for an edge, 3 for a corner.
        let max_sq = i64::from(self.nn());
        let mut out = Vec::new();
        for dk in -1..=1_i64 {
            for dj in -1..=1_i64 {
                for di in -1..=1_i64 {
                    let sq = di * di + dj * dj + dk * dk;
                    if sq > 0 && sq <= max_sq {
                        out.push([di, dj, dk]);
                    }
                }
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// Thresholds and parameters
// ---------------------------------------------------------------------------

/// Which voxels survive the voxelwise threshold. All ends are inclusive, as in
/// `3dClusterize`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VoxelThreshold {
    /// `value >= t` (`-1sided RIGHT_TAIL t`).
    RightTail(f64),
    /// `value <= t` (`-1sided LEFT_TAIL t`).
    LeftTail(f64),
    /// `value <= left_upper || value >= right_lower` (`-2sided` / `-bisided`).
    /// Whether a cluster may contain both tails is the `tails` setting of
    /// [`VolumeClusterParams`].
    TwoSided {
        /// Upper bound of the left (usually negative) tail.
        left_upper: f64,
        /// Lower bound of the right (usually positive) tail.
        right_lower: f64,
    },
    /// `lo <= value <= hi` (`-within_range lo hi`).
    WithinRange {
        /// Lower end.
        lo: f64,
        /// Upper end.
        hi: f64,
    },
}

/// Which tail a surviving voxel belongs to, for [`Tails::Separate`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    Left,
    Right,
}

impl VoxelThreshold {
    /// Check the numbers: all finite, and for a two-sided threshold
    /// `left_upper <= right_lower`; for a range `lo <= hi`.
    pub fn validate(&self) -> Result<()> {
        match *self {
            Self::RightTail(t) | Self::LeftTail(t) => ensure_finite("threshold", t).map(|_| ()),
            Self::TwoSided {
                left_upper,
                right_lower,
            } => {
                ensure_finite("left tail bound", left_upper)?;
                ensure_finite("right tail bound", right_lower)?;
                if left_upper > right_lower {
                    return Err(Error::InvalidParameter {
                        name: "two-sided threshold".into(),
                        reason: format!(
                            "the left tail bound {left_upper} is above the right tail bound {right_lower}"
                        ),
                    });
                }
                Ok(())
            }
            Self::WithinRange { lo, hi } => {
                ensure_finite("range low end", lo)?;
                ensure_finite("range high end", hi)?;
                if lo > hi {
                    return Err(Error::InvalidParameter {
                        name: "range".into(),
                        reason: format!("low end {lo} is above high end {hi}"),
                    });
                }
                Ok(())
            }
        }
    }

    /// Whether `value` survives, and on which side. NaN and infinities never do.
    fn classify(&self, value: f64) -> Option<Side> {
        if !value.is_finite() {
            return None;
        }
        match *self {
            Self::RightTail(t) => (value >= t).then_some(Side::Right),
            Self::LeftTail(t) => (value <= t).then_some(Side::Left),
            Self::TwoSided {
                left_upper,
                right_lower,
            } => {
                if value <= left_upper {
                    Some(Side::Left)
                } else if value >= right_lower {
                    Some(Side::Right)
                } else {
                    None
                }
            }
            Self::WithinRange { lo, hi } => (value >= lo && value <= hi).then_some(Side::Right),
        }
    }

    /// Whether `value` survives.
    pub fn passes(&self, value: f64) -> bool {
        self.classify(value).is_some()
    }
}

/// How the two tails of a two-sided threshold are treated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Tails {
    /// One cluster may contain voxels of both tails (`-2sided`).
    Merged,
    /// A cluster contains voxels of one tail only (`-bisided`); the tails are
    /// clustered independently. The default, because it is what is almost always
    /// wanted for a signed statistic (and matches sumaru).
    #[default]
    Separate,
}

/// How surviving clusters are ordered (and so numbered, rank 1 first).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum VolumeSort {
    /// Most voxels first; equal sizes keep their discovery order (`3dClusterize`).
    #[default]
    Size,
    /// The order clusters were found: by their lowest voxel index (for
    /// `Tails::Separate`, every right-tail cluster before any left-tail one).
    Discovery,
}

/// Settings for [`cluster_volume`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VolumeClusterParams {
    /// Which voxels touch.
    pub connectivity: VoxelConnectivity,
    /// Which voxels survive the threshold.
    pub threshold: VoxelThreshold,
    /// How tails combine (only matters for [`VoxelThreshold::TwoSided`]).
    pub tails: Tails,
    /// Keep a cluster only if it has at least this many voxels.
    pub min_voxels: Option<usize>,
    /// Keep a cluster only if its volume is at least this, in cubic world units
    /// (needs an affine on the domain). When both limits are given, a cluster must
    /// satisfy both.
    pub min_volume: Option<f64>,
    /// Treat a voxel whose data value is exactly zero as not in any cluster
    /// (what `3dClusterize` does). Turn off to let a zero-valued voxel join a
    /// cluster when it passes the threshold.
    pub exclude_zero_data: bool,
    /// Order of the surviving clusters.
    pub sort: VolumeSort,
}

impl VolumeClusterParams {
    /// Parameters with the given connectivity and threshold and `3dClusterize`'s
    /// other defaults: separate tails, no size limit, zero data excluded, size order.
    pub fn new(connectivity: VoxelConnectivity, threshold: VoxelThreshold) -> Self {
        Self {
            connectivity,
            threshold,
            tails: Tails::Separate,
            min_voxels: None,
            min_volume: None,
            exclude_zero_data: true,
            sort: VolumeSort::Size,
        }
    }

    fn validate(&self) -> Result<()> {
        self.threshold.validate()?;
        if let Some(v) = self.min_volume {
            ensure_finite("minimum cluster volume", v)?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Input and output
// ---------------------------------------------------------------------------

/// The numbers to cluster, borrowed from the caller. Every slice has one entry
/// per voxel in `VolumeDomain` order.
#[derive(Debug, Clone, Copy)]
pub struct VolumeClusterInput<'a> {
    /// The grid (and optional voxel-to-world matrix).
    pub domain: &'a VolumeDomain,
    /// The values the threshold is applied to (a statistic, usually).
    pub threshold_values: &'a [f64],
    /// The values reported for each cluster (mean, peak, center of mass). `None`
    /// means "use the threshold values", like `3dClusterize` without `-idat`.
    pub data_values: Option<&'a [f64]>,
    /// Voxels to consider; `false` removes a voxel. `None` means all voxels.
    pub mask: Option<&'a [bool]>,
}

/// One surviving cluster.
#[derive(Debug, Clone, PartialEq)]
pub struct VolumeCluster {
    /// Rank, starting at 1 (the value written in a cluster map).
    pub label: u32,
    /// The linear index of the cluster's first voxel (lowest index).
    pub seed_voxel: usize,
    /// Number of voxels.
    pub voxel_count: usize,
    /// Volume in cubic world units (`voxel_count * voxel_volume`); `None` if the
    /// domain has no affine.
    pub volume: Option<f64>,
    /// Smallest and largest `(i, j, k)` over the voxels (the bounding box in grid
    /// coordinates).
    pub bounds_ijk: ([usize; 3], [usize; 3]),
    /// Mean grid position of the voxels.
    pub centroid_ijk: [f64; 3],
    /// Grid position of the center of mass: voxel positions weighted by the
    /// absolute data value.
    pub center_of_mass_ijk: [f64; 3],
    /// `centroid_ijk` in world coordinates (`None` without an affine).
    pub centroid_world: Option<[f64; 3]>,
    /// `center_of_mass_ijk` in world coordinates (`None` without an affine).
    pub center_of_mass_world: Option<[f64; 3]>,
    /// Smallest and largest world coordinate over the voxel centers, per axis
    /// (`None` without an affine). 3dClusterize's `minRL maxRL ...` columns.
    pub bounds_world: Option<([f64; 3], [f64; 3])>,
    /// Mean of the (signed) data values.
    pub mean: f64,
    /// Mean of the absolute data values.
    pub mean_abs: f64,
    /// Standard error of the mean: `sqrt(sample variance / n)`; 0 for one voxel.
    pub std_error: f64,
    /// The peak: the voxel with the largest absolute data value (the first one met
    /// when values tie), as its linear index and signed value.
    pub peak: (usize, f64),
    /// The peak voxel's world position (`None` without an affine).
    pub peak_world: Option<[f64; 3]>,
}

/// The result of clustering a volume.
#[derive(Debug, Clone, PartialEq)]
pub struct VolumeClusters {
    /// One entry per voxel: the cluster rank (1-based), or 0 for a voxel that is
    /// not in a surviving cluster. This is `3dClusterize`'s `-pref_map` volume.
    pub labels: Vec<u32>,
    /// The surviving clusters, rank 1 first.
    pub clusters: Vec<VolumeCluster>,
}

impl VolumeClusters {
    /// The linear indices of the voxels in cluster `label`, in ascending order
    /// (empty for an unknown label).
    pub fn voxels_for(&self, label: u32) -> Vec<usize> {
        self.labels
            .iter()
            .enumerate()
            .filter(|(_, &l)| l == label && label != 0)
            .map(|(i, _)| i)
            .collect()
    }

    /// `true` for each voxel that is in a surviving cluster (a binary cluster map).
    pub fn survivor_mask(&self) -> Vec<bool> {
        self.labels.iter().map(|&l| l != 0).collect()
    }
}

// ---------------------------------------------------------------------------
// The algorithm
// ---------------------------------------------------------------------------

/// Cluster the voxels of a volume.
///
/// Steps: (1) decide which voxels are active (mask, finite values, threshold, and
/// non-zero data); (2) group active voxels that touch, one flood fill per cluster,
/// seeding from the lowest unvisited voxel index; (3) drop clusters below the size
/// limits; (4) rank the survivors; (5) summarize each.
pub fn cluster_volume(
    input: &VolumeClusterInput<'_>,
    params: &VolumeClusterParams,
) -> Result<VolumeClusters> {
    params.validate()?;
    let domain = input.domain;
    let n = domain.voxel_count();
    check_len("threshold values", input.threshold_values.len(), n)?;
    if let Some(d) = input.data_values {
        check_len("data values", d.len(), n)?;
    }
    if let Some(m) = input.mask {
        check_len("mask", m.len(), n)?;
    }
    let data: &[f64] = input.data_values.unwrap_or(input.threshold_values);

    // A volume limit needs the voxel volume, which needs an affine.
    let voxel_volume = domain.voxel_volume();
    if params.min_volume.is_some() && voxel_volume.is_none() {
        return Err(Error::InvalidParameter {
            name: "min_volume".into(),
            reason: "a minimum volume needs the domain's voxel-to-world matrix".into(),
        });
    }

    // Step 1: which voxels are active, and which tail each belongs to.
    // `None` = not active.
    let side_of: Vec<Option<Side>> = (0..n)
        .map(|v| {
            if input.mask.is_some_and(|m| !m[v]) {
                return None;
            }
            let value = data[v];
            // A non-finite data value cannot be summarized; a zero is optionally
            // excluded (3dClusterize clusters the non-zero thresholded data).
            if !value.is_finite() || (params.exclude_zero_data && value == 0.0) {
                return None;
            }
            params.threshold.classify(input.threshold_values[v])
        })
        .collect();

    // Step 2: flood fill. With separate tails, a left voxel only joins left voxels
    // and a right voxel only right ones, and every right-tail cluster is found
    // before any left-tail one (the order 3dClusterize runs its two passes).
    let merged = !matches!(params.threshold, VoxelThreshold::TwoSided { .. })
        || params.tails == Tails::Merged;
    let passes: Vec<Option<Side>> = if merged {
        vec![None]
    } else {
        vec![Some(Side::Right), Some(Side::Left)]
    };
    let offsets = params.connectivity.offsets();
    let dims = domain.dims();
    let mut taken = vec![false; n];
    let mut found: Vec<Vec<usize>> = Vec::new(); // each cluster's voxels in discovery order
    for pass_side in passes {
        // Is voxel `v` active in this pass?
        let usable = |v: usize| match (side_of[v], pass_side) {
            (None, _) => false,
            (Some(_), None) => true,
            (Some(s), Some(want)) => s == want,
        };
        for seed in 0..n {
            if taken[seed] || !usable(seed) {
                continue;
            }
            taken[seed] = true;
            let mut members = vec![seed];
            // `members` grows while we walk it: a breadth-first search, visiting
            // neighbors in AFNI's order.
            let mut at = 0;
            while at < members.len() {
                let [i, j, k] = unflatten(members[at], dims);
                at += 1;
                for off in &offsets {
                    let (ni, nj, nk) = (i + off[0], j + off[1], k + off[2]);
                    if ni < 0
                        || nj < 0
                        || nk < 0
                        || ni >= dims[0] as i64
                        || nj >= dims[1] as i64
                        || nk >= dims[2] as i64
                    {
                        continue;
                    }
                    let nv = (ni as usize) + dims[0] * ((nj as usize) + dims[1] * (nk as usize));
                    if !taken[nv] && usable(nv) {
                        taken[nv] = true;
                        members.push(nv);
                    }
                }
            }
            found.push(members);
        }
    }

    // Step 3: size limits.
    found.retain(|members| {
        let count_ok = params.min_voxels.map_or(true, |m| members.len() >= m);
        let volume_ok = match (params.min_volume, voxel_volume) {
            (Some(limit), Some(vv)) => members.len() as f64 * vv >= limit,
            _ => true,
        };
        count_ok && volume_ok
    });

    // Step 4: rank. `sort_by` is stable, so equal sizes keep discovery order.
    if params.sort == VolumeSort::Size {
        found.sort_by_key(|members| std::cmp::Reverse(members.len()));
    }

    // Step 5: label map and summaries.
    let mut labels = vec![0_u32; n];
    let mut clusters = Vec::with_capacity(found.len());
    for (rank, members) in found.iter().enumerate() {
        let label = (rank + 1) as u32;
        for &v in members {
            labels[v] = label;
        }
        clusters.push(summarize(label, members, data, domain, voxel_volume));
    }
    Ok(VolumeClusters { labels, clusters })
}

/// Error unless a slice has one entry per voxel.
fn check_len(what: &str, found: usize, expected: usize) -> Result<()> {
    if found == expected {
        Ok(())
    } else {
        Err(Error::LengthMismatch {
            what: what.into(),
            expected,
            found,
        })
    }
}

/// `(i, j, k)` of a linear index, as signed numbers so neighbor offsets can be
/// added without wrapping.
fn unflatten(index: usize, dims: [usize; 3]) -> [i64; 3] {
    [
        (index % dims[0]) as i64,
        ((index / dims[0]) % dims[1]) as i64,
        (index / (dims[0] * dims[1])) as i64,
    ]
}

/// Statistics for one cluster. `members` is in discovery order, which is what
/// decides the peak when two voxels have equal absolute values.
fn summarize(
    label: u32,
    members: &[usize],
    data: &[f64],
    domain: &VolumeDomain,
    voxel_volume: Option<f64>,
) -> VolumeCluster {
    let dims = domain.dims();
    let count = members.len() as f64;
    let (mut sum, mut sum_abs) = (0.0_f64, 0.0_f64);
    let (mut centroid, mut weighted) = ([0.0_f64; 3], [0.0_f64; 3]);
    let (mut lo, mut hi) = ([usize::MAX; 3], [0_usize; 3]);
    let mut peak = (members[0], data[members[0]]);
    let mut peak_abs = f64::NEG_INFINITY;
    let seed = members.iter().copied().min().unwrap_or(members[0]);
    for &v in members {
        let value = data[v];
        sum += value;
        sum_abs += value.abs();
        let ijk = unflatten(v, dims);
        for axis in 0..3 {
            let c = ijk[axis] as f64;
            centroid[axis] += c;
            weighted[axis] += value.abs() * c;
            lo[axis] = lo[axis].min(ijk[axis] as usize);
            hi[axis] = hi[axis].max(ijk[axis] as usize);
        }
        // Strictly greater: on a tie the first voxel met keeps the peak.
        if value.abs() > peak_abs {
            peak_abs = value.abs();
            peak = (v, value);
        }
    }
    let mean = sum / count;
    // Sample variance by the two-pass formula (stable when values are large and
    // close together, unlike sum-of-squares).
    let variance = if members.len() > 1 {
        members
            .iter()
            .map(|&v| (data[v] - mean) * (data[v] - mean))
            .sum::<f64>()
            / (count - 1.0)
    } else {
        0.0
    };
    let centroid_ijk = centroid.map(|c| c / count);
    // `sum_abs == 0` cannot happen with zero data excluded, but can when it is not.
    let center_of_mass_ijk = weighted.map(|w| {
        if sum_abs == 0.0 {
            f64::NAN
        } else {
            w / sum_abs
        }
    });

    // World-space quantities exist only with an affine.
    let world = |p: [f64; 3]| domain.ijk_to_world(p).ok();
    let bounds_world = domain.affine().map(|_| {
        let (mut wlo, mut whi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for &v in members {
            let ijk = unflatten(v, dims);
            if let Some(w) = world([ijk[0] as f64, ijk[1] as f64, ijk[2] as f64]) {
                for axis in 0..3 {
                    wlo[axis] = wlo[axis].min(w[axis]);
                    whi[axis] = whi[axis].max(w[axis]);
                }
            }
        }
        (wlo, whi)
    });
    let peak_ijk = unflatten(peak.0, dims);

    VolumeCluster {
        label,
        seed_voxel: seed,
        voxel_count: members.len(),
        volume: voxel_volume.map(|vv| count * vv),
        bounds_ijk: (lo, hi),
        centroid_ijk,
        center_of_mass_ijk,
        centroid_world: world(centroid_ijk),
        center_of_mass_world: world(center_of_mass_ijk),
        bounds_world,
        mean,
        mean_abs: sum_abs / count,
        std_error: (variance / count).sqrt(),
        peak,
        peak_world: world([peak_ijk[0] as f64, peak_ijk[1] as f64, peak_ijk[2] as f64]),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A grid with 2 x 2.5 x 3 mm voxels and no rotation, origin at zero.
    fn domain(dims: [usize; 3]) -> VolumeDomain {
        let affine = [
            [2.0, 0.0, 0.0, 0.0],
            [0.0, 2.5, 0.0, 0.0],
            [0.0, 0.0, 3.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        VolumeDomain::new(None, dims, Some(affine)).unwrap()
    }

    /// A zero volume with the given `(i, j, k, value)` voxels set.
    fn volume(dims: [usize; 3], points: &[(usize, usize, usize, f64)]) -> Vec<f64> {
        let mut v = vec![0.0; dims[0] * dims[1] * dims[2]];
        for &(i, j, k, value) in points {
            v[i + dims[0] * (j + dims[1] * k)] = value;
        }
        v
    }

    fn run(dims: [usize; 3], values: &[f64], params: &VolumeClusterParams) -> VolumeClusters {
        let d = domain(dims);
        cluster_volume(
            &VolumeClusterInput {
                domain: &d,
                threshold_values: values,
                data_values: None,
                mask: None,
            },
            params,
        )
        .unwrap()
    }

    fn right(t: f64, nn: u8) -> VolumeClusterParams {
        VolumeClusterParams::new(
            VoxelConnectivity::from_nn(nn).unwrap(),
            VoxelThreshold::RightTail(t),
        )
    }

    #[test]
    fn neighbor_counts_and_order() {
        for (nn, count) in [(1, 6), (2, 18), (3, 26)] {
            let offs = VoxelConnectivity::from_nn(nn).unwrap().offsets();
            assert_eq!(offs.len(), count, "NN{nn}");
            assert!(!offs.contains(&[0, 0, 0]));
        }
        // AFNI's order: k slowest, i fastest, so the first neighbor is below in k.
        assert_eq!(VoxelConnectivity::Faces.offsets()[0], [0, 0, -1]);
        assert!(VoxelConnectivity::from_nn(4).is_err());
    }

    #[test]
    fn diagonal_voxels_join_only_at_the_right_level() {
        let dims = [4, 4, 4];
        // Two voxels touching along a face, an edge, and a corner respectively.
        let face = volume(dims, &[(1, 1, 1, 5.0), (2, 1, 1, 5.0)]);
        let edge = volume(dims, &[(1, 1, 1, 5.0), (2, 2, 1, 5.0)]);
        let corner = volume(dims, &[(1, 1, 1, 5.0), (2, 2, 2, 5.0)]);
        for (name, v, joins_at) in [
            ("face", &face, 1),
            ("edge", &edge, 2),
            ("corner", &corner, 3),
        ] {
            for nn in 1..=3 {
                let found = run(dims, v, &right(1.0, nn));
                let expected = if nn >= joins_at { 1 } else { 2 };
                assert_eq!(found.clusters.len(), expected, "{name} at NN{nn}");
            }
        }
    }

    #[test]
    fn thresholds_include_their_ends() {
        let dims = [8, 2, 2];
        let v = volume(
            dims,
            &[
                (0, 0, 0, 5.0),
                (2, 0, 0, 4.999),
                (4, 0, 0, -5.0),
                (6, 0, 0, -4.999),
            ],
        );
        let count = |t: VoxelThreshold| {
            run(
                dims,
                &v,
                &VolumeClusterParams::new(VoxelConnectivity::Faces, t),
            )
            .clusters
            .len()
        };
        assert_eq!(count(VoxelThreshold::RightTail(5.0)), 1);
        assert_eq!(count(VoxelThreshold::LeftTail(-5.0)), 1);
        assert_eq!(
            count(VoxelThreshold::TwoSided {
                left_upper: -5.0,
                right_lower: 5.0
            }),
            2
        );
        assert_eq!(count(VoxelThreshold::WithinRange { lo: 4.999, hi: 5.0 }), 2);
    }

    #[test]
    fn merged_tails_join_signs_and_separate_tails_do_not() {
        let dims = [6, 2, 2];
        // 5, 6, -7 in a row: one cluster when merged, two when separate.
        let v = volume(dims, &[(1, 0, 0, 5.0), (2, 0, 0, 6.0), (3, 0, 0, -7.0)]);
        let t = VoxelThreshold::TwoSided {
            left_upper: -5.0,
            right_lower: 5.0,
        };
        let mut p = VolumeClusterParams::new(VoxelConnectivity::Faces, t);
        p.tails = Tails::Merged;
        let merged = run(dims, &v, &p);
        assert_eq!(merged.clusters.len(), 1);
        assert_eq!(merged.clusters[0].voxel_count, 3);
        // The peak is the largest ABSOLUTE value, reported with its sign.
        assert_eq!(merged.clusters[0].peak.1, -7.0);
        p.tails = Tails::Separate;
        let split = run(dims, &v, &p);
        assert_eq!(
            split
                .clusters
                .iter()
                .map(|c| c.voxel_count)
                .collect::<Vec<_>>(),
            vec![2, 1]
        );
    }

    #[test]
    fn clusters_are_ranked_by_size_then_discovery() {
        let dims = [10, 3, 3];
        let v = volume(
            dims,
            &[
                (0, 0, 0, 3.0), // size 1, found first
                (2, 0, 0, 3.0), // size 2
                (3, 0, 0, 3.0),
                (6, 0, 0, 3.0), // size 1, found third
            ],
        );
        let found = run(dims, &v, &right(1.0, 1));
        let sizes: Vec<usize> = found.clusters.iter().map(|c| c.voxel_count).collect();
        assert_eq!(sizes, vec![2, 1, 1]);
        // The two size-1 clusters keep their discovery order.
        assert!(found.clusters[1].seed_voxel < found.clusters[2].seed_voxel);
        // Discovery order leaves them as found.
        let mut p = right(1.0, 1);
        p.sort = VolumeSort::Discovery;
        let found = run(dims, &v, &p);
        let sizes: Vec<usize> = found.clusters.iter().map(|c| c.voxel_count).collect();
        assert_eq!(sizes, vec![1, 2, 1]);
    }

    #[test]
    fn zero_data_is_excluded_unless_allowed() {
        let dims = [4, 2, 2];
        let thr = volume(dims, &[(0, 0, 0, 1.0), (1, 0, 0, 1.0), (2, 0, 0, 1.0)]);
        // The data are zero in the middle voxel: with the 3dClusterize rule it splits.
        let data = volume(dims, &[(0, 0, 0, 2.0), (2, 0, 0, 2.0)]);
        let d = domain(dims);
        let input = VolumeClusterInput {
            domain: &d,
            threshold_values: &thr,
            data_values: Some(&data),
            mask: None,
        };
        let mut p = right(1.0, 1);
        assert_eq!(cluster_volume(&input, &p).unwrap().clusters.len(), 2);
        p.exclude_zero_data = false;
        let all = cluster_volume(&input, &p).unwrap();
        assert_eq!(all.clusters.len(), 1);
        assert_eq!(all.clusters[0].voxel_count, 3);
    }

    #[test]
    fn mask_and_size_limits() {
        let dims = [6, 2, 2];
        let v = volume(dims, &[(0, 0, 0, 3.0), (1, 0, 0, 3.0), (4, 0, 0, 3.0)]);
        let d = domain(dims);
        let mut mask = vec![true; v.len()];
        mask[1] = false; // splits the pair
        let input = VolumeClusterInput {
            domain: &d,
            threshold_values: &v,
            data_values: None,
            mask: Some(&mask),
        };
        assert_eq!(
            cluster_volume(&input, &right(1.0, 1))
                .unwrap()
                .clusters
                .len(),
            2
        );
        // Without the mask: a pair and a single; two voxels is the limit.
        let free = run(dims, &v, &{
            let mut p = right(1.0, 1);
            p.min_voxels = Some(2);
            p
        });
        assert_eq!(free.clusters.len(), 1);
        // A volume limit: voxels are 15 cubic mm, so 30 keeps the pair, 31 drops it.
        let mut p = right(1.0, 1);
        p.min_volume = Some(30.0);
        assert_eq!(run(dims, &v, &p).clusters.len(), 1);
        p.min_volume = Some(31.0);
        assert_eq!(run(dims, &v, &p).clusters.len(), 0);
    }

    #[test]
    fn a_volume_limit_needs_an_affine() {
        let dims = [3, 3, 3];
        let d = VolumeDomain::new(None, dims, None).unwrap();
        let v = vec![1.0; 27];
        let mut p = right(0.5, 1);
        p.min_volume = Some(1.0);
        let input = VolumeClusterInput {
            domain: &d,
            threshold_values: &v,
            data_values: None,
            mask: None,
        };
        assert!(cluster_volume(&input, &p).is_err());
        // Without the limit it works, and the world fields are simply absent.
        p.min_volume = None;
        let found = cluster_volume(&input, &p).unwrap();
        assert_eq!(found.clusters[0].voxel_count, 27);
        assert!(found.clusters[0].volume.is_none() && found.clusters[0].bounds_world.is_none());
    }

    #[test]
    fn summary_numbers_and_world_coordinates() {
        let dims = [4, 4, 4];
        // Values 4, -2, 1 at (1,1,1), (2,1,1), (1,2,1): same numbers as the surface test.
        let v = volume(dims, &[(1, 1, 1, 4.0), (2, 1, 1, -2.0), (1, 2, 1, 1.0)]);
        let mut p = VolumeClusterParams::new(
            VoxelConnectivity::Faces,
            VoxelThreshold::TwoSided {
                left_upper: -1.5,
                right_lower: 0.5,
            },
        );
        p.tails = Tails::Merged;
        let found = run(dims, &v, &p);
        assert_eq!(found.clusters.len(), 1);
        let c = &found.clusters[0];
        assert_eq!((c.label, c.voxel_count, c.volume), (1, 3, Some(45.0)));
        assert_eq!(c.bounds_ijk, ([1, 1, 1], [2, 2, 1]));
        assert!((c.mean - 1.0).abs() < 1e-12 && (c.mean_abs - 7.0 / 3.0).abs() < 1e-12);
        // Sample variance 9, so SEM = sqrt(9 / 3).
        assert!((c.std_error - 3.0_f64.sqrt()).abs() < 1e-12);
        assert_eq!(c.peak.1, 4.0);
        // Center of mass weights by |value|: x = (4*1 + 2*2 + 1*1) / 7.
        assert!((c.center_of_mass_ijk[0] - 9.0 / 7.0).abs() < 1e-12);
        // World = grid * (2, 2.5, 3): the centroid x is (1 + 2 + 1)/3 * 2.
        let cw = c.centroid_world.unwrap();
        assert!((cw[0] - 8.0 / 3.0).abs() < 1e-12 && (cw[2] - 3.0).abs() < 1e-12);
        let (lo, hi) = c.bounds_world.unwrap();
        assert_eq!((lo, hi), ([2.0, 2.5, 3.0], [4.0, 5.0, 3.0]));
        assert_eq!(found.voxels_for(1).len(), 3);
        assert_eq!(found.survivor_mask().iter().filter(|&&s| s).count(), 3);
    }

    #[test]
    fn peak_ties_go_to_the_first_voxel_met() {
        let dims = [4, 2, 2];
        // Equal absolute values: the seed (lowest index) is visited first.
        let v = volume(dims, &[(0, 0, 0, 3.0), (1, 0, 0, -3.0)]);
        let mut p = VolumeClusterParams::new(
            VoxelConnectivity::Faces,
            VoxelThreshold::TwoSided {
                left_upper: -1.0,
                right_lower: 1.0,
            },
        );
        p.tails = Tails::Merged;
        assert_eq!(run(dims, &v, &p).clusters[0].peak, (0, 3.0));
    }

    #[test]
    fn bad_inputs_are_errors() {
        let dims = [2, 2, 2];
        let d = domain(dims);
        let v = vec![1.0; 8];
        let wrong = vec![1.0; 7];
        assert!(cluster_volume(
            &VolumeClusterInput {
                domain: &d,
                threshold_values: &wrong,
                data_values: None,
                mask: None
            },
            &right(0.0, 1)
        )
        .is_err());
        assert!(VoxelThreshold::TwoSided {
            left_upper: 1.0,
            right_lower: -1.0
        }
        .validate()
        .is_err());
        assert!(VoxelThreshold::WithinRange { lo: 2.0, hi: 1.0 }
            .validate()
            .is_err());
        assert!(VoxelThreshold::RightTail(f64::NAN).validate().is_err());
        // NaN values are simply never active.
        let mut nan = v.clone();
        nan[0] = f64::NAN;
        let found = run(dims, &nan, &right(0.5, 3));
        assert_eq!(found.clusters[0].voxel_count, 7);
    }
}
