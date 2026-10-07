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
// "Domains": the set of places a dataset has values for. A surface domain is a
// set of mesh nodes; a volume domain is a 3D voxel grid. A domain says WHAT the
// samples are and how many there are. It does not hold the values.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `mapping.rs` validates row-to-sample lists against a domain's sample count.
// * `dataset.rs` pairs a domain with columns of values.
// * Later: Phase 6 builds surface topology on `SurfaceDomain`, Phase 7 builds
//   voxel neighborhoods on `VolumeDomain`.
// * Surfaces and volumes are defined together on purpose, so the API does not
//   bake in one viewer's assumptions (roadmap Phase 1).
// ---------------------------------------------------------------------------

//! Sample domains: surface node sets and volume voxel grids.

use crate::error::{Error, Result};

/// A stable identifier for a domain, normally an AFNI ID code
/// (`domain_parent_idcode` in a NIML dataset, `IDCODE_STRING` in a `.HEAD`).
///
/// Two domains with the same id are *claimed* to be the same underlying
/// surface or grid. Equal node counts alone never imply that (roadmap
/// Phase 10): only a shared id, or an explicit caller policy, does.
///
/// This is a "newtype": a one-field struct that gives a plain `String` its own
/// type, so an id cannot be mixed up with an arbitrary string by accident.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct DomainId(String);

impl DomainId {
    /// Wrap an id string. Surrounding whitespace is trimmed; an empty id is an
    /// error (use `Option<DomainId>` for "no id").
    pub fn new(id: impl Into<String>) -> Result<Self> {
        let id = id.into();
        let trimmed = id.trim();
        if trimmed.is_empty() {
            return Err(Error::Empty("domain id".into()));
        }
        Ok(Self(trimmed.to_owned()))
    }

    /// The id text.
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// A set of surface nodes `0..node_count`.
///
/// Topology (which nodes connect) arrives in Phase 6; for now a surface domain
/// is just a count plus an optional identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SurfaceDomain {
    id: Option<DomainId>,
    node_count: usize,
}

impl SurfaceDomain {
    /// A surface with `node_count` nodes. A surface must have at least one.
    pub fn new(id: Option<DomainId>, node_count: usize) -> Result<Self> {
        if node_count == 0 {
            return Err(Error::Empty("surface domain".into()));
        }
        Ok(Self { id, node_count })
    }

    /// The surface's identity, if known.
    pub fn id(&self) -> Option<&DomainId> {
        self.id.as_ref()
    }

    /// Number of nodes.
    pub fn node_count(&self) -> usize {
        self.node_count
    }
}

/// Convert a position between AFNI's DICOM world frame (what FATCAT, `3dClusterize`
/// and AFNI files use: +x toward the patient's left, +y posterior, which AFNI calls
/// "RAI") and the RAS frame NIfTI and GIFTI use (+x right, +y anterior). The two
/// differ by the sign of x and y, so the same function converts both ways.
///
/// Only the axis signs change; nothing is resampled. A position that is already in
/// the other frame comes out wrong, so keep track of which one you hold.
pub fn flip_dicom_ras([x, y, z]: [f32; 3]) -> [f32; 3] {
    [-x, -y, z]
}

/// A 3D voxel grid. Voxel `(i, j, k)` has linear index `i + nx * (j + ny * k)`,
/// the same order AFNI and NIfTI use.
#[derive(Debug, Clone, PartialEq)]
pub struct VolumeDomain {
    id: Option<DomainId>,
    dims: [usize; 3],
    /// Voxel index `(i, j, k, 1)` to world `(x, y, z)`, row-major 4x4, if known.
    affine: Option<[[f64; 4]; 4]>,
}

impl VolumeDomain {
    /// A grid of `dims = [nx, ny, nz]` voxels. Every dimension must be at least 1,
    /// the voxel count must not overflow, and any affine entry must be finite.
    pub fn new(
        id: Option<DomainId>,
        dims: [usize; 3],
        affine: Option<[[f64; 4]; 4]>,
    ) -> Result<Self> {
        if dims.contains(&0) {
            return Err(Error::Empty("volume domain".into()));
        }
        // `checked_mul` returns None on overflow instead of wrapping.
        dims[0]
            .checked_mul(dims[1])
            .and_then(|n| n.checked_mul(dims[2]))
            .ok_or_else(|| Error::InvalidParameter {
                name: "dims".into(),
                reason: format!("{dims:?} voxels overflow usize"),
            })?;
        if let Some(m) = &affine {
            for v in m.iter().flatten() {
                crate::numeric::ensure_finite("affine entry", *v)?;
            }
        }
        Ok(Self { id, dims, affine })
    }

    /// The grid's identity, if known.
    pub fn id(&self) -> Option<&DomainId> {
        self.id.as_ref()
    }

    /// Grid size `[nx, ny, nz]`.
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// Voxel-to-world matrix, if the source supplied one.
    pub fn affine(&self) -> Option<&[[f64; 4]; 4]> {
        self.affine.as_ref()
    }

    /// World position `(x, y, z)` of a voxel center given by (possibly fractional)
    /// grid coordinates `(i, j, k)`: the affine applied to `(i, j, k, 1)`.
    ///
    /// Errors if the domain has no affine. The point does not have to lie inside
    /// the grid (this is plain arithmetic, so fractional and out-of-range
    /// coordinates are fine).
    pub fn ijk_to_world(&self, ijk: [f64; 3]) -> Result<[f64; 3]> {
        let m = self
            .affine
            .as_ref()
            .ok_or_else(|| Error::InvalidParameter {
                name: "affine".into(),
                reason: "this volume domain has no voxel-to-world matrix".into(),
            })?;
        // Row r of the matrix dotted with (i, j, k, 1).
        Ok(std::array::from_fn(|r| {
            m[r][0] * ijk[0] + m[r][1] * ijk[1] + m[r][2] * ijk[2] + m[r][3]
        }))
    }

    /// Fractional voxel coordinates corresponding to a world position.
    ///
    /// This is the inverse of [`ijk_to_world`](Self::ijk_to_world) for an
    /// invertible affine. The returned coordinates are not rounded or clipped;
    /// negative and out-of-grid locations are valid geometric results. A
    /// missing, singular, or numerically degenerate affine is an error.
    pub fn world_to_ijk(&self, world: [f64; 3]) -> Result<[f64; 3]> {
        let affine = self
            .affine
            .as_ref()
            .ok_or_else(|| Error::InvalidParameter {
                name: "affine".into(),
                reason: "this volume domain has no voxel-to-world matrix".into(),
            })?;
        let inverse = inverse_linear_part(affine)?;
        let translated: [f64; 3] = std::array::from_fn(|axis| world[axis] - affine[axis][3]);
        Ok(std::array::from_fn(|row| {
            inverse[row][0] * translated[0]
                + inverse[row][1] * translated[1]
                + inverse[row][2] * translated[2]
        }))
    }

    /// Describe a cropped sub-grid without moving or resampling values.
    ///
    /// `start` is the source-grid voxel that becomes `(0, 0, 0)` and `dims` is
    /// the new size. The affine translation is shifted so world coordinates of
    /// retained voxels remain unchanged. Because the set of voxels changes,
    /// the returned domain intentionally has no domain identifier.
    pub fn crop(&self, start: [usize; 3], dims: [usize; 3]) -> Result<Self> {
        for axis in 0..3 {
            let end =
                start[axis]
                    .checked_add(dims[axis])
                    .ok_or_else(|| Error::InvalidParameter {
                        name: "crop bounds".into(),
                        reason: format!("axis {axis} start plus length overflows"),
                    })?;
            if dims[axis] == 0 || end > self.dims[axis] {
                return Err(Error::InvalidParameter {
                    name: "crop bounds".into(),
                    reason: format!(
                        "axis {axis} range {}..{} is outside 0..{}",
                        start[axis], end, self.dims[axis]
                    ),
                });
            }
        }
        let affine = self
            .affine
            .as_ref()
            .map(|affine| shifted_affine(affine, start.map(|value| value as f64)));
        Self::new(None, dims, affine)
    }

    /// Describe a zero-padded grid without allocating or copying values.
    ///
    /// `before` and `after` give added voxels on the low and high side of each
    /// axis. The old `(0, 0, 0)` voxel moves to `before`, and the affine origin
    /// is shifted so that voxel keeps the same world position. The returned
    /// domain has no domain identifier because its sample set has changed.
    pub fn pad(&self, before: [usize; 3], after: [usize; 3]) -> Result<Self> {
        let mut dims = [0; 3];
        for axis in 0..3 {
            dims[axis] = self.dims[axis]
                .checked_add(before[axis])
                .and_then(|size| size.checked_add(after[axis]))
                .ok_or_else(|| Error::InvalidParameter {
                    name: "padding".into(),
                    reason: format!("axis {axis} size overflows"),
                })?;
        }
        let offset = before.map(|value| -(value as f64));
        let affine = self
            .affine
            .as_ref()
            .map(|affine| shifted_affine(affine, offset));
        Self::new(None, dims, affine)
    }

    /// Volume of one voxel in cubic world units (cubic millimetres for an affine
    /// in millimetres): the absolute determinant of the affine's 3x3 part. `None`
    /// if the domain has no affine. For a rotated or oblique grid this is still the
    /// true voxel volume, because rotation does not change a determinant's size.
    pub fn voxel_volume(&self) -> Option<f64> {
        let m = self.affine.as_ref()?;
        // Determinant of the upper-left 3x3 block, expanded along the first row.
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        Some(det.abs())
    }

    /// Total number of voxels.
    pub fn voxel_count(&self) -> usize {
        // Cannot overflow: checked in `new`.
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    /// Linear index of voxel `(i, j, k)`, or an error if it is outside the grid.
    ///
    /// Takes signed coordinates so a negative neighbor offset reports
    /// [`Error::IndexOutOfRange`] instead of wrapping.
    pub fn linear_index(&self, i: i64, j: i64, k: i64) -> Result<usize> {
        let [nx, ny, nz] = self.dims;
        // Check each axis separately so the error names the offending
        // coordinate and its axis length. `checked_index` rejects negatives
        // before converting, so nothing can wrap.
        let i = crate::numeric::checked_index(i, nx)?;
        let j = crate::numeric::checked_index(j, ny)?;
        let k = crate::numeric::checked_index(k, nz)?;
        let idx = i + nx * (j + ny * k);
        Ok(idx)
    }

    /// Inverse of [`linear_index`](Self::linear_index): `(i, j, k)` of a linear
    /// index.
    pub fn ijk(&self, index: usize) -> Result<[usize; 3]> {
        if index >= self.voxel_count() {
            return Err(Error::IndexOutOfRange {
                index: i64::try_from(index).unwrap_or(i64::MAX),
                len: self.voxel_count(),
            });
        }
        let [nx, ny, _] = self.dims;
        Ok([index % nx, (index / nx) % ny, index / (nx * ny)])
    }
}

fn shifted_affine(affine: &[[f64; 4]; 4], offset: [f64; 3]) -> [[f64; 4]; 4] {
    let mut shifted = *affine;
    for row in 0..3 {
        shifted[row][3] = affine[row][3]
            + affine[row][0] * offset[0]
            + affine[row][1] * offset[1]
            + affine[row][2] * offset[2];
    }
    shifted
}

fn inverse_linear_part(affine: &[[f64; 4]; 4]) -> Result<[[f64; 3]; 3]> {
    let scale = affine[..3]
        .iter()
        .flat_map(|row| row[..3].iter())
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    if scale == 0.0 || !scale.is_finite() {
        return Err(singular_affine_error());
    }
    let a: [[f64; 3]; 3] =
        std::array::from_fn(|row| std::array::from_fn(|column| affine[row][column] / scale));
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if !det.is_finite() || det.abs() <= 64.0 * f64::EPSILON {
        return Err(singular_affine_error());
    }
    let inverse_scale = 1.0 / (det * scale);
    Ok([
        [
            (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * inverse_scale,
            (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * inverse_scale,
            (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inverse_scale,
        ],
        [
            (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * inverse_scale,
            (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inverse_scale,
            (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * inverse_scale,
        ],
        [
            (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * inverse_scale,
            (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * inverse_scale,
            (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inverse_scale,
        ],
    ])
}

fn singular_affine_error() -> Error {
    Error::InvalidParameter {
        name: "affine".into(),
        reason: "voxel-to-world linear transform is singular or numerically degenerate".into(),
    }
}

/// The domain a dataset lives on: a surface or a volume.
#[derive(Debug, Clone, PartialEq)]
pub enum Domain {
    /// A surface mesh's nodes.
    Surface(SurfaceDomain),
    /// A voxel grid.
    Volume(VolumeDomain),
}

impl Domain {
    /// How many samples the domain holds (nodes or voxels).
    pub fn sample_count(&self) -> usize {
        match self {
            Domain::Surface(s) => s.node_count(),
            Domain::Volume(v) => v.voxel_count(),
        }
    }

    /// The domain's identity, if known.
    pub fn id(&self) -> Option<&DomainId> {
        match self {
            Domain::Surface(s) => s.id(),
            Domain::Volume(v) => v.id(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids_are_trimmed_and_nonempty() {
        assert_eq!(DomainId::new("  XYZ_1 ").unwrap().as_str(), "XYZ_1");
        assert!(DomainId::new("   ").is_err());
    }

    #[test]
    fn surface_needs_nodes() {
        assert!(SurfaceDomain::new(None, 0).is_err());
        assert_eq!(SurfaceDomain::new(None, 42).unwrap().node_count(), 42);
    }

    #[test]
    fn volume_validates_dims_and_affine() {
        assert!(VolumeDomain::new(None, [4, 0, 6], None).is_err());
        assert!(VolumeDomain::new(None, [usize::MAX, 2, 2], None).is_err());
        let mut bad = [[0.0; 4]; 4];
        bad[0][0] = f64::NAN;
        assert!(VolumeDomain::new(None, [1, 1, 1], Some(bad)).is_err());
        assert_eq!(
            VolumeDomain::new(None, [4, 5, 6], None)
                .unwrap()
                .voxel_count(),
            120
        );
    }

    #[test]
    fn linear_index_and_ijk_are_inverses_and_checked() {
        let v = VolumeDomain::new(None, [4, 5, 6], None).unwrap();
        assert_eq!(v.linear_index(0, 0, 0).unwrap(), 0);
        assert_eq!(v.linear_index(3, 0, 0).unwrap(), 3);
        assert_eq!(v.linear_index(0, 1, 0).unwrap(), 4);
        assert_eq!(v.linear_index(0, 0, 1).unwrap(), 20);
        for idx in 0..v.voxel_count() {
            let [i, j, k] = v.ijk(idx).unwrap();
            assert_eq!(v.linear_index(i as i64, j as i64, k as i64).unwrap(), idx);
        }
        assert!(v.linear_index(4, 0, 0).is_err());
        assert!(v.linear_index(0, -1, 0).is_err());
        assert!(v.ijk(120).is_err());
    }

    #[test]
    fn domain_enum_reports_sample_counts() {
        let s = Domain::Surface(SurfaceDomain::new(None, 10).unwrap());
        let v = Domain::Volume(VolumeDomain::new(None, [2, 3, 4], None).unwrap());
        assert_eq!((s.sample_count(), v.sample_count()), (10, 24));
    }

    #[test]
    fn world_positions_and_voxel_volume() {
        // 2 x 2.5 x 3 mm voxels, origin (10, 20, 30).
        let affine = [
            [2.0, 0.0, 0.0, 10.0],
            [0.0, 2.5, 0.0, 20.0],
            [0.0, 0.0, 3.0, 30.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let d = VolumeDomain::new(None, [4, 4, 4], Some(affine)).unwrap();
        assert_eq!(d.ijk_to_world([1.0, 2.0, 3.0]).unwrap(), [12.0, 25.0, 39.0]);
        // Fractional and out-of-grid positions are plain arithmetic.
        assert_eq!(
            d.ijk_to_world([-1.0, 0.5, 0.0]).unwrap(),
            [8.0, 21.25, 30.0]
        );
        assert_eq!(d.voxel_volume(), Some(15.0));
        // A rotation (here 90 degrees about z, axes swapped and one flipped) keeps
        // the volume, and a mirror (negative determinant) still gives a positive one.
        let rotated = [
            [0.0, -2.5, 0.0, 0.0],
            [2.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, -3.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let r = VolumeDomain::new(None, [4, 4, 4], Some(rotated)).unwrap();
        assert_eq!(r.voxel_volume(), Some(15.0));
        // Without a matrix there is no world position and no volume.
        let none = VolumeDomain::new(None, [2, 2, 2], None).unwrap();
        assert!(none.ijk_to_world([0.0; 3]).is_err());
        assert_eq!(none.voxel_volume(), None);
    }

    #[test]
    fn world_to_ijk_inverts_oblique_affines_and_rejects_singular_ones() {
        let affine = [
            [2.0, 0.1, 0.0, -12.0],
            [0.0, -3.0, 0.2, 9.0],
            [0.0, 0.0, 4.0, 2.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let domain = VolumeDomain::new(None, [8, 9, 10], Some(affine)).unwrap();
        let ijk = [4.25, -0.5, 6.75];
        let back = domain
            .world_to_ijk(domain.ijk_to_world(ijk).unwrap())
            .unwrap();
        for (actual, expected) in back.into_iter().zip(ijk) {
            assert!((actual - expected).abs() < 1.0e-12);
        }

        let singular = VolumeDomain::new(
            None,
            [2, 2, 2],
            Some([
                [1.0, 2.0, 3.0, 0.0],
                [2.0, 4.0, 6.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]),
        )
        .unwrap();
        assert!(singular.world_to_ijk([0.0; 3]).is_err());
    }

    #[test]
    fn crop_and_pad_keep_retained_voxels_at_the_same_world_positions() {
        let affine = [
            [2.0, 0.25, 0.0, 10.0],
            [0.0, -3.0, 0.1, 20.0],
            [0.0, 0.0, 4.0, 30.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let source = VolumeDomain::new(
            Some(DomainId::new("source").unwrap()),
            [10, 11, 12],
            Some(affine),
        )
        .unwrap();

        let cropped = source.crop([2, 3, 4], [5, 6, 7]).unwrap();
        assert_eq!(cropped.dims(), [5, 6, 7]);
        assert_eq!(cropped.id(), None);
        assert_eq!(
            cropped.ijk_to_world([0.0; 3]).unwrap(),
            source.ijk_to_world([2.0, 3.0, 4.0]).unwrap()
        );
        assert!(source.crop([8, 0, 0], [3, 1, 1]).is_err());

        let padded = source.pad([2, 3, 4], [1, 2, 3]).unwrap();
        assert_eq!(padded.dims(), [13, 16, 19]);
        assert_eq!(padded.id(), None);
        assert_eq!(
            padded.ijk_to_world([2.0, 3.0, 4.0]).unwrap(),
            source.ijk_to_world([0.0; 3]).unwrap()
        );
    }
}
