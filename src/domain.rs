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
}
