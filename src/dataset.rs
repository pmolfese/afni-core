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
// `Dataset`: a domain (surface or volume), a row-to-sample mapping, and a list
// of typed columns. This is the format-neutral object the rest of afni-core
// (statistics, overlays, clusters, ROIs) operates on. Constructors validate, so
// a `Dataset` that exists is internally consistent.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * Combines `domain.rs` (where), `mapping.rs` (which samples have rows) and
//   `column.rs` (the values and their metadata).
// * Ported from `sumaru/src/dataset.rs`, generalised so volumes work too and
//   with NaN/duplicate/out-of-domain policies made explicit.
// * `afni-io` converts files into this type (see `afni_io::adapt`); it never
//   appears in file syntax here.
//
// INVARIANTS ENFORCED BY `Dataset::new`
//
//   1. at least one column;
//   2. every column has the same number of rows, and that equals the map's;
//   3. a dense map has exactly one row per domain sample;
//   4. an indexed map's indices are in-domain and unique;
//   5. time step, if present, is finite and positive; time start is finite.
// ---------------------------------------------------------------------------

//! The format-neutral dataset model.

use crate::column::{DataColumn, MissingFill};
use crate::domain::Domain;
use crate::error::{Error, Result};
use crate::mapping::SampleMap;
use crate::numeric::ensure_finite;

/// What sort of data a dataset holds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DatasetKind {
    /// Ordinary scalar values (statistics, thickness, ...).
    Scalar,
    /// Integer label keys.
    Label,
    /// A time series, one column per time point.
    TimeSeries,
    /// Region-of-interest membership.
    Roi,
    /// The source did not say.
    Unknown,
    /// A kind this crate has no name for; the source's text is kept.
    Other(String),
}

/// Identifiers that tie a dataset to the things it came from.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ParentIds {
    /// The dataset's own id.
    pub self_id: Option<String>,
    /// Id of the surface/grid the rows refer to.
    pub domain_parent: Option<String>,
    /// Id of the surface whose geometry the dataset came from.
    pub geometry_parent: Option<String>,
}

/// A validated dataset. See the module docs for the invariants.
#[derive(Debug, Clone, PartialEq)]
pub struct Dataset {
    kind: DatasetKind,
    domain: Domain,
    map: SampleMap,
    columns: Vec<DataColumn>,
    time_step_seconds: Option<f64>,
    time_start_seconds: Option<f64>,
    parent_ids: ParentIds,
}

impl Dataset {
    /// Build and validate a dataset.
    pub fn new(
        kind: DatasetKind,
        domain: Domain,
        map: SampleMap,
        columns: Vec<DataColumn>,
    ) -> Result<Self> {
        if columns.is_empty() {
            return Err(Error::Empty("dataset (no columns)".into()));
        }
        // Invariant 2: every column matches the map's row count.
        let rows = map.row_count();
        for c in &columns {
            if c.len() != rows {
                return Err(Error::LengthMismatch {
                    what: format!("column '{}' rows", c.label()),
                    expected: rows,
                    found: c.len(),
                });
            }
        }
        // Invariants 3 and 4. An indexed map is re-validated here because the
        // map may have been built against a different domain size.
        match &map {
            SampleMap::Dense { rows } => {
                if *rows != domain.sample_count() {
                    return Err(Error::LengthMismatch {
                        what: "dense dataset rows vs domain samples".into(),
                        expected: domain.sample_count(),
                        found: *rows,
                    });
                }
            }
            SampleMap::Indexed { indices } => {
                SampleMap::indexed(indices.clone(), domain.sample_count())?;
            }
        }
        Ok(Self {
            kind,
            domain,
            map,
            columns,
            time_step_seconds: None,
            time_start_seconds: None,
            parent_ids: ParentIds::default(),
        })
    }

    /// Convenience: a dense dataset (one row per domain sample).
    pub fn dense(kind: DatasetKind, domain: Domain, columns: Vec<DataColumn>) -> Result<Self> {
        let map = SampleMap::dense(domain.sample_count());
        Self::new(kind, domain, map, columns)
    }

    /// Convenience: a sparse dataset listing the samples that have rows.
    pub fn indexed(
        kind: DatasetKind,
        domain: Domain,
        indices: Vec<u32>,
        columns: Vec<DataColumn>,
    ) -> Result<Self> {
        let map = SampleMap::indexed(indices, domain.sample_count())?;
        Self::new(kind, domain, map, columns)
    }

    /// Set the time step in seconds. Must be finite and greater than zero.
    pub fn with_time_step_seconds(mut self, step: Option<f64>) -> Result<Self> {
        if let Some(s) = step {
            ensure_finite("time step", s)?;
            if s <= 0.0 {
                return Err(Error::InvalidParameter {
                    name: "time step".into(),
                    reason: format!("must be positive, got {s}"),
                });
            }
        }
        self.time_step_seconds = step;
        Ok(self)
    }

    /// Set the time of the first time-series column, in seconds. Must be finite.
    pub fn with_time_start_seconds(mut self, start: Option<f64>) -> Result<Self> {
        if let Some(s) = start {
            ensure_finite("time start", s)?;
        }
        self.time_start_seconds = start;
        Ok(self)
    }

    /// Attach parent identifiers.
    pub fn with_parent_ids(mut self, ids: ParentIds) -> Self {
        self.parent_ids = ids;
        self
    }

    /// What kind of data this is.
    pub fn kind(&self) -> &DatasetKind {
        &self.kind
    }
    /// The domain the rows live on.
    pub fn domain(&self) -> &Domain {
        &self.domain
    }
    /// The row-to-sample mapping.
    pub fn map(&self) -> &SampleMap {
        &self.map
    }
    /// The columns, in order.
    pub fn columns(&self) -> &[DataColumn] {
        &self.columns
    }
    /// Number of rows.
    pub fn row_count(&self) -> usize {
        self.map.row_count()
    }
    /// Seconds between time-series columns, if known.
    pub fn time_step_seconds(&self) -> Option<f64> {
        self.time_step_seconds
    }
    /// Time of the first time-series column, if known.
    pub fn time_start_seconds(&self) -> Option<f64> {
        self.time_start_seconds
    }
    /// Parent identifiers.
    pub fn parent_ids(&self) -> &ParentIds {
        &self.parent_ids
    }

    /// Whether rows are stored sparsely.
    pub fn is_sparse(&self) -> bool {
        self.map.is_indexed()
    }

    /// The domain sample described by `row`.
    pub fn sample_for_row(&self, row: usize) -> Option<u32> {
        self.map.sample_for_row(row)
    }

    /// The first column whose label is exactly `label`. (Labels need not be
    /// unique; AFNI permits duplicates.)
    pub fn column(&self, label: &str) -> Option<&DataColumn> {
        self.columns.iter().find(|c| c.label() == label)
    }

    /// All columns with the given role.
    pub fn columns_with_role<'a>(
        &'a self,
        role: &'a crate::column::ColumnRole,
    ) -> impl Iterator<Item = &'a DataColumn> {
        self.columns.iter().filter(move |c| c.role() == role)
    }

    /// Equality of the whole dataset with NaN equal to NaN.
    ///
    /// Plain `==` is IEEE-faithful, so a dataset containing missing (NaN)
    /// samples is never `==` to its own clone. This compares everything else
    /// exactly and values with [`ColumnData::eq_nan_aware`](crate::column::ColumnData::eq_nan_aware).
    pub fn eq_nan_aware(&self, other: &Self) -> bool {
        self.kind == other.kind
            && self.domain == other.domain
            && self.map == other.map
            && self.time_step_seconds == other.time_step_seconds
            && self.time_start_seconds == other.time_start_seconds
            && self.parent_ids == other.parent_ids
            && self.columns.len() == other.columns.len()
            && self
                .columns
                .iter()
                .zip(&other.columns)
                .all(|(a, b)| a.eq_nan_aware(b))
    }

    /// An equivalent dataset with a dense map: one row per domain sample, with
    /// `fill` standing in for samples that had no row. A dense dataset is
    /// returned as an unchanged clone.
    ///
    /// Row order of the result is sample order, so round-tripping a sparse
    /// dataset through `to_dense` and reading the original rows back by sample
    /// yields identical values.
    pub fn to_dense(&self, fill: &MissingFill) -> Result<Self> {
        if !self.is_sparse() {
            return Ok(self.clone());
        }
        let n = self.domain.sample_count();
        let columns = self
            .columns
            .iter()
            .map(|c| c.with_values(c.values().expand(&self.map, n, fill)?))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            map: SampleMap::dense(n),
            columns,
            ..self.clone()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::column::{ColumnData, ColumnRole};
    use crate::domain::{SurfaceDomain, VolumeDomain};

    fn surface(n: usize) -> Domain {
        Domain::Surface(SurfaceDomain::new(None, n).unwrap())
    }

    fn col(label: &str, v: Vec<f32>) -> DataColumn {
        DataColumn::new(label, ColumnRole::Generic, ColumnData::Float32(v)).unwrap()
    }

    #[test]
    fn dense_dataset_must_cover_the_domain() {
        assert!(Dataset::dense(
            DatasetKind::Scalar,
            surface(3),
            vec![col("a", vec![1.0; 3])]
        )
        .is_ok());
        assert!(Dataset::dense(
            DatasetKind::Scalar,
            surface(4),
            vec![col("a", vec![1.0; 3])]
        )
        .is_err());
    }

    #[test]
    fn columns_must_agree_and_exist() {
        assert!(Dataset::dense(DatasetKind::Scalar, surface(3), vec![]).is_err());
        let err = Dataset::dense(
            DatasetKind::Scalar,
            surface(3),
            vec![col("a", vec![1.0; 3]), col("b", vec![1.0; 2])],
        );
        assert!(matches!(err, Err(Error::LengthMismatch { .. })));
    }

    #[test]
    fn sparse_validates_indices() {
        let ok = Dataset::indexed(
            DatasetKind::Scalar,
            surface(10),
            vec![9, 0],
            vec![col("a", vec![1.0, 2.0])],
        );
        assert!(ok.unwrap().is_sparse());
        let oob = Dataset::indexed(
            DatasetKind::Scalar,
            surface(10),
            vec![10],
            vec![col("a", vec![1.0])],
        );
        assert!(matches!(oob, Err(Error::IndexOutOfRange { .. })));
        let dup = Dataset::indexed(
            DatasetKind::Scalar,
            surface(10),
            vec![1, 1],
            vec![col("a", vec![1.0, 2.0])],
        );
        assert!(matches!(dup, Err(Error::DuplicateIndex { index: 1 })));
        let short = Dataset::indexed(
            DatasetKind::Scalar,
            surface(10),
            vec![1, 2],
            vec![col("a", vec![1.0])],
        );
        assert!(matches!(short, Err(Error::LengthMismatch { .. })));
    }

    #[test]
    fn time_parameters_are_checked_not_silently_dropped() {
        let d = || {
            Dataset::dense(
                DatasetKind::TimeSeries,
                surface(2),
                vec![col("t0", vec![0.0; 2])],
            )
            .unwrap()
        };
        assert_eq!(
            d().with_time_step_seconds(Some(2.0))
                .unwrap()
                .time_step_seconds(),
            Some(2.0)
        );
        assert!(d().with_time_step_seconds(Some(0.0)).is_err());
        assert!(d().with_time_step_seconds(Some(f64::NAN)).is_err());
        assert!(d().with_time_start_seconds(Some(f64::INFINITY)).is_err());
        assert_eq!(
            d().with_time_start_seconds(Some(-4.0))
                .unwrap()
                .time_start_seconds(),
            Some(-4.0)
        );
    }

    #[test]
    fn dense_and_sparse_are_equivalent() {
        // Sparse rows for samples 3 and 1 of a 5-node surface.
        let sparse = Dataset::indexed(
            DatasetKind::Scalar,
            surface(5),
            vec![3, 1],
            vec![col("a", vec![30.0, 10.0])],
        )
        .unwrap();
        let dense = sparse.to_dense(&MissingFill::default()).unwrap();
        assert!(!dense.is_sparse());
        assert_eq!(dense.row_count(), 5);
        let v = dense.columns()[0].values();
        // Every row of the sparse set appears at its sample in the dense one.
        for row in 0..sparse.row_count() {
            let sample = sparse.sample_for_row(row).unwrap() as usize;
            assert_eq!(sparse.columns()[0].values().get_f64(row), v.get_f64(sample));
        }
        // Samples without rows are NaN under the default fill.
        assert!(v.get_f64(0).unwrap().is_nan());
        // Densifying a dense dataset changes nothing. Derived `==` follows IEEE
        // (NaN != NaN), so use `eq_nan_aware` for data with missing samples.
        let zero = MissingFill {
            float: 0.0,
            ..MissingFill::default()
        };
        let dense0 = sparse.to_dense(&zero).unwrap();
        assert_eq!(dense0.to_dense(&zero).unwrap(), dense0);
    }

    #[test]
    fn volumes_use_the_same_model() {
        let dom = Domain::Volume(VolumeDomain::new(None, [2, 2, 2], None).unwrap());
        let d = Dataset::dense(DatasetKind::Scalar, dom, vec![col("vol0", vec![0.0; 8])]).unwrap();
        assert_eq!(d.row_count(), 8);
        assert_eq!(d.column("vol0").unwrap().label(), "vol0");
        assert!(d.column("nope").is_none());
    }

    #[test]
    fn role_filter_finds_columns() {
        let a =
            DataColumn::new("a", ColumnRole::Intensity, ColumnData::Float32(vec![1.0])).unwrap();
        let b =
            DataColumn::new("b", ColumnRole::Threshold, ColumnData::Float32(vec![2.0])).unwrap();
        let d = Dataset::dense(DatasetKind::Scalar, surface(1), vec![a, b]).unwrap();
        let hits: Vec<_> = d
            .columns_with_role(&ColumnRole::Threshold)
            .map(|c| c.label())
            .collect();
        assert_eq!(hits, ["b"]);
    }
}
