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
// A file-neutral boolean mask over every sample in a Domain. One entry always
// means one volume voxel or surface node; it never means "one stored row".
//
// WHY DOMAIN-SIZED MATTERS
//
// A sparse Dataset contains rows for only some domain samples, often in an
// arbitrary order. Passing a raw row-sized &[bool] between algorithms makes it
// easy to apply a mask to the wrong voxel/node. SampleMask owns the complete
// domain-sized vector and maps it into dataset row order only through `rows`.
// Missing sparse rows therefore remain false rather than shifting indices.
//
// DOMAIN IDENTITY
//
// The Domain is retained with the values. Boolean operations and dataset-row
// views require exact structural Domain equality, including a volume affine
// and any domain identifier. This catches mismatched identified surfaces and
// grids. Two anonymous surfaces with the same node count are structurally
// indistinguishable, so callers should preserve their DomainId when identity
// matters. A caller that deliberately establishes another correspondence can
// construct a new mask on the destination domain explicitly.
// ---------------------------------------------------------------------------

//! Validated masks over volume voxels or surface nodes.

use std::iter::FusedIterator;

use crate::column::{ColumnData, DataColumn};
use crate::dataset::Dataset;
use crate::domain::Domain;
use crate::error::{Error, Result};
use crate::mapping::SampleMap;

/// Rule for turning numeric values into mask membership.
///
/// Both rules exclude NaN and infinity. This matches AFNI-style mask creation,
/// where non-finite values are not valid selected samples.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MaskRule {
    /// Select every finite value other than zero.
    #[default]
    NonZero,
    /// Select every finite value strictly greater than zero.
    Positive,
}

impl MaskRule {
    /// Whether one numeric value is selected by this rule.
    pub fn selects(self, value: f64) -> bool {
        value.is_finite()
            && match self {
                Self::NonZero => value != 0.0,
                Self::Positive => value > 0.0,
            }
    }
}

/// A validated boolean mask with one value per sample in a [`Domain`].
///
/// The values always use domain order: volume linear voxel order or ascending
/// surface-node index. Dataset row order is a separate concept, especially for
/// sparse data; use [`rows`](Self::rows) when consuming a mask with a dataset.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleMask {
    domain: Domain,
    values: Vec<bool>,
}

impl SampleMask {
    /// Construct a mask from one boolean per domain sample.
    pub fn new(domain: Domain, values: Vec<bool>) -> Result<Self> {
        if values.len() != domain.sample_count() {
            return Err(Error::LengthMismatch {
                what: "sample mask values".into(),
                expected: domain.sample_count(),
                found: values.len(),
            });
        }
        Ok(Self { domain, values })
    }

    /// Construct a domain-sized mask from typed numeric values.
    ///
    /// `values` must have one entry per domain sample. Text values are rejected
    /// because their truth semantics would otherwise have to be guessed.
    pub fn from_numeric(domain: Domain, values: &ColumnData, rule: MaskRule) -> Result<Self> {
        Self::new(domain, numeric_membership(values, rule)?)
    }

    /// Construct a domain-sized mask from a dataset-row column.
    ///
    /// The column length must equal the dataset row count. Dense rows map
    /// directly to samples. Sparse rows are placed at their named domain sample,
    /// and every sample absent from the dataset is unselected. The column may be
    /// borrowed from the dataset or computed separately with the same row layout.
    pub fn from_column(dataset: &Dataset, column: &DataColumn, rule: MaskRule) -> Result<Self> {
        if column.len() != dataset.row_count() {
            return Err(Error::LengthMismatch {
                what: format!("mask column '{}' rows", column.label()),
                expected: dataset.row_count(),
                found: column.len(),
            });
        }
        let row_values = numeric_membership(column.values(), rule)?;
        let mut values = vec![false; dataset.domain().sample_count()];
        for (row, selected) in row_values.into_iter().enumerate() {
            let sample = dataset
                .sample_for_row(row)
                .ok_or_else(|| Error::InvalidParameter {
                    name: "dataset map".into(),
                    reason: format!("row {row} has no domain sample"),
                })? as usize;
            values[sample] = selected;
        }
        Self::new(dataset.domain().clone(), values)
    }

    /// Domain whose voxels or nodes these values describe.
    pub fn domain(&self) -> &Domain {
        &self.domain
    }

    /// Boolean membership in domain-sample order.
    pub fn values(&self) -> &[bool] {
        &self.values
    }

    /// Number of domain samples represented by the mask.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether the domain contains no samples.
    ///
    /// Valid core domains are nonempty, so a successfully constructed mask
    /// always returns `false`. The method is provided alongside [`len`](Self::len)
    /// for conventional collection-style APIs.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Number of selected samples.
    pub fn count(&self) -> usize {
        self.values.iter().filter(|&&selected| selected).count()
    }

    /// Whether at least one sample is selected.
    pub fn any(&self) -> bool {
        self.values.iter().any(|&selected| selected)
    }

    /// Whether every domain sample is selected.
    pub fn all(&self) -> bool {
        self.values.iter().all(|&selected| selected)
    }

    /// Whether domain sample `sample` is selected.
    pub fn is_selected(&self, sample: usize) -> Result<bool> {
        self.values
            .get(sample)
            .copied()
            .ok_or(Error::IndexOutOfRange {
                index: i64::try_from(sample).unwrap_or(i64::MAX),
                len: self.values.len(),
            })
    }

    /// Indices of selected domain samples, in ascending order.
    pub fn selected_samples(&self) -> impl Iterator<Item = usize> + '_ {
        self.values
            .iter()
            .enumerate()
            .filter_map(|(sample, &selected)| selected.then_some(sample))
    }

    /// Require exact structural compatibility with a domain.
    ///
    /// Equality includes sample count and type, any domain identifier, and a
    /// volume affine when present. Anonymous surfaces with the same node count
    /// cannot be distinguished; preserve their [`DomainId`](crate::domain::DomainId)
    /// when identity matters. This is still stricter than checking only the
    /// length of a raw boolean slice.
    pub fn require_domain(&self, domain: &Domain) -> Result<()> {
        if &self.domain == domain {
            Ok(())
        } else {
            Err(Error::DomainMismatch {
                what: "sample mask and requested domain differ".into(),
            })
        }
    }

    /// Iterate mask membership in a dataset's stored row order.
    ///
    /// The dataset must have the same domain. For sparse data this follows the
    /// dataset's index list, so the iterator length is the stored row count, not
    /// the full mask length.
    pub fn rows<'a>(&'a self, dataset: &'a Dataset) -> Result<SampleMaskRows<'a>> {
        self.require_domain(dataset.domain())?;
        Ok(SampleMaskRows {
            values: &self.values,
            map: dataset.map(),
            next: 0,
            rows: dataset.row_count(),
        })
    }

    /// Intersection: selected only where both masks are selected.
    pub fn and(&self, other: &Self) -> Result<Self> {
        self.combine(other, "intersection", |left, right| left && right)
    }

    /// Union: selected where either mask is selected.
    pub fn or(&self, other: &Self) -> Result<Self> {
        self.combine(other, "union", |left, right| left || right)
    }

    /// Complement over this mask's complete domain.
    pub fn inverted(&self) -> Self {
        Self {
            domain: self.domain.clone(),
            values: self.values.iter().map(|selected| !selected).collect(),
        }
    }

    /// Shared checked implementation for two-mask boolean operations.
    fn combine(
        &self,
        other: &Self,
        operation: &str,
        combine: impl Fn(bool, bool) -> bool,
    ) -> Result<Self> {
        if self.domain != other.domain {
            return Err(Error::DomainMismatch {
                what: format!("cannot compute mask {operation} across different domains"),
            });
        }
        Ok(Self {
            domain: self.domain.clone(),
            values: self
                .values
                .iter()
                .copied()
                .zip(other.values.iter().copied())
                .map(|(left, right)| combine(left, right))
                .collect(),
        })
    }
}

/// Mask membership in one dataset's stored row order.
///
/// Returned by [`SampleMask::rows`]. It is exact-sized even when the dataset is
/// sparse or its rows are not in ascending sample order.
#[derive(Debug, Clone)]
pub struct SampleMaskRows<'a> {
    values: &'a [bool],
    map: &'a SampleMap,
    next: usize,
    rows: usize,
}

impl Iterator for SampleMaskRows<'_> {
    type Item = bool;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next == self.rows {
            return None;
        }
        let row = self.next;
        self.next += 1;
        // Dataset construction validated every row-to-sample mapping, while
        // SampleMask::rows required the same Domain before constructing us.
        let sample = self
            .map
            .sample_for_row(row)
            .expect("validated dataset row mapping") as usize;
        Some(self.values[sample])
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.rows - self.next;
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for SampleMaskRows<'_> {}
impl FusedIterator for SampleMaskRows<'_> {}

/// Apply a rule to every value of a numeric typed column.
fn numeric_membership(values: &ColumnData, rule: MaskRule) -> Result<Vec<bool>> {
    if !values.is_numeric() {
        return Err(Error::InvalidParameter {
            name: "mask values".into(),
            reason: "text values cannot define numeric mask membership".into(),
        });
    }
    Ok((0..values.len())
        .map(|row| {
            rule.selects(
                values
                    .get_f64(row)
                    .expect("validated numeric column and in-range row"),
            )
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::column::ColumnRole;
    use crate::domain::{DomainId, SurfaceDomain, VolumeDomain};

    fn surface(samples: usize) -> Domain {
        Domain::Surface(SurfaceDomain::new(None, samples).unwrap())
    }

    #[test]
    fn construction_checks_domain_length_and_sample_indices() {
        let domain = surface(4);
        assert!(matches!(
            SampleMask::new(domain.clone(), vec![true; 3]),
            Err(Error::LengthMismatch {
                expected: 4,
                found: 3,
                ..
            })
        ));

        let mask = SampleMask::new(domain, vec![false, true, false, true]).unwrap();
        assert_eq!((mask.len(), mask.count()), (4, 2));
        assert!(!mask.is_empty());
        assert!(mask.any());
        assert!(!mask.all());
        assert!(mask.is_selected(1).unwrap());
        assert_eq!(mask.selected_samples().collect::<Vec<_>>(), [1, 3]);
        assert!(matches!(
            mask.is_selected(4),
            Err(Error::IndexOutOfRange { index: 4, len: 4 })
        ));
    }

    #[test]
    fn numeric_rules_exclude_nonfinite_values_and_reject_text() {
        let values = ColumnData::Float64(vec![0.0, -2.0, 3.0, f64::NAN, f64::INFINITY]);
        let nonzero = SampleMask::from_numeric(surface(5), &values, MaskRule::NonZero).unwrap();
        assert_eq!(nonzero.values(), [false, true, true, false, false]);
        let positive = SampleMask::from_numeric(surface(5), &values, MaskRule::Positive).unwrap();
        assert_eq!(positive.values(), [false, false, true, false, false]);

        let text = ColumnData::Text(vec!["yes".into()]);
        assert!(matches!(
            SampleMask::from_numeric(surface(1), &text, MaskRule::NonZero),
            Err(Error::InvalidParameter { ref name, .. }) if name == "mask values"
        ));
    }

    #[test]
    fn sparse_columns_expand_to_domain_samples_and_rows_map_back() {
        let dataset = Dataset::indexed(
            crate::dataset::DatasetKind::Scalar,
            surface(7),
            vec![5, 1, 3],
            vec![DataColumn::new(
                "threshold",
                ColumnRole::Generic,
                ColumnData::Float32(vec![2.0, 0.0, -1.0]),
            )
            .unwrap()],
        )
        .unwrap();
        let mask =
            SampleMask::from_column(&dataset, &dataset.columns()[0], MaskRule::NonZero).unwrap();

        // Samples 5 and 3 are selected; absent samples and row 1's zero are not.
        assert_eq!(
            mask.values(),
            [false, false, false, true, false, true, false]
        );
        // Dataset row order is [sample 5, sample 1, sample 3].
        assert_eq!(
            mask.rows(&dataset).unwrap().collect::<Vec<_>>(),
            [true, false, true]
        );
    }

    #[test]
    fn domain_identity_is_checked_even_when_lengths_match() {
        let surface_mask = SampleMask::new(surface(2), vec![true, false]).unwrap();
        let volume = Domain::Volume(VolumeDomain::new(None, [2, 1, 1], None).unwrap());
        let volume_mask = SampleMask::new(volume.clone(), vec![false, true]).unwrap();

        assert!(matches!(
            surface_mask.or(&volume_mask),
            Err(Error::DomainMismatch { .. })
        ));
        assert!(matches!(
            surface_mask.require_domain(&volume),
            Err(Error::DomainMismatch { .. })
        ));

        let left_surface =
            Domain::Surface(SurfaceDomain::new(Some(DomainId::new("left").unwrap()), 2).unwrap());
        let right_surface =
            Domain::Surface(SurfaceDomain::new(Some(DomainId::new("right").unwrap()), 2).unwrap());
        let identified_mask = SampleMask::new(left_surface, vec![true, false]).unwrap();
        assert!(matches!(
            identified_mask.require_domain(&right_surface),
            Err(Error::DomainMismatch { .. })
        ));

        let identity = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let shifted = [
            [1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let identity_volume =
            Domain::Volume(VolumeDomain::new(None, [2, 1, 1], Some(identity)).unwrap());
        let shifted_volume =
            Domain::Volume(VolumeDomain::new(None, [2, 1, 1], Some(shifted)).unwrap());
        let affine_mask = SampleMask::new(identity_volume, vec![true, false]).unwrap();
        assert!(matches!(
            affine_mask.require_domain(&shifted_volume),
            Err(Error::DomainMismatch { .. })
        ));
    }

    #[test]
    fn boolean_operations_preserve_the_domain() {
        let domain = surface(4);
        let left = SampleMask::new(domain.clone(), vec![true, true, false, false]).unwrap();
        let right = SampleMask::new(domain.clone(), vec![true, false, true, false]).unwrap();

        assert_eq!(
            left.and(&right).unwrap().values(),
            [true, false, false, false]
        );
        assert_eq!(left.or(&right).unwrap().values(), [true, true, true, false]);
        assert_eq!(left.inverted().values(), [false, false, true, true]);
        assert_eq!(left.and(&right).unwrap().domain(), &domain);
    }
}
