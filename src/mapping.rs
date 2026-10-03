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
// `SampleMap`: how the ROWS of a dataset relate to the SAMPLES of its domain.
//
//   Dense    row i is sample i                 (every sample has a row)
//   Indexed  row r is sample indices[r]        (only some samples have rows)
//
// Indexed ("sparse") data is common for surfaces: a cluster map may list only
// the 300 nodes that survive thresholding. We keep it as a short index list
// and expand to the full domain only when a consumer asks (`expand`).
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `dataset.rs` owns one `SampleMap` and checks column lengths against it.
// * `domain.rs` supplies the sample count the indices are validated against.
// * `numeric.rs` supplies the checked conversions used here.
//
// POLICIES (documented here because they are decisions, not accidents)
//
// * Out-of-domain index  -> error (`IndexOutOfRange`).
// * Duplicate index      -> error (`DuplicateIndex`); never "last one wins".
// * Unsorted indices     -> allowed; row order is preserved exactly.
// * Missing samples      -> samples with no row are simply absent; `expand`
//                           fills them with a caller-chosen value.
// ---------------------------------------------------------------------------

//! Dense and indexed row-to-sample mappings.

use crate::error::{Error, Result};

/// How dataset rows map onto domain samples.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SampleMap {
    /// Row `i` is sample `i`, for `i` in `0..rows`. `rows` equals the domain's
    /// sample count (checked when a [`Dataset`](crate::dataset::Dataset) is built).
    Dense {
        /// Number of rows.
        rows: usize,
    },
    /// Row `r` is sample `indices[r]`.
    Indexed {
        /// Sample index for each row (unique, each below the domain size).
        indices: Vec<u32>,
    },
}

impl SampleMap {
    /// A dense map of `rows` rows.
    pub fn dense(rows: usize) -> Self {
        SampleMap::Dense { rows }
    }

    /// An indexed map, validated against a domain of `domain_len` samples.
    ///
    /// Fails with [`Error::IndexOutOfRange`] if any index is `>= domain_len`
    /// and [`Error::DuplicateIndex`] if any index repeats.
    pub fn indexed(indices: Vec<u32>, domain_len: usize) -> Result<Self> {
        // A bit-set over the domain: `seen[n]` is true once index n is used.
        // This finds duplicates in O(rows + domain) without sorting, so the
        // caller's row order is untouched.
        let mut seen = vec![false; domain_len];
        for &n in &indices {
            let slot = seen.get_mut(n as usize).ok_or(Error::IndexOutOfRange {
                index: i64::from(n),
                len: domain_len,
            })?;
            if std::mem::replace(slot, true) {
                return Err(Error::DuplicateIndex { index: n });
            }
        }
        Ok(SampleMap::Indexed { indices })
    }

    /// Number of rows.
    pub fn row_count(&self) -> usize {
        match self {
            SampleMap::Dense { rows } => *rows,
            SampleMap::Indexed { indices } => indices.len(),
        }
    }

    /// Whether this is an indexed (sparse) map.
    pub fn is_indexed(&self) -> bool {
        matches!(self, SampleMap::Indexed { .. })
    }

    /// Whether row `i` is sample `i` for every row and the rows cover a domain of
    /// `domain_len` samples exactly. True for a dense map of that length and for
    /// an indexed map listing `0, 1, 2, ...` in order.
    ///
    /// AFNI writes even fully populated surface datasets with an index list,
    /// so this identifies data that is dense in all but representation.
    pub fn is_identity(&self, domain_len: usize) -> bool {
        match self {
            SampleMap::Dense { rows } => *rows == domain_len,
            SampleMap::Indexed { indices } => {
                indices.len() == domain_len
                    && indices
                        .iter()
                        .enumerate()
                        .all(|(row, &n)| n as usize == row)
            }
        }
    }

    /// The sample that row `row` describes, or `None` if `row` is out of range.
    pub fn sample_for_row(&self, row: usize) -> Option<u32> {
        match self {
            SampleMap::Dense { rows } => {
                // `try_from` instead of `as`: a row beyond u32 is "no sample".
                (row < *rows).then(|| u32::try_from(row).ok()).flatten()
            }
            SampleMap::Indexed { indices } => indices.get(row).copied(),
        }
    }

    /// The row that holds `sample`, or `None` if that sample has no row.
    ///
    /// Linear time for an indexed map; if you need many lookups, build a dense
    /// view with [`expand`](Self::expand) or [`row_lookup`](Self::row_lookup).
    pub fn row_for_sample(&self, sample: u32) -> Option<usize> {
        match self {
            SampleMap::Dense { rows } => ((sample as usize) < *rows).then_some(sample as usize),
            SampleMap::Indexed { indices } => indices.iter().position(|&n| n == sample),
        }
    }

    /// A table `lookup[sample] = Some(row)` over a domain of `domain_len`
    /// samples, for fast repeated lookups.
    pub fn row_lookup(&self, domain_len: usize) -> Vec<Option<usize>> {
        let mut lookup = vec![None; domain_len];
        for row in 0..self.row_count() {
            if let Some(sample) = self.sample_for_row(row) {
                if let Some(slot) = lookup.get_mut(sample as usize) {
                    *slot = Some(row);
                }
            }
        }
        lookup
    }

    /// Expand per-row `values` into one value per domain sample, using `fill`
    /// for samples that have no row. This is the dense view of sparse data.
    ///
    /// `values.len()` must equal [`row_count`](Self::row_count).
    pub fn expand<T: Clone>(&self, values: &[T], domain_len: usize, fill: T) -> Result<Vec<T>> {
        if values.len() != self.row_count() {
            return Err(Error::LengthMismatch {
                what: "values to expand".into(),
                expected: self.row_count(),
                found: values.len(),
            });
        }
        let mut out = vec![fill; domain_len];
        for (row, value) in values.iter().enumerate() {
            let sample = self.sample_for_row(row).map(|n| n as usize);
            // Dense maps are validated against the domain elsewhere; stay safe
            // here by reporting rather than panicking on a bad index.
            match sample.and_then(|n| out.get_mut(n)) {
                Some(slot) => *slot = value.clone(),
                None => {
                    return Err(Error::IndexOutOfRange {
                        index: sample.map_or(-1, |n| n as i64),
                        len: domain_len,
                    })
                }
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dense_maps_rows_to_themselves() {
        let m = SampleMap::dense(3);
        assert_eq!(m.row_count(), 3);
        assert!(!m.is_indexed());
        assert_eq!(m.sample_for_row(2), Some(2));
        assert_eq!(m.sample_for_row(3), None);
        assert_eq!(m.row_for_sample(1), Some(1));
        assert_eq!(m.row_for_sample(3), None);
    }

    #[test]
    fn indexed_preserves_order_and_looks_up_both_ways() {
        let m = SampleMap::indexed(vec![7, 2, 5], 10).unwrap();
        assert_eq!(m.row_count(), 3);
        assert_eq!(m.sample_for_row(0), Some(7));
        assert_eq!(m.row_for_sample(5), Some(2));
        assert_eq!(m.row_for_sample(3), None);
        let lookup = m.row_lookup(10);
        assert_eq!(lookup[7], Some(0));
        assert_eq!(lookup[0], None);
    }

    #[test]
    fn indexed_rejects_out_of_domain_and_duplicates() {
        assert_eq!(
            SampleMap::indexed(vec![0, 10], 10),
            Err(Error::IndexOutOfRange { index: 10, len: 10 })
        );
        assert_eq!(
            SampleMap::indexed(vec![1, 4, 1], 10),
            Err(Error::DuplicateIndex { index: 1 })
        );
        // An empty index list is a legal (empty) map; the dataset layer
        // decides whether zero rows are acceptable.
        assert_eq!(SampleMap::indexed(vec![], 10).unwrap().row_count(), 0);
    }

    #[test]
    fn identity_detects_dense_in_disguise() {
        assert!(SampleMap::dense(3).is_identity(3));
        assert!(!SampleMap::dense(3).is_identity(4));
        assert!(SampleMap::indexed(vec![0, 1, 2], 3).unwrap().is_identity(3));
        assert!(!SampleMap::indexed(vec![0, 1, 2], 5).unwrap().is_identity(5));
        assert!(!SampleMap::indexed(vec![1, 0, 2], 3).unwrap().is_identity(3));
    }

    #[test]
    fn expand_fills_missing_samples() {
        let m = SampleMap::indexed(vec![3, 1], 5).unwrap();
        let dense = m.expand(&[30.0, 10.0], 5, f64::NAN).unwrap();
        assert_eq!(dense[1], 10.0);
        assert_eq!(dense[3], 30.0);
        assert!(dense[0].is_nan() && dense[2].is_nan() && dense[4].is_nan());
        assert!(m.expand(&[1.0], 5, 0.0).is_err());
    }

    #[test]
    fn dense_expand_is_identity() {
        let m = SampleMap::dense(3);
        assert_eq!(m.expand(&[1, 2, 3], 3, 0).unwrap(), vec![1, 2, 3]);
        // A dense map longer than the domain is reported, not a panic.
        assert!(m.expand(&[1, 2, 3], 2, 0).is_err());
    }
}
