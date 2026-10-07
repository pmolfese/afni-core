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
// Reusable descriptive reductions over one numeric Dataset column. This is the
// common implementation for the loops that AFNI C programs repeatedly write to
// count values, find extrema, and compute means and variances inside a mask.
//
// WHY THE DATASET IS PART OF THE API
//
// A sparse dataset row is not necessarily the same number as its surface node
// or volume voxel. Accepting Dataset rather than just &[f64] lets extrema report
// both coordinates and lets SampleMask remain domain-sized without guessing how
// it maps to stored rows.
// ---------------------------------------------------------------------------

//! Mask-aware descriptive reductions over typed dataset columns.

use crate::dataset::Dataset;
use crate::error::{Error, Result};
use crate::mask::SampleMask;
use crate::numeric::NonFinitePolicy;

/// Denominator used for the variance in a [`ColumnSummary`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VarianceNormalization {
    /// Divide by `n`, treating the included values as the complete population.
    Population,
    /// Divide by `n - 1`, estimating variance from a sample.
    ///
    /// A summary containing fewer than two included values has no sample
    /// variance and returns `None` for that field.
    Sample,
}

/// Explicit policies controlling a column summary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ColumnSummaryOptions {
    /// How NaN and infinity are handled.
    ///
    /// [`NonFinitePolicy::Skip`] omits them from arithmetic but still counts
    /// them. [`NonFinitePolicy::Reject`] returns an error at the first one.
    /// [`NonFinitePolicy::Propagate`] includes them: NaN poisons every numeric
    /// result, infinity participates in the sum and extrema, and variance is
    /// NaN because it is not finite.
    pub non_finite: NonFinitePolicy,
    /// Whether variance uses a population or sample denominator.
    pub variance: VarianceNormalization,
}

impl ColumnSummaryOptions {
    /// Construct options while requiring both scientific choices explicitly.
    pub const fn new(non_finite: NonFinitePolicy, variance: VarianceNormalization) -> Self {
        Self {
            non_finite,
            variance,
        }
    }
}

/// One column extremum and where it occurred.
///
/// Ties use the first value in stored row order. For a dense dataset `row` and
/// `sample` are equal; for sparse data `sample` is the actual voxel or node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ColumnExtremum {
    /// Numeric value at the extremum.
    pub value: f64,
    /// Zero-based row in each stored dataset column.
    pub row: usize,
    /// Zero-based sample in the complete volume or surface domain.
    pub sample: usize,
}

/// Counts and descriptive statistics for one numeric dataset column.
///
/// `selected_count` counts stored rows admitted by the optional mask. A domain
/// sample absent from a sparse dataset has no value and is therefore not
/// counted even if the mask selects it. `included_count` is the number that
/// actually entered arithmetic after applying the non-finite policy.
#[derive(Debug, Clone, PartialEq)]
pub struct ColumnSummary {
    /// Stored rows selected by the mask, before non-finite handling.
    pub selected_count: usize,
    /// Selected finite values.
    pub finite_count: usize,
    /// Selected NaN or infinite values.
    pub non_finite_count: usize,
    /// Values included in arithmetic after applying the non-finite policy.
    pub included_count: usize,
    /// Selected finite values equal to positive or negative zero.
    pub zero_count: usize,
    /// Selected finite values unequal to zero.
    pub nonzero_count: usize,
    /// Sum of included values, or `None` when none were included.
    pub sum: Option<f64>,
    /// Arithmetic mean of included values, or `None` when none were included.
    pub mean: Option<f64>,
    /// Variance using the requested normalization.
    ///
    /// This is `None` for no included values and for sample variance with fewer
    /// than two. It is NaN when propagated non-finite input makes variance
    /// undefined.
    pub variance: Option<f64>,
    /// Smallest included value and its location.
    pub min: Option<ColumnExtremum>,
    /// Largest included value and its location.
    pub max: Option<ColumnExtremum>,
}

impl ColumnSummary {
    /// Standard deviation derived from [`variance`](Self::variance).
    pub fn standard_deviation(&self) -> Option<f64> {
        self.variance.map(f64::sqrt)
    }
}

/// Summarize one numeric dataset column, optionally inside a domain mask.
///
/// The column and mask are validated before iteration. Values are processed in
/// stored row order, and Welford's online algorithm is used for finite variance
/// to avoid the cancellation in `mean(x*x) - mean(x)^2`. Integer storage is
/// converted through [`ColumnData::get_f64`](crate::column::ColumnData::get_f64),
/// so an `i64` outside the exactly representable `f64` range may lose precision.
pub fn summarize_column(
    dataset: &Dataset,
    column_index: usize,
    mask: Option<&SampleMask>,
    options: ColumnSummaryOptions,
) -> Result<ColumnSummary> {
    let column = dataset.column_at(column_index)?;
    if let Some(mask) = mask {
        mask.require_domain(dataset.domain())?;
    }
    if !column.values().is_numeric() {
        return Err(Error::InvalidParameter {
            name: "summary column".into(),
            reason: format!("column '{}' is text, not numeric", column.label()),
        });
    }

    let mut accumulator = SummaryAccumulator::new(options);
    for row in 0..dataset.row_count() {
        let sample = dataset
            .sample_for_row(row)
            .expect("validated dataset maps every stored row") as usize;
        if mask.is_some_and(|mask| !mask.values()[sample]) {
            continue;
        }
        let value = column
            .values()
            .get_f64(row)
            .expect("validated numeric column and in-range row");
        accumulator.push(value, row, sample, column.label())?;
    }
    Ok(accumulator.finish())
}

/// Internal one-pass state kept separate from the public result.
#[derive(Debug)]
struct SummaryAccumulator {
    options: ColumnSummaryOptions,
    selected_count: usize,
    finite_count: usize,
    non_finite_count: usize,
    included_count: usize,
    zero_count: usize,
    nonzero_count: usize,
    finite_sum: f64,
    mean: f64,
    squared_deviations: f64,
    min: Option<ColumnExtremum>,
    max: Option<ColumnExtremum>,
    first_nan: Option<ColumnExtremum>,
    positive_infinity: bool,
    negative_infinity: bool,
}

impl SummaryAccumulator {
    fn new(options: ColumnSummaryOptions) -> Self {
        Self {
            options,
            selected_count: 0,
            finite_count: 0,
            non_finite_count: 0,
            included_count: 0,
            zero_count: 0,
            nonzero_count: 0,
            finite_sum: 0.0,
            mean: 0.0,
            squared_deviations: 0.0,
            min: None,
            max: None,
            first_nan: None,
            positive_infinity: false,
            negative_infinity: false,
        }
    }

    fn push(&mut self, value: f64, row: usize, sample: usize, label: &str) -> Result<()> {
        self.selected_count += 1;
        let location = ColumnExtremum { value, row, sample };

        if !value.is_finite() {
            self.non_finite_count += 1;
            match self.options.non_finite {
                NonFinitePolicy::Reject => {
                    return Err(Error::NonFinite {
                        what: format!("column '{label}' row {row}"),
                        value,
                    });
                }
                NonFinitePolicy::Skip => return Ok(()),
                NonFinitePolicy::Propagate => {
                    self.included_count += 1;
                    if value.is_nan() {
                        self.first_nan.get_or_insert(location);
                    } else if value.is_sign_positive() {
                        self.positive_infinity = true;
                        update_extrema(&mut self.min, location, true);
                        update_extrema(&mut self.max, location, false);
                    } else {
                        self.negative_infinity = true;
                        update_extrema(&mut self.min, location, true);
                        update_extrema(&mut self.max, location, false);
                    }
                    return Ok(());
                }
            }
        }

        self.finite_count += 1;
        self.included_count += 1;
        if value == 0.0 {
            self.zero_count += 1;
        } else {
            self.nonzero_count += 1;
        }

        // Welford update: stable even when the mean is much larger than the
        // variation. Only finite values reach this state.
        let finite_n = self.finite_count as f64;
        let delta = value - self.mean;
        self.mean += delta / finite_n;
        let delta_after = value - self.mean;
        self.squared_deviations += delta * delta_after;
        self.finite_sum += value;

        update_extrema(&mut self.min, location, true);
        update_extrema(&mut self.max, location, false);
        Ok(())
    }

    fn finish(self) -> ColumnSummary {
        let propagated_non_finite =
            self.options.non_finite == NonFinitePolicy::Propagate && self.non_finite_count != 0;

        let (sum, mean, variance, min, max) = if self.included_count == 0 {
            (None, None, None, None, None)
        } else if propagated_non_finite {
            let sum =
                if self.first_nan.is_some() || (self.positive_infinity && self.negative_infinity) {
                    f64::NAN
                } else if self.positive_infinity {
                    f64::INFINITY
                } else if self.negative_infinity {
                    f64::NEG_INFINITY
                } else {
                    self.finite_sum
                };
            let mean = sum / self.included_count as f64;
            let variance =
                variance_exists(self.options.variance, self.included_count).then_some(f64::NAN);
            let (min, max) = self.first_nan.map_or((self.min, self.max), |nan| {
                // NaN is unordered, so propagation makes both extrema NaN and
                // records the first row that caused the indeterminate result.
                (Some(nan), Some(nan))
            });
            (Some(sum), Some(mean), variance, min, max)
        } else {
            let denominator = match self.options.variance {
                VarianceNormalization::Population => self.included_count,
                VarianceNormalization::Sample => self.included_count.saturating_sub(1),
            };
            let variance =
                (denominator != 0).then_some(self.squared_deviations / denominator as f64);
            (
                Some(self.finite_sum),
                Some(self.mean),
                variance,
                self.min,
                self.max,
            )
        };

        ColumnSummary {
            selected_count: self.selected_count,
            finite_count: self.finite_count,
            non_finite_count: self.non_finite_count,
            included_count: self.included_count,
            zero_count: self.zero_count,
            nonzero_count: self.nonzero_count,
            sum,
            mean,
            variance,
            min,
            max,
        }
    }
}

fn variance_exists(normalization: VarianceNormalization, count: usize) -> bool {
    match normalization {
        VarianceNormalization::Population => count != 0,
        VarianceNormalization::Sample => count >= 2,
    }
}

/// Update one extremum without replacing the first row on a tie.
fn update_extrema(current: &mut Option<ColumnExtremum>, candidate: ColumnExtremum, minimum: bool) {
    let replace = match current {
        None => true,
        Some(old) if minimum => candidate.value < old.value,
        Some(old) => candidate.value > old.value,
    };
    if replace {
        *current = Some(candidate);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::column::{ColumnData, ColumnRole, DataColumn};
    use crate::dataset::DatasetKind;
    use crate::domain::{Domain, SurfaceDomain, VolumeDomain};

    fn surface(samples: usize) -> Domain {
        Domain::Surface(SurfaceDomain::new(None, samples).unwrap())
    }

    fn numbers(values: Vec<f64>) -> DataColumn {
        DataColumn::new("numbers", ColumnRole::Generic, ColumnData::Float64(values)).unwrap()
    }

    fn options(
        non_finite: NonFinitePolicy,
        variance: VarianceNormalization,
    ) -> ColumnSummaryOptions {
        ColumnSummaryOptions::new(non_finite, variance)
    }

    #[test]
    fn dense_population_summary_reports_values_counts_and_first_extrema() {
        let dataset = Dataset::dense(
            DatasetKind::Scalar,
            surface(5),
            vec![numbers(vec![2.0, -1.0, 0.0, 5.0, -1.0])],
        )
        .unwrap();
        let summary = summarize_column(
            &dataset,
            0,
            None,
            options(NonFinitePolicy::Skip, VarianceNormalization::Population),
        )
        .unwrap();

        assert_eq!(summary.selected_count, 5);
        assert_eq!(summary.finite_count, 5);
        assert_eq!(summary.non_finite_count, 0);
        assert_eq!(summary.included_count, 5);
        assert_eq!(summary.zero_count, 1);
        assert_eq!(summary.nonzero_count, 4);
        assert_eq!(summary.sum, Some(5.0));
        assert_eq!(summary.mean, Some(1.0));
        assert!((summary.variance.unwrap() - 5.2).abs() < 1.0e-12);
        assert!((summary.standard_deviation().unwrap() - 5.2_f64.sqrt()).abs() < 1.0e-12);
        assert_eq!(
            summary.min,
            Some(ColumnExtremum {
                value: -1.0,
                row: 1,
                sample: 1,
            })
        );
        assert_eq!(
            summary.max,
            Some(ColumnExtremum {
                value: 5.0,
                row: 3,
                sample: 3,
            })
        );
    }

    #[test]
    fn sample_variance_and_empty_selection_are_explicit() {
        let dataset = Dataset::dense(
            DatasetKind::Scalar,
            surface(3),
            vec![numbers(vec![1.0, 2.0, 3.0])],
        )
        .unwrap();
        let sample = summarize_column(
            &dataset,
            0,
            None,
            options(NonFinitePolicy::Skip, VarianceNormalization::Sample),
        )
        .unwrap();
        assert_eq!(sample.variance, Some(1.0));

        let none = SampleMask::new(surface(3), vec![false; 3]).unwrap();
        let empty = summarize_column(
            &dataset,
            0,
            Some(&none),
            options(NonFinitePolicy::Skip, VarianceNormalization::Sample),
        )
        .unwrap();
        assert_eq!(empty.selected_count, 0);
        assert_eq!(empty.sum, None);
        assert_eq!(empty.mean, None);
        assert_eq!(empty.variance, None);
        assert_eq!(empty.min, None);
        assert_eq!(empty.max, None);
    }

    #[test]
    fn sparse_masking_uses_domain_samples_and_reports_both_indices() {
        let domain = surface(6);
        let dataset = Dataset::indexed(
            DatasetKind::Scalar,
            domain.clone(),
            vec![4, 1, 5],
            vec![numbers(vec![40.0, 10.0, 50.0])],
        )
        .unwrap();
        // Select domain samples 1, 2, and 4. Sample 2 has no stored row and
        // therefore contributes no value to the summary.
        let mask = SampleMask::new(domain, vec![false, true, true, false, true, false]).unwrap();
        let summary = summarize_column(
            &dataset,
            0,
            Some(&mask),
            options(NonFinitePolicy::Skip, VarianceNormalization::Population),
        )
        .unwrap();

        assert_eq!(summary.selected_count, 2);
        assert_eq!(summary.mean, Some(25.0));
        assert_eq!(summary.variance, Some(225.0));
        assert_eq!(
            summary.min,
            Some(ColumnExtremum {
                value: 10.0,
                row: 1,
                sample: 1,
            })
        );
        assert_eq!(summary.max.unwrap().sample, 4);
    }

    #[test]
    fn non_finite_policies_skip_reject_or_propagate() {
        let dataset = Dataset::dense(
            DatasetKind::Scalar,
            surface(4),
            vec![numbers(vec![1.0, f64::NAN, f64::INFINITY, 3.0])],
        )
        .unwrap();
        let skipped = summarize_column(
            &dataset,
            0,
            None,
            options(NonFinitePolicy::Skip, VarianceNormalization::Population),
        )
        .unwrap();
        assert_eq!(skipped.selected_count, 4);
        assert_eq!(skipped.finite_count, 2);
        assert_eq!(skipped.non_finite_count, 2);
        assert_eq!(skipped.included_count, 2);
        assert_eq!(skipped.mean, Some(2.0));
        assert_eq!(skipped.variance, Some(1.0));

        assert!(matches!(
            summarize_column(
                &dataset,
                0,
                None,
                options(
                    NonFinitePolicy::Reject,
                    VarianceNormalization::Population
                )
            ),
            Err(Error::NonFinite { ref what, value })
                if what == "column 'numbers' row 1" && value.is_nan()
        ));

        let propagated = summarize_column(
            &dataset,
            0,
            None,
            options(
                NonFinitePolicy::Propagate,
                VarianceNormalization::Population,
            ),
        )
        .unwrap();
        assert_eq!(propagated.included_count, 4);
        assert!(propagated.sum.unwrap().is_nan());
        assert!(propagated.mean.unwrap().is_nan());
        assert!(propagated.variance.unwrap().is_nan());
        assert!(propagated.min.unwrap().value.is_nan());
        assert_eq!(propagated.min.unwrap().row, 1);
        assert!(propagated.max.unwrap().value.is_nan());
    }

    #[test]
    fn text_bad_indices_and_wrong_domains_are_rejected_before_reduction() {
        let text = DataColumn::new(
            "words",
            ColumnRole::Generic,
            ColumnData::Text(vec!["one".into(), "two".into()]),
        )
        .unwrap();
        let dataset = Dataset::dense(DatasetKind::Scalar, surface(2), vec![text]).unwrap();
        let opts = options(NonFinitePolicy::Skip, VarianceNormalization::Population);
        assert!(matches!(
            summarize_column(&dataset, 0, None, opts),
            Err(Error::InvalidParameter { ref name, .. }) if name == "summary column"
        ));
        assert!(matches!(
            summarize_column(&dataset, 1, None, opts),
            Err(Error::ColumnIndexOutOfRange { index: 1, len: 1 })
        ));

        let volume = Domain::Volume(VolumeDomain::new(None, [2, 1, 1], None).unwrap());
        let wrong = SampleMask::new(volume, vec![true, true]).unwrap();
        assert!(matches!(
            summarize_column(&dataset, 0, Some(&wrong), opts),
            Err(Error::DomainMismatch { .. })
        ));
    }

    #[test]
    fn finite_variance_is_stable_around_a_large_offset() {
        let dataset = Dataset::dense(
            DatasetKind::Scalar,
            surface(3),
            vec![numbers(vec![1.0e12 - 1.0, 1.0e12, 1.0e12 + 1.0])],
        )
        .unwrap();
        let summary = summarize_column(
            &dataset,
            0,
            None,
            options(NonFinitePolicy::Skip, VarianceNormalization::Population),
        )
        .unwrap();
        assert_eq!(summary.mean, Some(1.0e12));
        assert!((summary.variance.unwrap() - 2.0 / 3.0).abs() < 1.0e-12);
    }
}
