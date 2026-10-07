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
// Reusable processing loops over Dataset rows and time series. These are the
// file-neutral counterparts of the repeated loops in AFNI C programs that:
//
// * fetch several sub-brick values at one voxel;
// * map a voxel's complete time series to one scalar; or
// * copy a time series into work memory, alter it, and write it back.
//
// WHY THIS IS SEPARATE FROM DATASET
//
// Dataset owns validated structure. TimeSeriesView owns the row-oriented view.
// SampleMask owns domain-sized selection. This module composes those pieces and
// owns only iteration policy: scratch-buffer reuse, sparse row/sample mapping,
// masked-out behavior, and atomic construction of validated outputs.
// ---------------------------------------------------------------------------

//! Checked voxel/node and time-series processing loops.

use crate::column::{ColumnData, ColumnRole, DataColumn, ValueMetadataPolicy};
use crate::dataset::{Dataset, DatasetKind};
use crate::error::{Error, Result};
use crate::mask::SampleMask;
use crate::numeric::NonFinitePolicy;
use crate::reduction::VarianceNormalization;

/// Location of one stored dataset row in both row and domain coordinates.
///
/// `row` indexes every [`DataColumn`] in the dataset. `sample` is the volume
/// voxel or surface-node index represented by that row. They are equal for a
/// dense dataset and may differ for a sparse dataset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SampleContext {
    /// Zero-based stored dataset row.
    pub row: usize,
    /// Zero-based sample index in the dataset's complete [`Domain`](crate::domain::Domain).
    pub sample: usize,
}

/// What a masked time-series transformation writes for unselected samples.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UnselectedSeries {
    /// Copy the original time series unchanged.
    Preserve,
    /// Write the same value at every time point.
    ///
    /// Zero is the common AFNI command-line convention; NaN is useful when
    /// masked-out samples should remain visibly missing.
    Fill(f64),
}

/// Description of one column produced by [`summarize_time_series_with`].
///
/// The label and role become ordinary [`DataColumn`] metadata. `unselected`
/// is written for a row excluded by the optional [`SampleMask`], so a
/// multi-output operation can choose a different outside-mask value for each
/// result.
#[derive(Debug, Clone, PartialEq)]
pub struct SummaryOutput {
    label: String,
    role: ColumnRole,
    unselected: f64,
}

impl SummaryOutput {
    /// Define one output column.
    ///
    /// An empty or whitespace-only label is rejected before any row callback
    /// runs.
    pub fn new(label: impl Into<String>, role: ColumnRole, unselected: f64) -> Result<Self> {
        let label = label.into();
        if label.trim().is_empty() {
            return Err(Error::Empty("summary output label".into()));
        }
        Ok(Self {
            label,
            role,
            unselected,
        })
    }

    /// Output column label.
    pub fn label(&self) -> &str {
        &self.label
    }

    /// Output column role.
    pub fn role(&self) -> &ColumnRole {
        &self.role
    }

    /// Value written for a row outside the mask.
    pub fn unselected_value(&self) -> f64 {
        self.unselected
    }
}

/// A built-in statistic computed across one voxel or node's time series.
///
/// These variants cover the common streaming summaries that do not need to
/// retain or sort the complete series. Use [`summarize_time_series_with`] for a
/// specialized statistic or for several results produced by one custom
/// calculation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TimeSeriesStatistic {
    /// Number of samples included after applying the non-finite policy.
    Count,
    /// Sum of included samples.
    Sum,
    /// Arithmetic mean of included samples.
    Mean,
    /// Sum of squared included samples.
    SumOfSquares,
    /// Euclidean length: the square root of the sum of squares.
    L2Norm,
    /// Variance with an explicit population or sample denominator.
    Variance(VarianceNormalization),
    /// Square root of the variance with an explicit denominator.
    StandardDeviation(VarianceNormalization),
    /// Square root of the mean squared value.
    RootMeanSquare,
    /// Least-squares slope per second, using the dataset's time step.
    ///
    /// This statistic requires [`Dataset::time_step_seconds`] to be present.
    /// Skipped non-finite samples retain their original temporal positions.
    Slope,
    /// Smallest included value.
    Minimum,
    /// Largest included value.
    Maximum,
}

impl TimeSeriesStatistic {
    /// Describe a generic output column for this statistic.
    ///
    /// `unselected` is written outside the optional spatial mask. Undefined
    /// statistics within a selected row—for example a mean after every sample
    /// was skipped—are represented by NaN instead.
    pub fn output(self, label: impl Into<String>, unselected: f64) -> Result<TimeSeriesSummary> {
        self.output_with_role(label, ColumnRole::Generic, unselected)
    }

    /// Describe an output column with an explicit semantic role.
    pub fn output_with_role(
        self,
        label: impl Into<String>,
        role: ColumnRole,
        unselected: f64,
    ) -> Result<TimeSeriesSummary> {
        TimeSeriesSummary::new(self, label, role, unselected)
    }
}

/// One requested output from [`summarize_time_series`].
#[derive(Debug, Clone, PartialEq)]
pub struct TimeSeriesSummary {
    statistic: TimeSeriesStatistic,
    output: SummaryOutput,
}

impl TimeSeriesSummary {
    /// Pair a built-in statistic with its output-column description.
    pub fn new(
        statistic: TimeSeriesStatistic,
        label: impl Into<String>,
        role: ColumnRole,
        unselected: f64,
    ) -> Result<Self> {
        Ok(Self {
            statistic,
            output: SummaryOutput::new(label, role, unselected)?,
        })
    }

    /// Statistic calculated for each selected time series.
    pub fn statistic(&self) -> TimeSeriesStatistic {
        self.statistic
    }

    /// Description of the resulting column.
    pub fn output_description(&self) -> &SummaryOutput {
        &self.output
    }
}

/// Combine one or more numeric columns into one derived `f64` column.
///
/// `column_indices` defines the order of the values passed to `transform` and
/// may repeat an index. All indices and numeric types are validated before the
/// first callback. The same scratch vector is reused for every selected row.
/// `mask`, when present, is checked against the dataset domain and interpreted
/// in domain-sample order; `transform` is not called for an unselected row and
/// `unselected_value` is written instead.
///
/// The output has one value per stored dataset row, so it can be passed directly
/// to [`Dataset::append_column`], [`Dataset::replace_column`], or
/// [`Dataset::replace_columns`].
pub fn combine_columns(
    dataset: &Dataset,
    column_indices: &[usize],
    mask: Option<&SampleMask>,
    unselected_value: f64,
    label: impl Into<String>,
    role: ColumnRole,
    mut transform: impl FnMut(SampleContext, &[f64]) -> Result<f64>,
) -> Result<DataColumn> {
    if column_indices.is_empty() {
        return Err(Error::Empty("voxelwise input columns".into()));
    }
    require_mask_domain(dataset, mask)?;

    // Resolve and validate every input before invoking user code, so a bad
    // later column cannot leave externally observed partial work.
    let columns = column_indices
        .iter()
        .map(|&index| dataset.column_at(index))
        .collect::<Result<Vec<_>>>()?;
    if let Some(column) = columns.iter().find(|column| !column.values().is_numeric()) {
        return Err(Error::InvalidParameter {
            name: "voxelwise input column".into(),
            reason: format!("column '{}' is text, not numeric", column.label()),
        });
    }

    let mut scratch = vec![0.0; columns.len()];
    let mut output = Vec::with_capacity(dataset.row_count());
    for row in 0..dataset.row_count() {
        let context = sample_context(dataset, row)?;
        if !sample_is_selected(mask, context.sample)? {
            output.push(unselected_value);
            continue;
        }
        for (value, column) in scratch.iter_mut().zip(&columns) {
            *value = column
                .values()
                .get_f64(row)
                .expect("validated numeric column and dataset row");
        }
        output.push(transform(context, &scratch)?);
    }
    DataColumn::new(label, role, ColumnData::Float64(output))
}

/// Reduce each selected voxel/node row across time into several output columns.
///
/// A dataset row is one spatial sample and its temporal values live in several
/// columns. This function presents those values as one reusable mutable slice
/// and invokes `summarize` once per selected row. Returning an array permits
/// several related results—such as low-band and total power—from one expensive
/// calculation. `N` must be greater than zero.
///
/// The input dataset must satisfy [`Dataset::time_series`]. Every output follows
/// the dataset's stored row order, including sparse maps. The optional mask is
/// checked before the first callback and interpreted in domain-sample order.
/// Masked-out rows do not invoke `summarize`; each receives the corresponding
/// [`SummaryOutput::unselected_value`]. Construction is atomic: an error
/// returns no output columns.
pub fn summarize_time_series_with<const N: usize>(
    dataset: &Dataset,
    mask: Option<&SampleMask>,
    outputs: [SummaryOutput; N],
    mut summarize: impl FnMut(SampleContext, &mut [f64]) -> Result<[f64; N]>,
) -> Result<[DataColumn; N]> {
    let columns = summarize_time_series_dynamic_with(
        dataset,
        mask,
        &outputs,
        |context, series, row_outputs| {
            row_outputs.copy_from_slice(&summarize(context, series)?);
            Ok(())
        },
    )?;
    columns
        .try_into()
        .map_err(|columns: Vec<DataColumn>| Error::LengthMismatch {
            what: "time-series summary output columns".into(),
            expected: N,
            found: columns.len(),
        })
}

/// Dynamically reduce each selected voxel/node time series into output columns.
///
/// This is the runtime-sized counterpart to [`summarize_time_series_with`]. It
/// is intended for command-line programs whose requested outputs are assembled
/// from user options. `outputs.len()` determines the length of `row_outputs`
/// passed to the callback, and the same series and output buffers are reused
/// for every selected row. The callback should write each requested result into
/// `row_outputs`; entries start as NaN for every selected row.
///
/// Masked-out rows do not invoke `summarize` and instead receive each
/// [`SummaryOutput::unselected_value`]. The returned columns follow `outputs`
/// order. An empty output list is rejected before any callback is invoked.
pub fn summarize_time_series_dynamic_with(
    dataset: &Dataset,
    mask: Option<&SampleMask>,
    outputs: &[SummaryOutput],
    mut summarize: impl FnMut(SampleContext, &mut [f64], &mut [f64]) -> Result<()>,
) -> Result<Vec<DataColumn>> {
    if outputs.is_empty() {
        return Err(Error::Empty("time-series summary outputs".into()));
    }
    let time = dataset.time_series()?;
    require_mask_domain(dataset, mask)?;

    let mut values: Vec<Vec<f64>> = (0..outputs.len())
        .map(|_| Vec::with_capacity(time.series_count()))
        .collect();
    let mut series = vec![0.0; time.time_point_count()];
    let mut row_outputs = vec![f64::NAN; outputs.len()];
    for row in 0..time.series_count() {
        let context = sample_context(dataset, row)?;
        if sample_is_selected(mask, context.sample)? {
            time.copy_series_into(row, &mut series)?;
            row_outputs.fill(f64::NAN);
            summarize(context, &mut series, &mut row_outputs)?;
        } else {
            for (value, output) in row_outputs.iter_mut().zip(outputs) {
                *value = output.unselected;
            }
        }
        for (column, &value) in values.iter_mut().zip(&row_outputs) {
            column.push(value);
        }
    }

    outputs
        .iter()
        .zip(values)
        .map(|(output, values)| {
            DataColumn::new(
                output.label.clone(),
                output.role.clone(),
                ColumnData::Float64(values),
            )
        })
        .collect()
}

/// Compute common statistics across every voxel or node's time series.
///
/// All requested statistics are accumulated together in one pass over each
/// selected series. The returned scalar dataset is derived from `dataset`: it
/// preserves the domain and dense-or-sparse row mapping, clears the time axis,
/// and contains one `Float64` column per request in request order.
///
/// `non_finite` applies before every statistic. `Skip` omits NaN and infinity,
/// `Reject` reports the first offending row and time point, and `Propagate`
/// includes them using IEEE-style results. An undefined selected-row result is
/// NaN; values outside `mask` use the per-output `unselected` value.
pub fn summarize_time_series<const N: usize>(
    dataset: &Dataset,
    mask: Option<&SampleMask>,
    summaries: [TimeSeriesSummary; N],
    non_finite: NonFinitePolicy,
) -> Result<Dataset> {
    summarize_time_series_dynamic(dataset, mask, &summaries, non_finite)
}

/// Compute a runtime-selected list of common time-series statistics.
///
/// This is the dynamic counterpart to [`summarize_time_series`]. It is useful
/// for AFNI-style programs in which flags such as `-mean`, `-stdev`, and
/// `-slope` build a `Vec<TimeSeriesSummary>` at runtime. All requested
/// statistics are still accumulated together in one traversal of each spatial
/// series, and output columns preserve request order.
pub fn summarize_time_series_dynamic(
    dataset: &Dataset,
    mask: Option<&SampleMask>,
    summaries: &[TimeSeriesSummary],
    non_finite: NonFinitePolicy,
) -> Result<Dataset> {
    if summaries.is_empty() {
        return Err(Error::Empty("time-series summaries".into()));
    }
    let statistics: Vec<TimeSeriesStatistic> =
        summaries.iter().map(TimeSeriesSummary::statistic).collect();
    let time_step = if statistics.contains(&TimeSeriesStatistic::Slope) {
        Some(
            dataset
                .time_step_seconds()
                .ok_or_else(|| Error::InvalidParameter {
                    name: "time-series slope".into(),
                    reason: "requires a dataset time step in seconds".into(),
                })?,
        )
    } else {
        None
    };
    let outputs: Vec<SummaryOutput> = summaries
        .iter()
        .map(|summary| summary.output.clone())
        .collect();
    let columns = summarize_time_series_dynamic_with(
        dataset,
        mask,
        &outputs,
        |context, series, row_outputs| {
            let summary = TemporalAccumulator::from_series(context, series, non_finite)?;
            for (output, &statistic) in row_outputs.iter_mut().zip(&statistics) {
                *output = summary.value(statistic, time_step);
            }
            Ok(())
        },
    )?;
    dataset.derive(DatasetKind::Scalar, columns)
}

/// Transform each selected spatial time series and return a validated dataset.
///
/// Selected temporal columns become `Float64`; non-temporal columns retain their
/// original type and values. Every selected row is copied into one reusable
/// buffer, passed to `transform` as a mutable slice of fixed length, and then
/// transposed back into temporal columns. Because the callback cannot resize
/// the slice, the number of time points cannot drift.
///
/// `unselected` explicitly chooses whether masked-out rows retain their input
/// or receive a fill value. `metadata` independently chooses whether each
/// transformed temporal column retains value-dependent metadata. Construction
/// is atomic: a callback error returns no new dataset and leaves `dataset`
/// unchanged.
pub fn transform_time_series(
    dataset: &Dataset,
    mask: Option<&SampleMask>,
    unselected: UnselectedSeries,
    metadata: ValueMetadataPolicy,
    mut transform: impl FnMut(SampleContext, &mut [f64]) -> Result<()>,
) -> Result<Dataset> {
    let time = dataset.time_series()?;
    require_mask_domain(dataset, mask)?;
    let time_columns: Vec<usize> = time.column_indices().collect();

    // Each inner vector becomes one temporal output column. Reserving rows once
    // avoids repeated allocation while keeping the callback's input contiguous.
    let mut output_values: Vec<Vec<f64>> = (0..time.time_point_count())
        .map(|_| Vec::with_capacity(time.series_count()))
        .collect();
    let mut scratch = vec![0.0; time.time_point_count()];

    for row in 0..time.series_count() {
        let context = sample_context(dataset, row)?;
        if sample_is_selected(mask, context.sample)? {
            time.copy_series_into(row, &mut scratch)?;
            transform(context, &mut scratch)?;
        } else {
            match unselected {
                UnselectedSeries::Preserve => time.copy_series_into(row, &mut scratch)?,
                UnselectedSeries::Fill(value) => scratch.fill(value),
            }
        }
        for (column, &value) in output_values.iter_mut().zip(&scratch) {
            column.push(value);
        }
    }

    let mut columns = dataset.columns().to_vec();
    for (index, values) in time_columns.into_iter().zip(output_values) {
        columns[index] = columns[index].with_values(ColumnData::Float64(values), metadata)?;
    }
    dataset.replace_columns(columns)
}

/// Validate an optional mask before processing begins.
fn require_mask_domain(dataset: &Dataset, mask: Option<&SampleMask>) -> Result<()> {
    if let Some(mask) = mask {
        mask.require_domain(dataset.domain())?;
    }
    Ok(())
}

/// Convert a valid dataset row to the public dual-index context.
fn sample_context(dataset: &Dataset, row: usize) -> Result<SampleContext> {
    let sample = dataset
        .sample_for_row(row)
        .ok_or_else(|| Error::InvalidParameter {
            name: "dataset map".into(),
            reason: format!("row {row} has no domain sample"),
        })? as usize;
    Ok(SampleContext { row, sample })
}

/// Test selection after domain compatibility has been established.
fn sample_is_selected(mask: Option<&SampleMask>, sample: usize) -> Result<bool> {
    mask.map_or(Ok(true), |mask| mask.is_selected(sample))
}

/// Streaming state shared by all built-in summaries for one temporal row.
#[derive(Debug, Clone, Copy)]
struct TemporalAccumulator {
    included_count: usize,
    finite_count: usize,
    finite_sum: f64,
    sum_of_squares: f64,
    mean: f64,
    squared_deviations: f64,
    mean_time_point: f64,
    time_point_deviations: f64,
    time_value_deviations: f64,
    minimum: Option<f64>,
    maximum: Option<f64>,
    saw_nan: bool,
    saw_positive_infinity: bool,
    saw_negative_infinity: bool,
    propagated_non_finite: bool,
}

impl TemporalAccumulator {
    fn from_series(
        context: SampleContext,
        series: &[f64],
        non_finite: NonFinitePolicy,
    ) -> Result<Self> {
        let mut summary = Self {
            included_count: 0,
            finite_count: 0,
            finite_sum: 0.0,
            sum_of_squares: 0.0,
            mean: 0.0,
            squared_deviations: 0.0,
            mean_time_point: 0.0,
            time_point_deviations: 0.0,
            time_value_deviations: 0.0,
            minimum: None,
            maximum: None,
            saw_nan: false,
            saw_positive_infinity: false,
            saw_negative_infinity: false,
            propagated_non_finite: false,
        };

        for (time_point, &value) in series.iter().enumerate() {
            if !value.is_finite() {
                match non_finite {
                    NonFinitePolicy::Reject => {
                        return Err(Error::NonFinite {
                            what: format!(
                                "time-series row {} (sample {}) time point {time_point}",
                                context.row, context.sample
                            ),
                            value,
                        });
                    }
                    NonFinitePolicy::Skip => continue,
                    NonFinitePolicy::Propagate => {
                        summary.included_count += 1;
                        summary.propagated_non_finite = true;
                        if value.is_nan() {
                            summary.saw_nan = true;
                        } else {
                            summary.update_extrema(value);
                            if value.is_sign_positive() {
                                summary.saw_positive_infinity = true;
                            } else {
                                summary.saw_negative_infinity = true;
                            }
                        }
                        continue;
                    }
                }
            }

            summary.included_count += 1;
            summary.finite_count += 1;
            summary.finite_sum += value;
            summary.sum_of_squares += value * value;
            summary.update_extrema(value);

            // Welford's update avoids the cancellation in mean(x*x)-mean(x)^2.
            let count = summary.finite_count as f64;
            let time_point = time_point as f64;
            let value_delta = value - summary.mean;
            let time_delta = time_point - summary.mean_time_point;
            summary.mean += value_delta / count;
            summary.mean_time_point += time_delta / count;
            summary.squared_deviations += value_delta * (value - summary.mean);
            summary.time_point_deviations += time_delta * (time_point - summary.mean_time_point);
            summary.time_value_deviations += time_delta * (value - summary.mean);
        }
        Ok(summary)
    }

    fn update_extrema(&mut self, value: f64) {
        if self.minimum.map_or(true, |current| value < current) {
            self.minimum = Some(value);
        }
        if self.maximum.map_or(true, |current| value > current) {
            self.maximum = Some(value);
        }
    }

    fn value(self, statistic: TimeSeriesStatistic, time_step: Option<f64>) -> f64 {
        match statistic {
            TimeSeriesStatistic::Count => self.included_count as f64,
            TimeSeriesStatistic::Sum => self.sum(),
            TimeSeriesStatistic::Mean => {
                if self.included_count == 0 {
                    f64::NAN
                } else {
                    self.sum() / self.included_count as f64
                }
            }
            TimeSeriesStatistic::SumOfSquares => self.sum_of_squares(),
            TimeSeriesStatistic::L2Norm => self.sum_of_squares().sqrt(),
            TimeSeriesStatistic::Variance(normalization) => self.variance(normalization),
            TimeSeriesStatistic::StandardDeviation(normalization) => {
                self.variance(normalization).sqrt()
            }
            TimeSeriesStatistic::RootMeanSquare => {
                if self.included_count == 0 || self.saw_nan {
                    f64::NAN
                } else if self.saw_positive_infinity || self.saw_negative_infinity {
                    f64::INFINITY
                } else {
                    (self.sum_of_squares / self.included_count as f64).sqrt()
                }
            }
            TimeSeriesStatistic::Slope => self.slope(
                time_step.expect("slope requests validate the dataset time step before rows"),
            ),
            TimeSeriesStatistic::Minimum => self.extremum(self.minimum),
            TimeSeriesStatistic::Maximum => self.extremum(self.maximum),
        }
    }

    fn sum(self) -> f64 {
        if self.included_count == 0
            || self.saw_nan
            || (self.saw_positive_infinity && self.saw_negative_infinity)
        {
            f64::NAN
        } else if self.saw_positive_infinity {
            f64::INFINITY
        } else if self.saw_negative_infinity {
            f64::NEG_INFINITY
        } else {
            self.finite_sum
        }
    }

    fn variance(self, normalization: VarianceNormalization) -> f64 {
        if self.propagated_non_finite {
            return f64::NAN;
        }
        match normalization {
            VarianceNormalization::Population if self.included_count != 0 => {
                self.squared_deviations / self.included_count as f64
            }
            VarianceNormalization::Sample if self.included_count > 1 => {
                self.squared_deviations / (self.included_count - 1) as f64
            }
            VarianceNormalization::Population | VarianceNormalization::Sample => f64::NAN,
        }
    }

    fn sum_of_squares(self) -> f64 {
        if self.included_count == 0 || self.saw_nan {
            f64::NAN
        } else if self.saw_positive_infinity || self.saw_negative_infinity {
            f64::INFINITY
        } else {
            self.sum_of_squares
        }
    }

    fn slope(self, time_step: f64) -> f64 {
        if self.propagated_non_finite || self.finite_count < 2 || self.time_point_deviations == 0.0
        {
            f64::NAN
        } else {
            self.time_value_deviations / self.time_point_deviations / time_step
        }
    }

    fn extremum(self, value: Option<f64>) -> f64 {
        if self.saw_nan {
            f64::NAN
        } else {
            value.unwrap_or(f64::NAN)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::column::{ColumnRange, RecordedRange};
    use crate::dataset::DatasetKind;
    use crate::domain::{Domain, SurfaceDomain, VolumeDomain};

    fn surface(samples: usize) -> Domain {
        Domain::Surface(SurfaceDomain::new(None, samples).unwrap())
    }

    fn column(label: &str, role: ColumnRole, values: Vec<f64>) -> DataColumn {
        DataColumn::new(label, role, ColumnData::Float64(values)).unwrap()
    }

    #[test]
    fn row_mapping_uses_numeric_columns_sparse_samples_and_a_domain_mask() {
        let dataset = Dataset::indexed(
            DatasetKind::Scalar,
            surface(6),
            vec![4, 1, 3],
            vec![
                column("a", ColumnRole::Intensity, vec![40.0, 10.0, 30.0]),
                DataColumn::new("b", ColumnRole::Threshold, ColumnData::Int32(vec![4, 1, 3]))
                    .unwrap(),
            ],
        )
        .unwrap();
        let mask = SampleMask::new(
            dataset.domain().clone(),
            vec![false, true, false, false, true, false],
        )
        .unwrap();
        let mut seen = Vec::new();

        let output = combine_columns(
            &dataset,
            &[0, 1, 0],
            Some(&mask),
            f64::NAN,
            "combined",
            ColumnRole::Generic,
            |context, values| {
                seen.push(context);
                Ok(values[0] + values[1] + values[2])
            },
        )
        .unwrap();

        assert_eq!(
            seen,
            [
                SampleContext { row: 0, sample: 4 },
                SampleContext { row: 1, sample: 1 }
            ]
        );
        let ColumnData::Float64(values) = output.values() else {
            panic!("row mapping must produce Float64")
        };
        assert_eq!(values[0], 84.0);
        assert_eq!(values[1], 21.0);
        assert!(values[2].is_nan());
    }

    #[test]
    fn row_mapping_validates_every_input_before_calling_the_callback() {
        let dataset = Dataset::dense(
            DatasetKind::Scalar,
            surface(2),
            vec![
                column("numbers", ColumnRole::Generic, vec![1.0, 2.0]),
                DataColumn::new(
                    "words",
                    ColumnRole::Generic,
                    ColumnData::Text(vec!["a".into(), "b".into()]),
                )
                .unwrap(),
            ],
        )
        .unwrap();
        let mut calls = 0;
        let result = combine_columns(
            &dataset,
            &[0, 1],
            None,
            0.0,
            "bad",
            ColumnRole::Generic,
            |_, _| {
                calls += 1;
                Ok(0.0)
            },
        );
        assert!(matches!(
            result,
            Err(Error::InvalidParameter { ref name, .. }) if name == "voxelwise input column"
        ));
        assert_eq!(calls, 0);
        assert!(combine_columns(
            &dataset,
            &[],
            None,
            0.0,
            "empty",
            ColumnRole::Generic,
            |_, _| Ok(0.0)
        )
        .is_err());
        assert!(matches!(
            combine_columns(
                &dataset,
                &[2],
                None,
                0.0,
                "bad index",
                ColumnRole::Generic,
                |_, _| Ok(0.0)
            ),
            Err(Error::ColumnIndexOutOfRange { index: 2, len: 2 })
        ));
    }

    #[test]
    fn one_output_custom_summary_reuses_row_order_and_skips_masked_callbacks() {
        let dataset = Dataset::indexed(
            DatasetKind::TimeSeries,
            surface(8),
            vec![6, 2],
            vec![
                column("t0", ColumnRole::TimePoint, vec![1.0, 10.0]),
                column("t1", ColumnRole::TimePoint, vec![3.0, 14.0]),
            ],
        )
        .unwrap();
        let mask = SampleMask::new(
            dataset.domain().clone(),
            vec![false, false, true, false, false, false, false, false],
        )
        .unwrap();
        let mut seen = Vec::new();
        let [means] = summarize_time_series_with(
            &dataset,
            Some(&mask),
            [SummaryOutput::new("mean", ColumnRole::Generic, -1.0).unwrap()],
            |context, series| {
                seen.push(context);
                Ok([series.iter().sum::<f64>() / series.len() as f64])
            },
        )
        .unwrap();

        assert_eq!(means.values(), &ColumnData::Float64(vec![-1.0, 12.0]));
        assert_eq!(seen, [SampleContext { row: 1, sample: 2 }]);
    }

    #[test]
    fn row_summaries_return_several_columns_from_one_mutable_series_callback() {
        let dataset = Dataset::indexed(
            DatasetKind::TimeSeries,
            surface(8),
            vec![6, 2],
            vec![
                column("t0", ColumnRole::TimePoint, vec![1.0, 10.0]),
                column("t1", ColumnRole::TimePoint, vec![3.0, 14.0]),
            ],
        )
        .unwrap();
        let mask = SampleMask::new(
            dataset.domain().clone(),
            vec![false, false, true, false, false, false, false, false],
        )
        .unwrap();
        let sum = SummaryOutput::new("sum", ColumnRole::Generic, -1.0).unwrap();
        let spread = SummaryOutput::new("spread", ColumnRole::Generic, -2.0).unwrap();
        assert_eq!(sum.label(), "sum");
        assert_eq!(sum.role(), &ColumnRole::Generic);
        assert_eq!(sum.unselected_value(), -1.0);

        let mut calls = 0;
        let [sums, spreads] =
            summarize_time_series_with(&dataset, Some(&mask), [sum, spread], |context, series| {
                calls += 1;
                assert_eq!(context, SampleContext { row: 1, sample: 2 });
                let sum = series.iter().sum();
                series[0] = 100.0;
                Ok([sum, series[0] - series[1]])
            })
            .unwrap();

        assert_eq!(calls, 1);
        assert_eq!(sums.values(), &ColumnData::Float64(vec![-1.0, 24.0]));
        assert_eq!(spreads.values(), &ColumnData::Float64(vec![-2.0, 86.0]));
        // The reusable scratch buffer is private to the operation; input data
        // is never changed by a mutating summary callback.
        assert_eq!(dataset.columns()[0].values().get_f64(1), Some(10.0));

        let empty: [SummaryOutput; 0] = [];
        assert!(summarize_time_series_with(&dataset, None, empty, |_, _| Ok([])).is_err());
        assert!(SummaryOutput::new(" ", ColumnRole::Generic, 0.0).is_err());
    }

    #[test]
    fn dynamic_custom_summaries_reuse_runtime_sized_buffers() {
        let dataset = Dataset::dense(
            DatasetKind::TimeSeries,
            surface(2),
            vec![
                column("t0", ColumnRole::TimePoint, vec![1.0, 10.0]),
                column("t1", ColumnRole::TimePoint, vec![3.0, 14.0]),
            ],
        )
        .unwrap();
        let outputs = vec![
            SummaryOutput::new("sum", ColumnRole::Generic, -1.0).unwrap(),
            SummaryOutput::new("spread", ColumnRole::Generic, -2.0).unwrap(),
            SummaryOutput::new("optional", ColumnRole::Generic, -3.0).unwrap(),
        ];
        let mut series_pointer = None;
        let mut output_pointer = None;
        let mut calls = 0;

        let columns = summarize_time_series_dynamic_with(
            &dataset,
            None,
            &outputs,
            |_, series, row_outputs| {
                calls += 1;
                assert_eq!(row_outputs.len(), outputs.len());
                match (series_pointer, output_pointer) {
                    (Some(previous_series), Some(previous_output)) => {
                        assert_eq!(previous_series, series.as_ptr());
                        assert_eq!(previous_output, row_outputs.as_ptr());
                    }
                    (None, None) => {
                        series_pointer = Some(series.as_ptr());
                        output_pointer = Some(row_outputs.as_ptr());
                    }
                    _ => unreachable!(),
                }
                row_outputs[0] = series.iter().sum();
                row_outputs[1] = series[1] - series[0];
                // An intentionally unwritten output remains NaN for this row,
                // rather than retaining the preceding row's value.
                Ok(())
            },
        )
        .unwrap();

        assert_eq!(calls, 2);
        assert_eq!(columns.len(), 3);
        assert_eq!(columns[0].values(), &ColumnData::Float64(vec![4.0, 24.0]));
        assert_eq!(columns[1].values(), &ColumnData::Float64(vec![2.0, 4.0]));
        assert!(columns[2].values().get_f64(0).unwrap().is_nan());
        assert!(columns[2].values().get_f64(1).unwrap().is_nan());

        let mut empty_calls = 0;
        assert!(
            summarize_time_series_dynamic_with(&dataset, None, &[], |_, _, _| {
                empty_calls += 1;
                Ok(())
            })
            .is_err()
        );
        assert_eq!(empty_calls, 0);
    }

    #[test]
    fn dynamic_built_in_summaries_accept_runtime_requests() {
        let dataset = Dataset::dense(
            DatasetKind::TimeSeries,
            surface(2),
            vec![
                column("t0", ColumnRole::TimePoint, vec![1.0, 10.0]),
                column("t1", ColumnRole::TimePoint, vec![3.0, 14.0]),
            ],
        )
        .unwrap();
        let requests = vec![
            TimeSeriesStatistic::Mean.output("mean", 0.0).unwrap(),
            TimeSeriesStatistic::Maximum.output("maximum", 0.0).unwrap(),
        ];

        let dynamic = dataset
            .summarize_time_series_dynamic(None, &requests, NonFinitePolicy::Skip)
            .unwrap();
        let fixed = dataset
            .summarize_time_series(
                None,
                [requests[0].clone(), requests[1].clone()],
                NonFinitePolicy::Skip,
            )
            .unwrap();

        assert!(dynamic.eq_nan_aware(&fixed));
        assert_eq!(dynamic.columns()[0].label(), "mean");
        assert_eq!(dynamic.columns()[1].label(), "maximum");
        assert!(dataset
            .summarize_time_series_dynamic(None, &[], NonFinitePolicy::Skip)
            .is_err());
    }

    #[test]
    fn built_in_time_series_statistics_share_one_checked_pass() {
        let dataset = Dataset::dense(
            DatasetKind::TimeSeries,
            surface(3),
            vec![
                column("t0", ColumnRole::TimePoint, vec![1.0, 10.0, 100.0]),
                column("t1", ColumnRole::TimePoint, vec![2.0, f64::NAN, 200.0]),
                column("t2", ColumnRole::TimePoint, vec![3.0, 30.0, 300.0]),
                column("t3", ColumnRole::TimePoint, vec![4.0, f64::INFINITY, 400.0]),
            ],
        )
        .unwrap()
        .with_time_step_seconds(Some(2.0))
        .unwrap();
        let mask = SampleMask::new(dataset.domain().clone(), vec![true, true, false]).unwrap();
        let requests = [
            TimeSeriesStatistic::Count.output("count", -1.0).unwrap(),
            TimeSeriesStatistic::Sum.output("sum", -2.0).unwrap(),
            TimeSeriesStatistic::Mean.output("mean", -3.0).unwrap(),
            TimeSeriesStatistic::SumOfSquares
                .output("sum_of_squares", -3.25)
                .unwrap(),
            TimeSeriesStatistic::L2Norm.output("l2_norm", -3.5).unwrap(),
            TimeSeriesStatistic::Variance(VarianceNormalization::Population)
                .output("variance", -4.0)
                .unwrap(),
            TimeSeriesStatistic::StandardDeviation(VarianceNormalization::Sample)
                .output("stdev", -5.0)
                .unwrap(),
            TimeSeriesStatistic::RootMeanSquare
                .output("rms", -6.0)
                .unwrap(),
            TimeSeriesStatistic::Slope.output("slope", -6.5).unwrap(),
            TimeSeriesStatistic::Minimum
                .output("minimum", -7.0)
                .unwrap(),
            TimeSeriesStatistic::Maximum
                .output("maximum", -8.0)
                .unwrap(),
        ];
        assert_eq!(requests[0].statistic(), TimeSeriesStatistic::Count);
        assert_eq!(requests[0].output_description().label(), "count");

        let summary = dataset
            .summarize_time_series(Some(&mask), requests, NonFinitePolicy::Skip)
            .unwrap();

        assert_eq!(summary.kind(), &DatasetKind::Scalar);
        assert_eq!(summary.domain(), dataset.domain());
        assert_eq!(summary.map(), dataset.map());
        assert_eq!(summary.time_step_seconds(), None);
        let expected = [
            [4.0, 2.0, -1.0],
            [10.0, 40.0, -2.0],
            [2.5, 20.0, -3.0],
            [30.0, 1_000.0, -3.25],
            [30.0_f64.sqrt(), 1_000.0_f64.sqrt(), -3.5],
            [1.25, 100.0, -4.0],
            [(5.0_f64 / 3.0).sqrt(), 200.0_f64.sqrt(), -5.0],
            [7.5_f64.sqrt(), 500.0_f64.sqrt(), -6.0],
            [0.5, 5.0, -6.5],
            [1.0, 10.0, -7.0],
            [4.0, 30.0, -8.0],
        ];
        for (column, expected) in summary.columns().iter().zip(expected) {
            let ColumnData::Float64(values) = column.values() else {
                panic!("temporal statistics must produce Float64")
            };
            for (&actual, expected) in values.iter().zip(expected) {
                assert!((actual - expected).abs() < 1.0e-12);
            }
        }

        let rejected = summarize_time_series(
            &dataset,
            None,
            [TimeSeriesStatistic::Mean.output("mean", 0.0).unwrap()],
            NonFinitePolicy::Reject,
        );
        assert!(
            matches!(rejected, Err(Error::NonFinite { ref what, .. }) if what.contains("time point 1"))
        );

        let propagated = summarize_time_series(
            &dataset,
            None,
            [
                TimeSeriesStatistic::Count.output("count", 0.0).unwrap(),
                TimeSeriesStatistic::Sum.output("sum", 0.0).unwrap(),
                TimeSeriesStatistic::Minimum.output("minimum", 0.0).unwrap(),
            ],
            NonFinitePolicy::Propagate,
        )
        .unwrap();
        assert_eq!(
            propagated.columns()[0].values(),
            &ColumnData::Float64(vec![4.0, 4.0, 4.0])
        );
        assert!(propagated.columns()[1]
            .values()
            .get_f64(1)
            .unwrap()
            .is_nan());
        assert!(propagated.columns()[2]
            .values()
            .get_f64(1)
            .unwrap()
            .is_nan());

        let without_time_step = Dataset::dense(
            DatasetKind::TimeSeries,
            surface(1),
            vec![column("t0", ColumnRole::TimePoint, vec![1.0])],
        )
        .unwrap();
        assert!(summarize_time_series(
            &without_time_step,
            None,
            [TimeSeriesStatistic::Slope.output("slope", 0.0).unwrap()],
            NonFinitePolicy::Reject,
        )
        .is_err());
    }

    #[test]
    fn time_series_transform_preserves_or_fills_unselected_rows_and_other_columns() {
        let recorded = RecordedRange {
            range: ColumnRange::new(1.0, 20.0).unwrap(),
            min_sample: None,
            max_sample: None,
        };
        let dataset = Dataset::dense(
            DatasetKind::TimeSeries,
            surface(2),
            vec![
                column("t0", ColumnRole::TimePoint, vec![1.0, 10.0])
                    .with_units(Some("signal".into()))
                    .with_recorded_range(Some(recorded)),
                DataColumn::new("quality", ColumnRole::Mask, ColumnData::Int32(vec![7, 8]))
                    .unwrap(),
                column("t1", ColumnRole::TimePoint, vec![2.0, 20.0])
                    .with_units(Some("signal".into()))
                    .with_recorded_range(Some(recorded)),
            ],
        )
        .unwrap()
        .with_time_step_seconds(Some(2.0))
        .unwrap();
        let mask = SampleMask::new(dataset.domain().clone(), vec![true, false]).unwrap();

        let transformed = transform_time_series(
            &dataset,
            Some(&mask),
            UnselectedSeries::Preserve,
            ValueMetadataPolicy::DiscardValueMetadata,
            |context, series| {
                for value in series {
                    *value += context.sample as f64 * 100.0 + 1.0;
                }
                Ok(())
            },
        )
        .unwrap();

        assert_eq!(
            transformed.columns()[0].values(),
            &ColumnData::Float64(vec![2.0, 10.0])
        );
        assert_eq!(
            transformed.columns()[2].values(),
            &ColumnData::Float64(vec![3.0, 20.0])
        );
        assert_eq!(
            transformed.columns()[1].values(),
            &ColumnData::Int32(vec![7, 8])
        );
        assert_eq!(transformed.columns()[0].units(), None);
        assert_eq!(transformed.columns()[0].range_report().recorded, None);
        assert_eq!(transformed.time_step_seconds(), Some(2.0));
        assert_eq!(dataset.columns()[0].values().get_f64(0), Some(1.0));

        let filled = transform_time_series(
            &dataset,
            Some(&mask),
            UnselectedSeries::Fill(f64::NAN),
            ValueMetadataPolicy::Preserve,
            |_, _| Ok(()),
        )
        .unwrap();
        assert!(filled.columns()[0].values().get_f64(1).unwrap().is_nan());
        assert!(filled.columns()[2].values().get_f64(1).unwrap().is_nan());
        assert_eq!(filled.columns()[0].units(), Some("signal"));
    }

    #[test]
    fn mask_domain_and_callback_errors_return_without_a_partial_dataset() {
        let dataset = Dataset::dense(
            DatasetKind::TimeSeries,
            surface(2),
            vec![
                column("t0", ColumnRole::TimePoint, vec![1.0, 2.0]),
                column("t1", ColumnRole::TimePoint, vec![3.0, 4.0]),
            ],
        )
        .unwrap();
        let volume_mask = SampleMask::new(
            Domain::Volume(VolumeDomain::new(None, [2, 1, 1], None).unwrap()),
            vec![true, true],
        )
        .unwrap();
        assert!(matches!(
            summarize_time_series_with(
                &dataset,
                Some(&volume_mask),
                [SummaryOutput::new("bad", ColumnRole::Generic, 0.0).unwrap()],
                |_, _| Ok([0.0])
            ),
            Err(Error::DomainMismatch { .. })
        ));

        let result = transform_time_series(
            &dataset,
            None,
            UnselectedSeries::Preserve,
            ValueMetadataPolicy::DiscardValueMetadata,
            |context, series| {
                series[0] = 999.0;
                if context.row == 1 {
                    Err(Error::InvalidParameter {
                        name: "callback".into(),
                        reason: "test failure".into(),
                    })
                } else {
                    Ok(())
                }
            },
        );
        assert!(matches!(
            result,
            Err(Error::InvalidParameter { ref name, .. }) if name == "callback"
        ));
        assert_eq!(dataset.columns()[0].values().get_f64(0), Some(1.0));
        assert_eq!(dataset.columns()[0].values().get_f64(1), Some(2.0));
    }
}
