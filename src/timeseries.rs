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
// A read-only, file-neutral view that turns Dataset's column-major storage
// (one column per time point) into the row/time-series access used by AFNI
// programs (one series per voxel or surface node).
//
// WHY IT IS A VIEW
//
// Dataset deliberately keeps each typed column intact. Physically transposing
// all values merely to inspect a series would duplicate a potentially large
// dataset. TimeSeriesView instead borrows the Dataset. Its `series` method
// returns an iterator that visits one value from each selected column, while
// `copy_series_into` fills a caller-owned contiguous work buffer for algorithms
// such as detrending and FFTs.
//
// COLUMN SELECTION
//
// A TimeSeries dataset normally marks each temporal column as TimePoint. Some
// real file formats do not carry that role and arrive as Generic columns. To
// match the established InstaCorr behavior, explicit TimePoint columns win; if
// none exist, numeric Generic columns are used. Text Generic columns and other
// roles are not silently treated as time points.
// ---------------------------------------------------------------------------

//! Row-oriented access to a [`Dataset`] time series.

use std::iter::FusedIterator;
use std::slice;

use crate::column::{ColumnRole, DataColumn};
use crate::dataset::{Dataset, DatasetKind};
use crate::error::{Error, Result};

/// Which dataset columns form the temporal axis.
///
/// This is private because callers should not need to reproduce the fallback
/// rule. Once a view exists, all of its methods use the same selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TimeColumnSelection {
    TimePoint,
    NumericGeneric,
}

impl TimeColumnSelection {
    /// Whether a column belongs to this view's temporal axis.
    fn includes(self, column: &DataColumn) -> bool {
        match self {
            Self::TimePoint => column.role() == &ColumnRole::TimePoint,
            Self::NumericGeneric => {
                column.role() == &ColumnRole::Generic && column.values().is_numeric()
            }
        }
    }
}

/// A borrowed row-oriented view of a time-series [`Dataset`].
///
/// A dataset stores one column per time point. This view presents one iterator
/// per spatial row without copying or transposing the dataset. Numeric storage
/// types are widened to `f64`, the common type used by `afni-core` signal
/// algorithms. An `i64` value outside the exactly representable `f64` range may
/// therefore lose integer precision; neuroimaging time-series amplitudes are
/// expected to be floating-point values.
///
/// Construct a view with [`Dataset::time_series`](Dataset::time_series) or
/// [`TimeSeriesView::new`].
#[derive(Debug, Clone, Copy)]
pub struct TimeSeriesView<'a> {
    dataset: &'a Dataset,
    selection: TimeColumnSelection,
    time_point_count: usize,
}

impl<'a> TimeSeriesView<'a> {
    /// Validate and borrow a time-series dataset.
    ///
    /// The dataset must have [`DatasetKind::TimeSeries`]. Explicit
    /// [`ColumnRole::TimePoint`] columns are used when present and must all be
    /// numeric. When there are no explicit time-point columns, numeric
    /// [`ColumnRole::Generic`] columns are used for file formats that cannot
    /// record the role. At least one usable time point is required.
    pub fn new(dataset: &'a Dataset) -> Result<Self> {
        if dataset.kind() != &DatasetKind::TimeSeries {
            return Err(Error::InvalidParameter {
                name: "dataset kind".into(),
                reason: "time-series access needs a TimeSeries dataset".into(),
            });
        }

        let explicit_count = dataset
            .columns()
            .iter()
            .filter(|column| column.role() == &ColumnRole::TimePoint)
            .count();
        let (selection, time_point_count) = if explicit_count == 0 {
            let count = dataset
                .columns()
                .iter()
                .filter(|column| {
                    column.role() == &ColumnRole::Generic && column.values().is_numeric()
                })
                .count();
            (TimeColumnSelection::NumericGeneric, count)
        } else {
            if let Some(column) = dataset.columns().iter().find(|column| {
                column.role() == &ColumnRole::TimePoint && !column.values().is_numeric()
            }) {
                return Err(Error::InvalidParameter {
                    name: "time-point column".into(),
                    reason: format!("column '{}' is text, not numeric", column.label()),
                });
            }
            (TimeColumnSelection::TimePoint, explicit_count)
        };

        if time_point_count == 0 {
            return Err(Error::InvalidParameter {
                name: "time points".into(),
                reason: "the dataset has no numeric TimePoint or Generic columns".into(),
            });
        }

        Ok(Self {
            dataset,
            selection,
            time_point_count,
        })
    }

    /// The dataset being viewed.
    pub fn dataset(&self) -> &'a Dataset {
        self.dataset
    }

    /// Number of spatial series, equal to the dataset's row count.
    ///
    /// For a sparse dataset this is the number of stored rows, not the full
    /// number of samples in its domain. Use [`sample_for_series`](Self::sample_for_series)
    /// to recover the voxel or surface-node index represented by a row.
    pub fn series_count(&self) -> usize {
        self.dataset.row_count()
    }

    /// Number of temporal values in every series.
    pub fn time_point_count(&self) -> usize {
        self.time_point_count
    }

    /// Dataset column indices that form this view's temporal axis, in order.
    ///
    /// This exposes the already-validated column selection for code that must
    /// construct replacement columns. Most read-only callers only need
    /// [`series`](Self::series) or [`copy_series_into`](Self::copy_series_into).
    pub fn column_indices(&self) -> impl Iterator<Item = usize> + '_ {
        self.dataset
            .columns()
            .iter()
            .enumerate()
            .filter_map(move |(index, column)| self.selection.includes(column).then_some(index))
    }

    /// Seconds between time points, when the dataset supplies it.
    pub fn time_step_seconds(&self) -> Option<f64> {
        self.dataset.time_step_seconds()
    }

    /// Time of the first point in seconds, when the dataset supplies it.
    pub fn time_start_seconds(&self) -> Option<f64> {
        self.dataset.time_start_seconds()
    }

    /// Domain sample represented by series row `row`.
    ///
    /// Dense data return the same number as `row`; sparse data return the
    /// corresponding entry from their [`SampleMap`](crate::mapping::SampleMap).
    pub fn sample_for_series(&self, row: usize) -> Result<u32> {
        self.check_row(row)?;
        self.dataset
            .sample_for_row(row)
            .ok_or(Error::IndexOverflow { value: row })
    }

    /// Iterate over one row's temporal values without allocating.
    ///
    /// Values follow dataset column order and are widened to `f64`. The
    /// iterator is exact-sized, so callers can inspect its remaining length or
    /// collect it efficiently when they do need an owned vector.
    pub fn series(&self, row: usize) -> Result<TimeSeriesIter<'a>> {
        self.check_row(row)?;
        Ok(self.series_unchecked(row))
    }

    /// Find and iterate the series for one domain sample.
    ///
    /// `Ok(None)` means a valid domain sample has no row in a sparse dataset.
    /// An out-of-domain sample is an error, keeping "missing" distinct from
    /// "invalid index".
    pub fn series_for_sample(&self, sample: u32) -> Result<Option<TimeSeriesIter<'a>>> {
        if sample as usize >= self.dataset.domain().sample_count() {
            return Err(Error::IndexOutOfRange {
                index: i64::from(sample),
                len: self.dataset.domain().sample_count(),
            });
        }
        Ok(self
            .dataset
            .map()
            .row_for_sample(sample)
            .map(|row| self.series_unchecked(row)))
    }

    /// Copy one series into a reusable contiguous work buffer.
    ///
    /// The destination length must equal [`time_point_count`](Self::time_point_count).
    /// Reusing one buffer is the preferred path for algorithms that mutate a
    /// series in place, such as detrending, filtering, and Fourier transforms.
    pub fn copy_series_into(&self, row: usize, output: &mut [f64]) -> Result<()> {
        self.check_row(row)?;
        if output.len() != self.time_point_count {
            return Err(Error::LengthMismatch {
                what: "time-series output buffer".into(),
                expected: self.time_point_count,
                found: output.len(),
            });
        }
        for (destination, value) in output.iter_mut().zip(self.series_unchecked(row)) {
            *destination = value;
        }
        Ok(())
    }

    /// Copy one series into a newly allocated vector.
    ///
    /// Prefer [`series`](Self::series) for read-only calculations and
    /// [`copy_series_into`](Self::copy_series_into) inside a loop that can reuse
    /// its buffer.
    pub fn copy_series(&self, row: usize) -> Result<Vec<f64>> {
        Ok(self.series(row)?.collect())
    }

    /// Report an out-of-range row consistently before constructing an iterator.
    fn check_row(&self, row: usize) -> Result<()> {
        if row >= self.series_count() {
            return Err(Error::IndexOutOfRange {
                index: i64::try_from(row).unwrap_or(i64::MAX),
                len: self.series_count(),
            });
        }
        Ok(())
    }

    /// Construct after the caller has established that `row` is valid.
    fn series_unchecked(&self, row: usize) -> TimeSeriesIter<'a> {
        let columns: &'a [DataColumn] = self.dataset.columns();
        TimeSeriesIter {
            columns: columns.iter(),
            selection: self.selection,
            row,
            remaining: self.time_point_count,
        }
    }
}

/// Allocation-free iterator over the temporal values of one dataset row.
///
/// This type is returned by [`TimeSeriesView::series`] and
/// [`TimeSeriesView::series_for_sample`]. Callers normally use it through the
/// standard [`Iterator`] methods rather than constructing it directly.
#[derive(Debug, Clone)]
pub struct TimeSeriesIter<'a> {
    columns: slice::Iter<'a, DataColumn>,
    selection: TimeColumnSelection,
    row: usize,
    remaining: usize,
}

impl Iterator for TimeSeriesIter<'_> {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        for column in self.columns.by_ref() {
            if !self.selection.includes(column) {
                continue;
            }
            // TimeSeriesView::new validated the selected columns as numeric,
            // and Dataset::new validated that every column contains this row.
            let value = column
                .values()
                .get_f64(self.row)
                .expect("validated numeric time-series column and row");
            self.remaining -= 1;
            return Some(value);
        }
        debug_assert_eq!(self.remaining, 0);
        None
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl ExactSizeIterator for TimeSeriesIter<'_> {}
impl FusedIterator for TimeSeriesIter<'_> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::column::ColumnData;
    use crate::domain::{Domain, SurfaceDomain};

    fn surface(samples: usize) -> Domain {
        Domain::Surface(SurfaceDomain::new(None, samples).unwrap())
    }

    fn column(label: &str, role: ColumnRole, values: ColumnData) -> DataColumn {
        DataColumn::new(label, role, values).unwrap()
    }

    #[test]
    fn explicit_time_points_form_all_dense_series_without_transposing() {
        let dataset = Dataset::dense(
            DatasetKind::TimeSeries,
            surface(2),
            vec![
                column(
                    "t0",
                    ColumnRole::TimePoint,
                    ColumnData::Float32(vec![1.5, 10.5]),
                ),
                column(
                    "ignored",
                    ColumnRole::Generic,
                    ColumnData::Float64(vec![99.0, 99.0]),
                ),
                column("t1", ColumnRole::TimePoint, ColumnData::Int32(vec![2, 20])),
            ],
        )
        .unwrap()
        .with_time_step_seconds(Some(2.0))
        .unwrap()
        .with_time_start_seconds(Some(0.5))
        .unwrap();

        let view = TimeSeriesView::new(&dataset).unwrap();
        assert_eq!((view.series_count(), view.time_point_count()), (2, 2));
        assert_eq!(
            (view.time_start_seconds(), view.time_step_seconds()),
            (Some(0.5), Some(2.0))
        );
        assert_eq!(view.series(0).unwrap().collect::<Vec<_>>(), [1.5, 2.0]);
        assert_eq!(view.copy_series(1).unwrap(), [10.5, 20.0]);
        assert_eq!(view.sample_for_series(1).unwrap(), 1);
    }

    #[test]
    fn numeric_generic_columns_are_the_file_format_fallback() {
        let dataset = Dataset::dense(
            DatasetKind::TimeSeries,
            surface(2),
            vec![
                column("g0", ColumnRole::Generic, ColumnData::UInt32(vec![1, 2])),
                column(
                    "text",
                    ColumnRole::Generic,
                    ColumnData::Text(vec!["a".into(), "b".into()]),
                ),
                column(
                    "mask",
                    ColumnRole::Mask,
                    ColumnData::Float32(vec![1.0, 1.0]),
                ),
                column(
                    "g1",
                    ColumnRole::Generic,
                    ColumnData::Float64(vec![3.0, 4.0]),
                ),
            ],
        )
        .unwrap();

        let view = dataset.time_series().unwrap();
        assert_eq!(view.time_point_count(), 2);
        assert_eq!(view.copy_series(0).unwrap(), [1.0, 3.0]);
        assert_eq!(view.copy_series(1).unwrap(), [2.0, 4.0]);
    }

    #[test]
    fn copy_into_reuses_an_exact_sized_buffer_and_checks_rows() {
        let dataset = Dataset::dense(
            DatasetKind::TimeSeries,
            surface(1),
            vec![
                column("t0", ColumnRole::TimePoint, ColumnData::Float64(vec![2.0])),
                column("t1", ColumnRole::TimePoint, ColumnData::Float64(vec![4.0])),
            ],
        )
        .unwrap();
        let view = dataset.time_series().unwrap();
        let mut buffer = vec![0.0; view.time_point_count()];
        view.copy_series_into(0, &mut buffer).unwrap();
        assert_eq!(buffer, [2.0, 4.0]);

        assert!(matches!(
            view.copy_series_into(0, &mut [0.0]),
            Err(Error::LengthMismatch {
                expected: 2,
                found: 1,
                ..
            })
        ));
        assert!(matches!(
            view.series(1),
            Err(Error::IndexOutOfRange { index: 1, len: 1 })
        ));
    }

    #[test]
    fn sparse_rows_retain_their_domain_sample_identity() {
        let dataset = Dataset::indexed(
            DatasetKind::TimeSeries,
            surface(8),
            vec![6, 2],
            vec![
                column(
                    "t0",
                    ColumnRole::TimePoint,
                    ColumnData::Float64(vec![60.0, 20.0]),
                ),
                column(
                    "t1",
                    ColumnRole::TimePoint,
                    ColumnData::Float64(vec![61.0, 21.0]),
                ),
            ],
        )
        .unwrap();
        let view = dataset.time_series().unwrap();

        assert_eq!(view.sample_for_series(0).unwrap(), 6);
        assert_eq!(
            view.series_for_sample(2)
                .unwrap()
                .unwrap()
                .collect::<Vec<_>>(),
            [20.0, 21.0]
        );
        assert!(view.series_for_sample(3).unwrap().is_none());
        assert!(matches!(
            view.series_for_sample(8),
            Err(Error::IndexOutOfRange { index: 8, len: 8 })
        ));
    }

    #[test]
    fn construction_rejects_wrong_kind_or_unusable_time_columns() {
        let scalar = Dataset::dense(
            DatasetKind::Scalar,
            surface(1),
            vec![column(
                "x",
                ColumnRole::Generic,
                ColumnData::Float64(vec![1.0]),
            )],
        )
        .unwrap();
        assert!(matches!(
            scalar.time_series(),
            Err(Error::InvalidParameter { ref name, .. }) if name == "dataset kind"
        ));

        let text_time = Dataset::dense(
            DatasetKind::TimeSeries,
            surface(1),
            vec![column(
                "words",
                ColumnRole::TimePoint,
                ColumnData::Text(vec!["not a number".into()]),
            )],
        )
        .unwrap();
        assert!(matches!(
            text_time.time_series(),
            Err(Error::InvalidParameter { ref name, .. }) if name == "time-point column"
        ));

        let no_time = Dataset::dense(
            DatasetKind::TimeSeries,
            surface(1),
            vec![column(
                "labels",
                ColumnRole::Label,
                ColumnData::Int32(vec![1]),
            )],
        )
        .unwrap();
        assert!(matches!(
            no_time.time_series(),
            Err(Error::InvalidParameter { ref name, .. }) if name == "time points"
        ));
    }
}
