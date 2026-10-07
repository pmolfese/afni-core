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
// One column of a dataset: typed values (`ColumnData`) plus typed metadata
// (`DataColumn`): a label, a role, units, a statistic (`StatSpec`), FDR/MDF
// curves, a label table, and the value range recorded in the source file.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `dataset.rs` holds a list of `DataColumn`s and checks they all have the
//   right number of rows.
// * `stat.rs`, `curve.rs`, `labels.rs` supply the metadata types stored here.
// * `numeric.rs` supplies the finiteness rules used for ranges.
//
// POLICIES
//
// * Values keep their own type. Label keys stay integers; doubles stay `f64`.
//   Nothing is coerced to `f32`.
// * NaN / infinity in float columns are kept as stored. They mean "no data".
//   `computed_range` and `finite_count` ignore them; nothing replaces them.
// * A column must have at least one row and a non-blank label.
// * The range RECORDED by the source file is stored separately from the range
//   COMPUTED from the data. `range_report` returns both so a disagreement
//   (stale header, edited data) is visible instead of silently "fixed".
// ---------------------------------------------------------------------------

//! Typed column storage and column metadata.

use crate::curve::ThresholdCurve;
use crate::error::{Error, Result};
use crate::labels::LabelTable;
use crate::mapping::SampleMap;
use crate::numeric::ensure_finite;
use crate::stat::StatSpec;

/// A closed interval `[min, max]` of data values.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ColumnRange {
    /// Smallest value.
    pub min: f64,
    /// Largest value.
    pub max: f64,
}

impl ColumnRange {
    /// A validated range: both ends finite and `min <= max`.
    pub fn new(min: f64, max: f64) -> Result<Self> {
        ensure_finite("range min", min)?;
        ensure_finite("range max", max)?;
        if min > max {
            return Err(Error::InvalidParameter {
                name: "range".into(),
                reason: format!("min {min} is greater than max {max}"),
            });
        }
        Ok(Self { min, max })
    }

    /// Whether `value` lies inside the closed interval (ends included).
    pub fn contains(&self, value: f64) -> bool {
        value >= self.min && value <= self.max
    }

    /// Position of `value` within the range as 0..=1, clamped. A zero-width
    /// range maps everything to 0.5 so callers never divide by zero.
    pub fn normalized(&self, value: f64) -> f64 {
        if (self.max - self.min).abs() <= f64::EPSILON {
            0.5
        } else {
            ((value - self.min) / (self.max - self.min)).clamp(0.0, 1.0)
        }
    }
}

/// A range as recorded in a file (AFNI `COLMS_RANGE`): the extremes plus, when
/// the file says so, the samples where they occur.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RecordedRange {
    /// The recorded extremes.
    pub range: ColumnRange,
    /// Domain sample holding the minimum, if recorded and valid.
    pub min_sample: Option<u32>,
    /// Domain sample holding the maximum, if recorded and valid.
    pub max_sample: Option<u32>,
}

/// The recorded range next to the freshly computed one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RangeReport {
    /// What the source file claimed, if anything.
    pub recorded: Option<RecordedRange>,
    /// Computed from the finite values now, if there are any (or the column is
    /// numeric).
    pub computed: Option<ColumnRange>,
}

impl RangeReport {
    /// True only when both exist and differ by more than `rel_tol` (relative
    /// to the larger magnitude, with an absolute floor of `rel_tol`).
    ///
    /// AFNI writes ranges as 32-bit floats, so use roughly `1e-6` when the
    /// data are `f64`.
    pub fn disagrees(&self, rel_tol: f64) -> bool {
        let (Some(rec), Some(comp)) = (self.recorded, self.computed) else {
            return false;
        };
        let differ = |a: f64, b: f64| (a - b).abs() > rel_tol * a.abs().max(b.abs()).max(1.0);
        differ(rec.range.min, comp.min) || differ(rec.range.max, comp.max)
    }
}

/// How a column is meant to be used.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ColumnRole {
    /// Node indices stored as data.
    NodeIndex,
    /// Values to color (an overlay intensity).
    Intensity,
    /// Values to threshold on.
    Threshold,
    /// Values that modulate brightness.
    Brightness,
    /// Integer label keys.
    Label,
    /// A statistic (see the column's [`StatSpec`]).
    Statistic,
    /// One time point of a time series.
    TimePoint,
    /// A 0/non-zero mask.
    Mask,
    /// Ordinary numbers with no stated purpose (AFNI `Generic_*`).
    Generic,
    /// A role this crate has no name for; the source's own text is kept.
    Other(String),
}

/// Typed values of one column.
#[derive(Debug, Clone, PartialEq)]
pub enum ColumnData {
    /// 32-bit signed integers.
    Int32(Vec<i32>),
    /// 32-bit unsigned integers.
    UInt32(Vec<u32>),
    /// 64-bit signed integers.
    Int64(Vec<i64>),
    /// 32-bit floats.
    Float32(Vec<f32>),
    /// 64-bit floats.
    Float64(Vec<f64>),
    /// Text.
    Text(Vec<String>),
}

/// What to do with metadata when replacing or transforming column values.
///
/// A recorded range, statistic, FDR/MDF curve, label table, or unit can become
/// false after arithmetic changes the values. Requiring this choice prevents a
/// generic transformation from silently carrying stale scientific meaning.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValueMetadataPolicy {
    /// Keep all metadata.
    ///
    /// Use this for representation-only changes such as sparse-to-dense
    /// expansion, or when the caller knows the replacement has identical
    /// scientific meaning.
    Preserve,
    /// Clear units, statistic metadata, FDR/MDF curves, label tables, and the
    /// recorded range.
    ///
    /// The column label and role remain because they identify the column's
    /// place in the dataset. Call the existing `with_*` builders to attach new
    /// metadata after the transformation.
    DiscardValueMetadata,
}

/// Values used for samples that have no row when expanding sparse data.
#[derive(Debug, Clone, PartialEq)]
pub struct MissingFill {
    /// For float columns. Default NaN ("no data").
    pub float: f64,
    /// For integer columns. Default 0; must fit the column's integer type.
    pub int: i64,
    /// For text columns. Default empty.
    pub text: String,
}

impl Default for MissingFill {
    fn default() -> Self {
        Self {
            float: f64::NAN,
            int: 0,
            text: String::new(),
        }
    }
}

impl ColumnData {
    /// Number of values.
    pub fn len(&self) -> usize {
        match self {
            Self::Int32(v) => v.len(),
            Self::UInt32(v) => v.len(),
            Self::Int64(v) => v.len(),
            Self::Float32(v) => v.len(),
            Self::Float64(v) => v.len(),
            Self::Text(v) => v.len(),
        }
    }

    /// Whether there are no values.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Value `row` as `f64`, or `None` if out of range or the column is text.
    /// Integers beyond 2^53 lose precision here; use the typed variants when
    /// exactness matters.
    pub fn get_f64(&self, row: usize) -> Option<f64> {
        match self {
            Self::Int32(v) => v.get(row).map(|&x| f64::from(x)),
            Self::UInt32(v) => v.get(row).map(|&x| f64::from(x)),
            Self::Int64(v) => v.get(row).map(|&x| x as f64),
            Self::Float32(v) => v.get(row).map(|&x| f64::from(x)),
            Self::Float64(v) => v.get(row).copied(),
            Self::Text(_) => None,
        }
    }

    /// Whether the column holds numbers (not text).
    pub fn is_numeric(&self) -> bool {
        !matches!(self, Self::Text(_))
    }

    /// How many values are finite. Integers always are; floats exclude
    /// NaN/infinity; text columns have none.
    pub fn finite_count(&self) -> usize {
        match self {
            Self::Text(_) => 0,
            Self::Int32(v) => v.len(),
            Self::UInt32(v) => v.len(),
            Self::Int64(v) => v.len(),
            Self::Float32(v) => v.iter().filter(|x| x.is_finite()).count(),
            Self::Float64(v) => v.iter().filter(|x| x.is_finite()).count(),
        }
    }

    /// Min and max over the finite values, or `None` if there are none (empty,
    /// all non-finite, or text).
    pub fn computed_range(&self) -> Option<ColumnRange> {
        let mut range: Option<ColumnRange> = None;
        for row in 0..self.len() {
            let Some(x) = self.get_f64(row).filter(|x| x.is_finite()) else {
                continue;
            };
            range = Some(match range {
                None => ColumnRange { min: x, max: x },
                Some(r) => ColumnRange {
                    min: r.min.min(x),
                    max: r.max.max(x),
                },
            });
        }
        range
    }

    /// Equality that treats NaN as equal to NaN.
    ///
    /// The derived `==` follows IEEE (`NaN != NaN`), so a column with missing
    /// samples is never `==` to its own clone. Use this when "same data,
    /// including the same missing samples" is what you mean. Floats are
    /// compared by value (so `0.0 == -0.0`), NaN payloads are ignored, and
    /// differing variants (e.g. `Float32` vs `Float64`) are never equal.
    pub fn eq_nan_aware(&self, other: &Self) -> bool {
        fn same<T: Copy + PartialEq>(a: &[T], b: &[T], is_nan: impl Fn(T) -> bool) -> bool {
            a.len() == b.len()
                && a.iter()
                    .zip(b)
                    .all(|(&x, &y)| x == y || (is_nan(x) && is_nan(y)))
        }
        match (self, other) {
            (Self::Float32(a), Self::Float32(b)) => same(a, b, f32::is_nan),
            (Self::Float64(a), Self::Float64(b)) => same(a, b, f64::is_nan),
            // Integers and text have no NaN, so ordinary equality is exact.
            _ => self == other,
        }
    }

    /// Rearrange rows into one value per domain sample, filling samples with no
    /// row from `fill`. See [`SampleMap::expand`].
    pub fn expand(&self, map: &SampleMap, domain_len: usize, fill: &MissingFill) -> Result<Self> {
        // Narrow the integer fill to the column's own type, with a clear error
        // if it does not fit (rather than truncating).
        fn narrow<T: TryFrom<i64>>(fill: i64) -> Result<T> {
            T::try_from(fill).map_err(|_| Error::InvalidParameter {
                name: "fill.int".into(),
                reason: format!("{fill} does not fit the column's integer type"),
            })
        }
        Ok(match self {
            Self::Int32(v) => Self::Int32(map.expand(v, domain_len, narrow(fill.int)?)?),
            Self::UInt32(v) => Self::UInt32(map.expand(v, domain_len, narrow(fill.int)?)?),
            Self::Int64(v) => Self::Int64(map.expand(v, domain_len, fill.int)?),
            Self::Float32(v) => Self::Float32(map.expand(
                v,
                domain_len,
                crate::numeric::narrow_to_f32(fill.float),
            )?),
            Self::Float64(v) => Self::Float64(map.expand(v, domain_len, fill.float)?),
            Self::Text(v) => Self::Text(map.expand(v, domain_len, fill.text.clone())?),
        })
    }
}

/// One dataset column: values plus metadata.
///
/// Fields are private so the invariants (non-blank label, at least one row)
/// hold for the column's whole life; read them through the accessors and set
/// optional metadata with the `with_*` builder methods.
#[derive(Debug, Clone, PartialEq)]
pub struct DataColumn {
    label: String,
    role: ColumnRole,
    units: Option<String>,
    stat: Option<StatSpec>,
    fdr_curve: Option<ThresholdCurve>,
    mdf_curve: Option<ThresholdCurve>,
    label_table: Option<LabelTable>,
    recorded_range: Option<RecordedRange>,
    values: ColumnData,
}

impl DataColumn {
    /// Create a column. Fails if `label` is blank or `values` is empty.
    pub fn new(label: impl Into<String>, role: ColumnRole, values: ColumnData) -> Result<Self> {
        let label = label.into();
        if label.trim().is_empty() {
            return Err(Error::Empty("column label".into()));
        }
        if values.is_empty() {
            return Err(Error::Empty(format!("column '{label}'")));
        }
        Ok(Self {
            label,
            role,
            units: None,
            stat: None,
            fdr_curve: None,
            mdf_curve: None,
            label_table: None,
            recorded_range: None,
            values,
        })
    }

    /// Set the units text (`None` clears it).
    pub fn with_units(mut self, units: Option<String>) -> Self {
        self.units = units.filter(|u| !u.trim().is_empty());
        self
    }

    /// Attach a statistic description.
    pub fn with_stat(mut self, stat: Option<StatSpec>) -> Self {
        self.stat = stat;
        self
    }

    /// Attach an FDR curve.
    pub fn with_fdr_curve(mut self, curve: Option<ThresholdCurve>) -> Self {
        self.fdr_curve = curve;
        self
    }

    /// Attach a missed-detection (MDF) curve.
    pub fn with_mdf_curve(mut self, curve: Option<ThresholdCurve>) -> Self {
        self.mdf_curve = curve;
        self
    }

    /// Attach a label table.
    pub fn with_label_table(mut self, table: Option<LabelTable>) -> Self {
        self.label_table = table;
        self
    }

    /// Record the range the source file claimed.
    pub fn with_recorded_range(mut self, range: Option<RecordedRange>) -> Self {
        self.recorded_range = range;
        self
    }

    /// The column label.
    pub fn label(&self) -> &str {
        &self.label
    }
    /// The column's role.
    pub fn role(&self) -> &ColumnRole {
        &self.role
    }
    /// Units text, if any.
    pub fn units(&self) -> Option<&str> {
        self.units.as_deref()
    }
    /// The statistic this column holds, if any.
    pub fn stat(&self) -> Option<&StatSpec> {
        self.stat.as_ref()
    }
    /// The FDR curve, if any.
    pub fn fdr_curve(&self) -> Option<&ThresholdCurve> {
        self.fdr_curve.as_ref()
    }
    /// The MDF curve, if any.
    pub fn mdf_curve(&self) -> Option<&ThresholdCurve> {
        self.mdf_curve.as_ref()
    }
    /// The label table, if any.
    pub fn label_table(&self) -> Option<&LabelTable> {
        self.label_table.as_ref()
    }
    /// The typed values.
    pub fn values(&self) -> &ColumnData {
        &self.values
    }

    /// Return a column with replacement values and an explicit metadata policy.
    ///
    /// This validates only the column-level invariant that values are nonempty.
    /// A standalone column does not know its dataset row count; use
    /// [`Dataset::replace_column`](crate::dataset::Dataset::replace_column) or
    /// [`Dataset::transform_column`](crate::dataset::Dataset::transform_column)
    /// to validate the replacement against a dataset atomically.
    pub fn with_values(&self, values: ColumnData, metadata: ValueMetadataPolicy) -> Result<Self> {
        if values.is_empty() {
            return Err(Error::Empty(format!("column '{}'", self.label)));
        }
        let mut output = Self {
            values,
            ..self.clone()
        };
        if metadata == ValueMetadataPolicy::DiscardValueMetadata {
            output.units = None;
            output.stat = None;
            output.fdr_curve = None;
            output.mdf_curve = None;
            output.label_table = None;
            output.recorded_range = None;
        }
        Ok(output)
    }

    /// Map every numeric value to an `f64` while choosing how metadata changes.
    ///
    /// Integer and `f32` input is widened through [`ColumnData::get_f64`]. An
    /// `i64` outside the exactly representable `f64` range may lose precision;
    /// callers that require exact integer arithmetic should match the typed
    /// [`ColumnData`] variant and pass the result to [`with_values`](Self::with_values).
    /// Text columns are rejected. Non-finite input or output is retained, just
    /// as it is in ordinary floating-point columns.
    pub fn map_numeric_to_f64(
        &self,
        metadata: ValueMetadataPolicy,
        mut transform: impl FnMut(f64) -> f64,
    ) -> Result<Self> {
        self.try_map_numeric_to_f64(metadata, |_, value| Ok(transform(value)))
    }

    /// Fallible, row-aware form of [`map_numeric_to_f64`](Self::map_numeric_to_f64).
    ///
    /// Evaluation stops at the first error and no column is returned. The
    /// zero-based row index lets transformations report or treat particular
    /// rows differently without maintaining an external counter.
    pub fn try_map_numeric_to_f64(
        &self,
        metadata: ValueMetadataPolicy,
        mut transform: impl FnMut(usize, f64) -> Result<f64>,
    ) -> Result<Self> {
        if !self.values.is_numeric() {
            return Err(Error::InvalidParameter {
                name: "column values".into(),
                reason: format!("column '{}' is text, not numeric", self.label),
            });
        }
        let values = (0..self.len())
            .map(|row| {
                let value = self
                    .values
                    .get_f64(row)
                    .expect("validated numeric column and in-range row");
                transform(row, value)
            })
            .collect::<Result<Vec<_>>>()?;
        self.with_values(ColumnData::Float64(values), metadata)
    }
    /// Number of rows.
    pub fn len(&self) -> usize {
        self.values.len()
    }
    /// Always false (a column has at least one row); provided to pair with `len`.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Equality of the whole column (metadata and values) with NaN equal to
    /// NaN; see [`ColumnData::eq_nan_aware`].
    pub fn eq_nan_aware(&self, other: &Self) -> bool {
        // Compare each metadata field with ordinary `==` (none of them hold
        // data samples), then the values NaN-aware.
        self.label == other.label
            && self.role == other.role
            && self.units == other.units
            && self.stat == other.stat
            && self.fdr_curve == other.fdr_curve
            && self.mdf_curve == other.mdf_curve
            && self.label_table == other.label_table
            && self.recorded_range == other.recorded_range
            && self.values.eq_nan_aware(&other.values)
    }

    /// The FDR q-value at `threshold` (an absolute statistic), from this column's
    /// stored FDR curve; `None` if the column has no curve. See
    /// [`crate::fdr::q_value_for_threshold`].
    pub fn q_for_threshold(&self, threshold: f64) -> Option<Result<crate::fdr::QValue>> {
        self.fdr_curve
            .as_ref()
            .map(|c| crate::fdr::q_value_for_threshold(c, threshold))
    }

    /// The statistic threshold that achieves FDR q-value `q`, from this column's
    /// stored FDR curve; `None` if the column has no curve. The largest absolute
    /// value in the column is used as AFNI uses it (for a `q` smaller than the
    /// curve can express). See [`crate::fdr::threshold_for_q`].
    pub fn threshold_for_q(&self, q: crate::fdr::QValue) -> Option<Result<f64>> {
        let curve = self.fdr_curve.as_ref()?;
        let max_abs = (0..self.values.len())
            .filter_map(|i| self.values.get_f64(i))
            .filter(|v| v.is_finite())
            .map(f64::abs)
            .fold(None, |m: Option<f64>, v| Some(m.map_or(v, |m| m.max(v))));
        Some(crate::fdr::threshold_for_q(curve, q, max_abs))
    }

    /// The missed-detection fraction at p-value `p`, from this column's MDF
    /// curve; `None` if the column has none. See
    /// [`crate::fdr::missed_detection_fraction`].
    pub fn missed_detection_fraction(&self, p: f64) -> Option<Result<f64>> {
        self.mdf_curve
            .as_ref()
            .map(|c| crate::fdr::missed_detection_fraction(c, p))
    }

    /// The recorded range next to the computed one.
    pub fn range_report(&self) -> RangeReport {
        RangeReport {
            recorded: self.recorded_range,
            computed: self.values.computed_range(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn col(values: ColumnData) -> DataColumn {
        DataColumn::new("c", ColumnRole::Generic, values).unwrap()
    }

    #[test]
    fn construction_validates_label_and_rows() {
        assert!(DataColumn::new(" ", ColumnRole::Generic, ColumnData::Int32(vec![1])).is_err());
        assert!(DataColumn::new("c", ColumnRole::Generic, ColumnData::Int32(vec![])).is_err());
    }

    #[test]
    fn computed_range_ignores_non_finite() {
        let c = col(ColumnData::Float32(vec![
            f32::NAN,
            2.0,
            -1.5,
            f32::INFINITY,
        ]));
        assert_eq!(
            c.values().computed_range(),
            Some(ColumnRange {
                min: -1.5,
                max: 2.0
            })
        );
        assert_eq!(c.values().finite_count(), 2);
        assert_eq!(
            col(ColumnData::Float64(vec![f64::NAN]))
                .values()
                .computed_range(),
            None
        );
        assert_eq!(
            col(ColumnData::Text(vec!["a".into()]))
                .values()
                .computed_range(),
            None
        );
    }

    #[test]
    fn integers_keep_exact_values() {
        let big = (1_i64 << 53) + 1;
        let c = col(ColumnData::Int64(vec![big]));
        assert_eq!(c.values(), &ColumnData::Int64(vec![big]));
    }

    #[test]
    fn range_report_preserves_both_when_they_disagree() {
        let rec = RecordedRange {
            range: ColumnRange::new(0.0, 10.0).unwrap(),
            min_sample: Some(0),
            max_sample: None,
        };
        let c = col(ColumnData::Float64(vec![1.0, 5.0])).with_recorded_range(Some(rec));
        let report = c.range_report();
        assert_eq!(report.recorded, Some(rec));
        assert_eq!(report.computed, Some(ColumnRange { min: 1.0, max: 5.0 }));
        assert!(report.disagrees(1e-6));
        // No recorded range: nothing to disagree with.
        assert!(!col(ColumnData::Float64(vec![1.0]))
            .range_report()
            .disagrees(1e-6));
    }

    #[test]
    fn nan_aware_equality() {
        let a = ColumnData::Float32(vec![1.0, f32::NAN]);
        assert_ne!(a, a.clone(), "derived == follows IEEE");
        assert!(a.eq_nan_aware(&a.clone()));
        assert!(!a.eq_nan_aware(&ColumnData::Float32(vec![1.0, 2.0])));
        assert!(!a.eq_nan_aware(&ColumnData::Float32(vec![1.0])));
        assert!(!a.eq_nan_aware(&ColumnData::Float64(vec![1.0, f64::NAN])));
        assert!(ColumnData::Int32(vec![1]).eq_nan_aware(&ColumnData::Int32(vec![1])));
        // Column level: metadata differences still count.
        let c = col(a.clone());
        assert!(c.eq_nan_aware(&c.clone()));
        assert!(!c.eq_nan_aware(&c.clone().with_units(Some("mm".into()))));
    }

    #[test]
    fn fdr_conveniences_need_a_curve() {
        let plain = col(ColumnData::Float64(vec![1.0, 2.0]));
        assert!(plain.q_for_threshold(1.0).is_none());
        assert!(plain
            .threshold_for_q(crate::fdr::QValue::new(0.05).unwrap())
            .is_none());
        assert!(plain.missed_detection_fraction(0.01).is_none());
        let curve = ThresholdCurve::new(0.0, 1.0, vec![0.0, 1.0, 2.0, 3.0]).unwrap();
        let with = col(ColumnData::Float64(vec![-9.0, 2.0])).with_fdr_curve(Some(curve.clone()));
        let q = with.q_for_threshold(2.0).unwrap().unwrap().get();
        assert!((q - crate::fdr::q_for_z(2.0).unwrap().get()).abs() < 1e-15);
        // A q smaller than the curve allows: the column's largest |value| (9) wins.
        let thr = with
            .threshold_for_q(crate::fdr::QValue::new(1e-12).unwrap())
            .unwrap()
            .unwrap();
        assert_eq!(thr, 9.0 * 1.000002);
        let mdf = col(ColumnData::Float64(vec![1.0])).with_mdf_curve(Some(
            ThresholdCurve::new(-3.0, 1.0, vec![1.0, 0.5, 0.1, 0.0]).unwrap(),
        ));
        assert_eq!(mdf.missed_detection_fraction(1.0).unwrap().unwrap(), 0.0);
    }

    #[test]
    fn range_helpers() {
        assert!(ColumnRange::new(2.0, 1.0).is_err());
        assert!(ColumnRange::new(f64::NAN, 1.0).is_err());
        let r = ColumnRange::new(0.0, 10.0).unwrap();
        assert!(r.contains(0.0) && r.contains(10.0) && !r.contains(10.1));
        assert_eq!(r.normalized(5.0), 0.5);
        assert_eq!(r.normalized(20.0), 1.0);
        assert_eq!(ColumnRange::new(3.0, 3.0).unwrap().normalized(9.0), 0.5);
    }

    #[test]
    fn expand_uses_fill_per_type() {
        let map = SampleMap::indexed(vec![2], 4).unwrap();
        let fill = MissingFill::default();
        let f = ColumnData::Float32(vec![1.5])
            .expand(&map, 4, &fill)
            .unwrap();
        assert_eq!(f.get_f64(2), Some(1.5));
        assert!(f.get_f64(0).unwrap().is_nan());
        let i = ColumnData::Int32(vec![9]).expand(&map, 4, &fill).unwrap();
        assert_eq!(i, ColumnData::Int32(vec![0, 0, 9, 0]));
        let bad = MissingFill {
            int: i64::MAX,
            ..MissingFill::default()
        };
        assert!(ColumnData::Int32(vec![9]).expand(&map, 4, &bad).is_err());
    }

    #[test]
    fn replacing_values_requires_an_explicit_metadata_policy() {
        let recorded = RecordedRange {
            range: ColumnRange::new(1.0, 2.0).unwrap(),
            min_sample: Some(0),
            max_sample: Some(1),
        };
        let curve = ThresholdCurve::new(0.0, 1.0, vec![0.0, 1.0]).unwrap();
        let source = DataColumn::new(
            "statistic",
            ColumnRole::Statistic,
            ColumnData::Float32(vec![1.0, 2.0]),
        )
        .unwrap()
        .with_units(Some("z".into()))
        .with_stat(Some(StatSpec::new(crate::stat::StatKind::Zscore, &[], 0.0)))
        .with_fdr_curve(Some(curve.clone()))
        .with_mdf_curve(Some(curve))
        .with_label_table(Some(LabelTable::default()))
        .with_recorded_range(Some(recorded));

        let preserved = source
            .with_values(
                ColumnData::Float64(vec![1.0, 2.0]),
                ValueMetadataPolicy::Preserve,
            )
            .unwrap();
        assert_eq!(preserved.label(), source.label());
        assert_eq!(preserved.role(), source.role());
        assert_eq!(preserved.units(), source.units());
        assert_eq!(preserved.stat(), source.stat());
        assert_eq!(preserved.fdr_curve(), source.fdr_curve());
        assert_eq!(preserved.mdf_curve(), source.mdf_curve());
        assert_eq!(preserved.label_table(), source.label_table());
        assert_eq!(preserved.range_report().recorded, Some(recorded));

        let discarded = source
            .with_values(
                ColumnData::Float64(vec![10.0, 20.0]),
                ValueMetadataPolicy::DiscardValueMetadata,
            )
            .unwrap();
        assert_eq!(discarded.label(), source.label());
        assert_eq!(discarded.role(), source.role());
        assert_eq!(discarded.units(), None);
        assert_eq!(discarded.stat(), None);
        assert_eq!(discarded.fdr_curve(), None);
        assert_eq!(discarded.mdf_curve(), None);
        assert_eq!(discarded.label_table(), None);
        assert_eq!(discarded.range_report().recorded, None);
        assert!(source
            .with_values(ColumnData::Float64(vec![]), ValueMetadataPolicy::Preserve)
            .is_err());
    }

    #[test]
    fn numeric_mapping_is_typed_fallible_and_row_aware() {
        let source = col(ColumnData::Int32(vec![-2, 3]));
        let doubled = source
            .map_numeric_to_f64(ValueMetadataPolicy::DiscardValueMetadata, |value| {
                value * 2.0
            })
            .unwrap();
        assert_eq!(doubled.values(), &ColumnData::Float64(vec![-4.0, 6.0]));

        let err = source.try_map_numeric_to_f64(
            ValueMetadataPolicy::DiscardValueMetadata,
            |row, value| {
                if row == 1 {
                    Err(Error::InvalidParameter {
                        name: "row".into(),
                        reason: "test failure".into(),
                    })
                } else {
                    Ok(value)
                }
            },
        );
        assert!(matches!(err, Err(Error::InvalidParameter { ref name, .. }) if name == "row"));
        assert_eq!(source.values(), &ColumnData::Int32(vec![-2, 3]));

        let text = col(ColumnData::Text(vec!["one".into()]));
        assert!(matches!(
            text.map_numeric_to_f64(ValueMetadataPolicy::Preserve, |value| value),
            Err(Error::InvalidParameter { ref name, .. }) if name == "column values"
        ));
    }
}
