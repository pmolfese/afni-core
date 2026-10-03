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
// Seed correlation ("InstaCorr") the way SUMA does it: clean every time series once
// (detrend, bandpass, remove regressors, scale to unit length), after which the
// correlation of any seed with every other series is just a dot product. This
// module owns the ORCHESTRATION; the arithmetic lives in `signal.rs`.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `signal.rs` supplies `bandpass_vectors`, `legendre_basis` and `normalize_l2`.
// * A `Dataset` of kind time series (one `TimePoint` column per sample, one row per
//   node) goes in through `prepare_dataset`; plain rows go in through
//   `prepare_rows`. A NIML file marks no column as a time point, so the numeric
//   generic columns of a time-series dataset count as the time points.
// * What comes out is a `SeedCorrelation`: the correlations, and the statistical
//   metadata that makes them a statistic: `Correl(samples, 1, removed_dof)`, where
//   `removed_dof` is how many dimensions the cleaning took out of each series. With
//   that, `stats.rs` turns a correlation into a p-value correctly; with the sample
//   count alone it would be too optimistic.
//
// THE PIPELINE (SUMA_dot.c: `SUMA_DotDetrendDset` and `SUMA_DotXform...`)
//
// * Whole dataset, with `normalize` on: a LINEAR detrend, an FFT bandpass of every
//   series, Legendre polynomial regressors `0..=polort` (plus any extra regressors)
//   filtered the same way and projected out, then every series scaled to unit L2
//   length. `removed_dof` is `THD_bandpass_vectors`' count. At least 9 samples are
//   needed even with no bandpass (AFNI's rule), and `polort + 1` must be less than
//   `samples - 3`.
// * `normalize` off: SUMA's "normalize_dset = no" switch skips everything: the raw
//   series are dotted, the result is NOT a correlation, and no statistic is attached.
// * A seed that is a series from outside the dataset is cleaned with a MEAN detrend
//   instead of a linear one (`SUMA_dot.c`, the single-vector path), then normalized.
// * A seed made of several rows (an ROI) is the mean of their PREPARED (unit-length)
//   rows, scaled to unit length again, the order AFNI's seed-blur uses. Averaging
//   unit vectors, not raw series, keeps a high-variance node from dominating.
//
// DIFFERENCES FROM AFNI (see the roadmap log)
//
// * Polynomial regressors under a linear detrend are partly annihilated by the
//   filter (the constant and linear columns shrink to rounding noise). AFNI's
//   pseudo-inverse keeps those noise directions and projects them out too, which
//   removes two extra, arbitrary dimensions that depend on 32-bit rounding. Core
//   drops them. The removed-dimension COUNT is still AFNI's, so the statistic's
//   degrees of freedom agree. The correlations differ by the effect of two arbitrary
//   directions, a few hundredths at 100 samples (measured in the conformance test).
// * A series that cleaning leaves with no variance (squared length `<= 1e-20`) is
//   marked invalid and correlates as NaN; AFNI leaves it unscaled, giving about 0.
// * A series with a non-finite value is marked invalid up front. (Because two
//   series share one FFT, a NaN would otherwise leak into its partner.)
// ---------------------------------------------------------------------------

//! Seed correlation over a set of time series.

use crate::column::{ColumnData, ColumnRole, DataColumn};
use crate::dataset::{Dataset, DatasetKind};
use crate::error::{Error, Result};
use crate::numeric::ensure_finite;
use crate::signal::{
    bandpass_vectors, legendre_basis, normalize_l2, BandSpec, Detrend, MIN_BANDPASS_SAMPLES,
};
use crate::stat::{StatKind, StatSpec};

/// How to clean the series before correlating.
#[derive(Debug, Clone, PartialEq)]
pub struct InstaCorrOptions {
    /// Seconds between samples. Needed when `band` is set; taken from the dataset's
    /// time step when `None`.
    pub tr_seconds: Option<f64>,
    /// The pass band `(low_hz, high_hz)`, or `None` for no filtering. Frequencies
    /// above Nyquist are clipped to it, as in AFNI.
    pub band: Option<(f64, f64)>,
    /// Highest Legendre polynomial order removed (`-1` for none; SUMA's default is 2).
    pub polort: i32,
    /// SUMA's `normalize_dset`: clean and scale the series. When `false` nothing is
    /// done and the output is a plain dot product (not a correlation).
    pub normalize: bool,
    /// Extra regressors to remove (SUMA's `-ort`), each as long as the series.
    pub orts: Vec<Vec<f64>>,
}

impl Default for InstaCorrOptions {
    /// SUMA's defaults: band 0.01 to 0.1 Hz, `polort` 2, normalize on.
    fn default() -> Self {
        Self {
            tr_seconds: None,
            band: Some((0.01, 0.1)),
            polort: 2,
            normalize: true,
            orts: Vec::new(),
        }
    }
}

impl InstaCorrOptions {
    /// Check the options against a series length. Mirrors SUMA's checks: with
    /// `normalize` off nothing else matters; otherwise at least 9 samples,
    /// `polort >= -1`, `polort + 1 < samples - 3`, and a band needs a positive TR
    /// and `0 <= low < high`.
    pub fn validate(&self, sample_count: usize) -> Result<()> {
        if sample_count < 2 {
            return Err(Error::InvalidParameter {
                name: "time points".into(),
                reason: format!("seed correlation needs at least 2, got {sample_count}"),
            });
        }
        if !self.normalize {
            return Ok(());
        }
        if sample_count < MIN_BANDPASS_SAMPLES {
            return Err(Error::InvalidParameter {
                name: "time points".into(),
                reason: format!(
                    "cleaning the series needs at least {MIN_BANDPASS_SAMPLES} samples, got {sample_count}"
                ),
            });
        }
        if self.polort < -1 {
            return Err(Error::InvalidParameter {
                name: "polort".into(),
                reason: format!("{} is below -1", self.polort),
            });
        }
        let nref = (self.polort + 1) as usize;
        if nref + 3 >= sample_count {
            return Err(Error::InvalidParameter {
                name: "polort".into(),
                reason: format!(
                    "{nref} baseline regressors are too many for {sample_count} samples (need fewer than samples - 3)"
                ),
            });
        }
        for (i, ort) in self.orts.iter().enumerate() {
            if ort.len() != sample_count {
                return Err(Error::LengthMismatch {
                    what: format!("ort {i}"),
                    expected: sample_count,
                    found: ort.len(),
                });
            }
        }
        if let Some((low, high)) = self.band {
            let tr = self.tr_seconds.ok_or_else(|| Error::InvalidParameter {
                name: "TR".into(),
                reason: "a band needs the time step in seconds".into(),
            })?;
            ensure_finite("TR", tr)?;
            if tr <= 0.0 {
                return Err(Error::InvalidParameter {
                    name: "TR".into(),
                    reason: format!("{tr} is not positive"),
                });
            }
            ensure_finite("low cutoff", low)?;
            if low < 0.0 || high.is_nan() || high <= low {
                return Err(Error::InvalidParameter {
                    name: "band".into(),
                    reason: format!("needs 0 <= low < high, got {low} and {high}"),
                });
            }
        }
        Ok(())
    }

    /// The regressors to remove: Legendre `0..=polort` and the extras.
    fn regressors(&self, len: usize) -> Result<Vec<Vec<f64>>> {
        let mut orts = if self.polort >= 0 {
            legendre_basis((self.polort + 1) as usize, len)?
        } else {
            Vec::new()
        };
        orts.extend(self.orts.iter().cloned());
        Ok(orts)
    }

    fn band_spec(&self) -> Option<BandSpec> {
        self.band.map(|(fbot, ftop)| BandSpec {
            dt: self.tr_seconds.unwrap_or(1.0),
            fbot,
            ftop,
        })
    }
}

/// A set of cleaned, unit-length series ready for correlation. Once built, changing
/// the seed costs one dot product per series.
#[derive(Debug, Clone, PartialEq)]
pub struct PreparedSeries {
    /// `row_count * sample_count` values, row after row (32-bit like SUMA's).
    rows: Vec<f32>,
    valid: Vec<bool>,
    sample_count: usize,
    removed_dof: usize,
    options: InstaCorrOptions,
}

/// Build prepared series from plain rows (one `Vec` of samples per node).
pub fn prepare_rows(rows: &[Vec<f64>], options: &InstaCorrOptions) -> Result<PreparedSeries> {
    let sample_count = rows.first().map_or(0, Vec::len);
    options.validate(sample_count)?;
    if let Some(bad) = rows.iter().position(|r| r.len() != sample_count) {
        return Err(Error::LengthMismatch {
            what: format!("series {bad}"),
            expected: sample_count,
            found: rows[bad].len(),
        });
    }
    // Rows with a non-finite value are set aside first (zeroed, so they cannot
    // disturb the series they share an FFT with) and stay invalid afterwards.
    let mut valid: Vec<bool> = rows
        .iter()
        .map(|r| r.iter().all(|v| v.is_finite()))
        .collect();
    let mut work: Vec<Vec<f64>> = rows
        .iter()
        .zip(&valid)
        .map(|(r, &ok)| {
            if ok {
                r.clone()
            } else {
                vec![0.0; sample_count]
            }
        })
        .collect();

    let mut removed_dof = 0;
    if options.normalize {
        let orts = options.regressors(sample_count)?;
        let band = options.band_spec();
        let report = bandpass_vectors(&mut work, band.as_ref(), Detrend::Linear, &orts)?;
        removed_dof = report.removed_dof;
        for (row, ok) in work.iter_mut().zip(valid.iter_mut()) {
            if *ok && normalize_l2(row).is_none() {
                *ok = false; // no variance left: nothing to correlate
            }
        }
    }
    let mut flat = Vec::with_capacity(work.len() * sample_count);
    for (row, &ok) in work.iter().zip(&valid) {
        if ok {
            flat.extend(row.iter().map(|&v| v as f32));
        } else {
            flat.extend(std::iter::repeat(0.0_f32).take(sample_count));
        }
    }
    Ok(PreparedSeries {
        rows: flat,
        valid,
        sample_count,
        removed_dof,
        options: options.clone(),
    })
}

/// Build prepared series from a time-series dataset: the `TimePoint` columns (or, if
/// there are none, its numeric `Generic` columns), in order, one row per dataset row.
/// The dataset's time step is used when the options give no TR. Sparse datasets keep their row layout; map rows back with
/// [`Dataset::sample_for_row`].
pub fn prepare_dataset(dataset: &Dataset, options: &InstaCorrOptions) -> Result<PreparedSeries> {
    if dataset.kind() != &DatasetKind::TimeSeries {
        return Err(Error::InvalidParameter {
            name: "dataset kind".into(),
            reason: "seed correlation needs a time-series dataset".into(),
        });
    }
    // The time points are the `TimePoint` columns. Files do not mark them (a NIML time
    // series stores its columns as generic numbers), so when there are none, the
    // numeric `Generic` columns of the time-series dataset are used.
    let mut columns: Vec<&DataColumn> = dataset.columns_with_role(&ColumnRole::TimePoint).collect();
    if columns.is_empty() {
        columns = dataset
            .columns_with_role(&ColumnRole::Generic)
            .filter(|c| c.values().is_numeric())
            .collect();
    }
    if columns.len() < 2 {
        return Err(Error::InvalidParameter {
            name: "time points".into(),
            reason: format!(
                "the dataset has {} time-point columns, need 2",
                columns.len()
            ),
        });
    }
    let mut opts = options.clone();
    if opts.tr_seconds.is_none() {
        opts.tr_seconds = dataset.time_step_seconds();
    }
    let rows: Vec<Vec<f64>> = (0..dataset.row_count())
        .map(|r| {
            columns
                .iter()
                .map(|c| {
                    c.values()
                        .get_f64(r)
                        .ok_or_else(|| Error::InvalidParameter {
                            name: "time-point column".into(),
                            reason: format!("column '{}' is not numeric", c.label()),
                        })
                })
                .collect::<Result<Vec<f64>>>()
        })
        .collect::<Result<_>>()?;
    prepare_rows(&rows, &opts)
}

/// Correlations of one seed with every series.
#[derive(Debug, Clone, PartialEq)]
pub struct SeedCorrelation {
    /// One value per row; `NaN` for a row that is invalid (non-finite input or no
    /// variance). Plain dot products, not correlations, when `normalize` was off.
    pub values: Vec<f32>,
    /// Time points per series.
    pub samples: usize,
    /// Dimensions the cleaning removed from each series (AFNI's `ndet`).
    pub removed_dof: usize,
    /// The seed rows used (empty for an external seed).
    pub seed_rows: Vec<usize>,
    /// `Correl(samples, 1, removed_dof)` when these are correlations with positive
    /// degrees of freedom; `None` otherwise (`normalize` off, or the cleaning took
    /// out as many dimensions as there are samples).
    pub stat: Option<StatSpec>,
}

impl SeedCorrelation {
    /// Degrees of freedom of each correlation: `samples - 1 - removed_dof`, if
    /// positive.
    pub fn degrees_of_freedom(&self) -> Option<usize> {
        self.samples
            .checked_sub(1 + self.removed_dof)
            .filter(|&d| d > 0)
    }

    /// The correlations as a `Float32` statistic column carrying the `Correl` spec.
    pub fn to_column(&self, label: &str) -> Result<DataColumn> {
        Ok(DataColumn::new(
            label,
            ColumnRole::Statistic,
            ColumnData::Float32(self.values.clone()),
        )?
        .with_stat(self.stat.clone()))
    }

    /// The correlations as a dataset with the same domain and rows as `source`.
    pub fn to_dataset(&self, source: &Dataset, label: &str) -> Result<Dataset> {
        let dataset = Dataset::new(
            DatasetKind::Scalar,
            source.domain().clone(),
            source.map().clone(),
            vec![self.to_column(label)?],
        )?;
        Ok(dataset.with_parent_ids(source.parent_ids().clone()))
    }
}

impl PreparedSeries {
    /// Number of series.
    pub fn row_count(&self) -> usize {
        self.valid.len()
    }

    /// Time points per series.
    pub fn sample_count(&self) -> usize {
        self.sample_count
    }

    /// Dimensions removed from each series by the cleaning.
    pub fn removed_dof(&self) -> usize {
        self.removed_dof
    }

    /// Whether row `row` can be correlated.
    pub fn is_valid(&self, row: usize) -> bool {
        self.valid.get(row).copied().unwrap_or(false)
    }

    /// The cleaned, unit-length series of a row.
    pub fn row(&self, row: usize) -> Option<&[f32]> {
        (row < self.row_count()).then(|| &self.rows[row * self.sample_count..][..self.sample_count])
    }

    fn finish(&self, seed: &[f64], seed_rows: Vec<usize>) -> SeedCorrelation {
        let n = self.sample_count;
        let values: Vec<f32> = (0..self.row_count())
            .map(|r| {
                if !self.valid[r] {
                    return f32::NAN;
                }
                let row = &self.rows[r * n..][..n];
                row.iter()
                    .zip(seed)
                    .map(|(&a, &b)| f64::from(a) * b)
                    .sum::<f64>() as f32
            })
            .collect();
        let stat = (self.options.normalize && n > 1 + self.removed_dof).then(|| {
            StatSpec::new(
                StatKind::Correl,
                &[n as f64, 1.0, self.removed_dof as f64],
                0.0,
            )
        });
        SeedCorrelation {
            values,
            samples: n,
            removed_dof: self.removed_dof,
            seed_rows,
            stat,
        }
    }

    /// Correlate every series with row `seed_row`. Errors if the row does not exist
    /// or is invalid.
    pub fn correlate_row(&self, seed_row: usize) -> Result<SeedCorrelation> {
        self.correlate_rows(&[seed_row])
    }

    /// Correlate every series with the mean of several rows (an ROI seed): the
    /// prepared rows are averaged and the average is scaled to unit length again
    /// (when `normalize` is on). Invalid rows are skipped; it is an error if no row
    /// is usable or a row does not exist.
    pub fn correlate_rows(&self, seed_rows: &[usize]) -> Result<SeedCorrelation> {
        let n = self.sample_count;
        if let Some(&bad) = seed_rows.iter().find(|&&r| r >= self.row_count()) {
            return Err(Error::IndexOutOfRange {
                index: bad as i64,
                len: self.row_count(),
            });
        }
        let used: Vec<usize> = seed_rows
            .iter()
            .copied()
            .filter(|&r| self.valid[r])
            .collect();
        if used.is_empty() {
            return Err(Error::Empty("valid seed rows".into()));
        }
        let mut seed = vec![0.0_f64; n];
        for &r in &used {
            for (s, &v) in seed.iter_mut().zip(&self.rows[r * n..][..n]) {
                *s += f64::from(v);
            }
        }
        seed.iter_mut().for_each(|s| *s /= used.len() as f64);
        if self.options.normalize && normalize_l2(&mut seed).is_none() {
            return Err(Error::InvalidParameter {
                name: "seed".into(),
                reason: "the seed rows cancel out: their average has no variance".into(),
            });
        }
        Ok(self.finish(&seed, used))
    }

    /// Correlate every series with a seed that is not one of them (for example a
    /// series loaded from a file). It is cleaned like SUMA's single-vector path (a
    /// MEAN detrend, the same band and regressors) and scaled to unit length.
    pub fn correlate_external(&self, series: &[f64]) -> Result<SeedCorrelation> {
        let n = self.sample_count;
        if series.len() != n {
            return Err(Error::LengthMismatch {
                what: "seed series".into(),
                expected: n,
                found: series.len(),
            });
        }
        if series.iter().any(|v| !v.is_finite()) {
            return Err(Error::NonFinite {
                what: "seed series".into(),
                value: f64::NAN,
            });
        }
        let mut seed = vec![series.to_vec()];
        if self.options.normalize {
            let orts = self.options.regressors(n)?;
            let band = self.options.band_spec();
            bandpass_vectors(&mut seed, band.as_ref(), Detrend::Mean, &orts)?;
            if normalize_l2(&mut seed[0]).is_none() {
                return Err(Error::InvalidParameter {
                    name: "seed".into(),
                    reason: "the seed has no variance after cleaning".into(),
                });
            }
        }
        Ok(self.finish(&seed[0], Vec::new()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::{Domain, SurfaceDomain};

    /// Rows made from two latent signals plus noise: rows 0..3 follow the first, 3..6
    /// the second, so correlations are strong within a group and weak between.
    fn rows(n: usize) -> Vec<Vec<f64>> {
        let a: Vec<f64> = (0..n).map(|t| (t as f64 * 0.31).sin()).collect();
        let b: Vec<f64> = (0..n).map(|t| (t as f64 * 0.17 + 1.0).cos()).collect();
        let mut state = 12345_u64;
        let mut noise = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 33) as f64 / (1u64 << 31) as f64 - 0.5) * 0.2
        };
        (0..6)
            .map(|r| {
                let base = if r < 3 { &a } else { &b };
                (0..n)
                    .map(|t| {
                        10.0 + 0.05 * t as f64 + 3.0 * base[t] * (1.0 + r as f64 * 0.1) + noise()
                    })
                    .collect()
            })
            .collect()
    }

    fn opts() -> InstaCorrOptions {
        InstaCorrOptions {
            tr_seconds: Some(1.0),
            band: Some((0.02, 0.3)),
            polort: 2,
            normalize: true,
            orts: vec![],
        }
    }

    fn pearson(x: &[f64], y: &[f64]) -> f64 {
        let n = x.len() as f64;
        let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
        let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
        for (a, b) in x.iter().zip(y) {
            sxy += (a - mx) * (b - my);
            sxx += (a - mx) * (a - mx);
            syy += (b - my) * (b - my);
        }
        sxy / (sxx * syy).sqrt()
    }

    #[test]
    fn prepared_rows_have_unit_length_and_the_self_correlation_is_one() {
        let p = prepare_rows(&rows(120), &opts()).unwrap();
        assert_eq!((p.row_count(), p.sample_count()), (6, 120));
        for r in 0..6 {
            let len: f64 = p
                .row(r)
                .unwrap()
                .iter()
                .map(|&v| f64::from(v).powi(2))
                .sum();
            assert!((len - 1.0).abs() < 1e-6, "row {r}: {len}");
        }
        let c = p.correlate_row(0).unwrap();
        assert!((c.values[0] - 1.0).abs() < 1e-6);
        // Same-group rows correlate strongly, other-group rows weakly.
        assert!(c.values[1] > 0.9 && c.values[2] > 0.9);
        assert!(c.values[4].abs() < 0.5);
        assert!(c.values.iter().all(|v| (-1.0001..=1.0001).contains(v)));
    }

    #[test]
    fn without_a_band_or_polort_it_is_pearson_after_a_linear_detrend() {
        // No filter, no regressors: the cleaning is just the linear detrend, so the
        // result is the Pearson correlation of linearly detrended series.
        let o = InstaCorrOptions {
            band: None,
            polort: -1,
            ..opts()
        };
        let data = rows(60);
        let p = prepare_rows(&data, &o).unwrap();
        let c = p.correlate_row(0).unwrap();
        let detrended: Vec<Vec<f64>> = data
            .iter()
            .map(|r| {
                let mut v = r.clone();
                Detrend::Linear.apply(&mut v);
                v
            })
            .collect();
        for r in 0..6 {
            let want = pearson(&detrended[0], &detrended[r]);
            assert!((f64::from(c.values[r]) - want).abs() < 1e-5, "row {r}");
        }
        // Detrend removed 1 dimension.
        assert_eq!(c.removed_dof, 1);
    }

    #[test]
    fn the_statistic_carries_the_removed_degrees_of_freedom() {
        let p = prepare_rows(&rows(120), &opts()).unwrap();
        let c = p.correlate_row(0).unwrap();
        let stat = c.stat.clone().unwrap();
        assert_eq!(stat.kind, StatKind::Correl);
        assert_eq!(stat.params[0], 120.0);
        assert_eq!(stat.params[1], 1.0);
        assert_eq!(stat.params[2], c.removed_dof as f64);
        // 1 (linear) + 2 + (2*jbot - 1) + (2*(60 - jtop) - 1) + 3 orts, all positive.
        assert!(c.removed_dof > 20 && c.degrees_of_freedom().unwrap() < 100);
        // The column keeps the spec so p-values use the right degrees of freedom.
        let col = c.to_column("corr").unwrap();
        assert_eq!(col.stat().unwrap(), &stat);
    }

    #[test]
    fn the_raw_mode_is_a_dot_product_with_no_statistic() {
        let o = InstaCorrOptions {
            normalize: false,
            band: None,
            ..opts()
        };
        let data = rows(40);
        let p = prepare_rows(&data, &o).unwrap();
        let c = p.correlate_row(1).unwrap();
        let dot: f64 = data[1].iter().zip(&data[2]).map(|(a, b)| a * b).sum();
        assert!((f64::from(c.values[2]) - dot).abs() / dot < 1e-6);
        assert!(c.stat.is_none() && c.removed_dof == 0);
        // Raw mode even runs on short series and needs no TR.
        assert!(prepare_rows(&[vec![1.0, 2.0, 3.0]], &o).is_ok());
    }

    #[test]
    fn roi_seeds_average_prepared_rows() {
        let p = prepare_rows(&rows(120), &opts()).unwrap();
        let roi = p.correlate_rows(&[0, 1, 2]).unwrap();
        assert_eq!(roi.seed_rows, vec![0, 1, 2]);
        // The mean of the group correlates strongly with each member and weakly with
        // the other group, and is at least as good as the worst single seed.
        assert!(roi.values[..3].iter().all(|&v| v > 0.9));
        assert!(roi.values[3..].iter().all(|v| v.abs() < 0.6));
        // Explicit check of the definition: normalize(mean of prepared rows).
        let mut mean = vec![0.0; 120];
        for r in 0..3 {
            for (m, &v) in mean.iter_mut().zip(p.row(r).unwrap()) {
                *m += f64::from(v) / 3.0;
            }
        }
        normalize_l2(&mut mean);
        let want: f64 = p
            .row(4)
            .unwrap()
            .iter()
            .zip(&mean)
            .map(|(&a, b)| f64::from(a) * b)
            .sum();
        assert!((f64::from(roi.values[4]) - want).abs() < 1e-6);
        assert!(p.correlate_rows(&[]).is_err());
        assert!(p.correlate_rows(&[0, 99]).is_err());
    }

    #[test]
    fn external_seeds_are_cleaned_with_a_mean_detrend() {
        let data = rows(120);
        let p = prepare_rows(&data, &opts()).unwrap();
        // The raw series of row 0 as an outside seed correlates ~ like row 0 itself
        // (not exactly: it was cleaned by the single-vector path, a mean detrend).
        let ext = p.correlate_external(&data[0]).unwrap();
        let own = p.correlate_row(0).unwrap();
        for r in 0..6 {
            assert!((ext.values[r] - own.values[r]).abs() < 0.05, "row {r}");
        }
        assert!(ext.seed_rows.is_empty());
        assert!(p.correlate_external(&[1.0; 5]).is_err());
        assert!(p.correlate_external(&vec![f64::NAN; 120]).is_err());
        assert!(p.correlate_external(&vec![3.0; 120]).is_err()); // no variance
    }

    #[test]
    fn bad_rows_are_invalid_not_fatal_and_do_not_leak() {
        let mut data = rows(120);
        data[2][10] = f64::NAN;
        data[3] = vec![5.0; 120]; // constant: no variance after cleaning
        let p = prepare_rows(&data, &opts()).unwrap();
        assert!(!p.is_valid(2) && !p.is_valid(3));
        assert!(p.is_valid(0) && p.is_valid(1));
        let c = p.correlate_row(0).unwrap();
        assert!(c.values[2].is_nan() && c.values[3].is_nan());
        // The NaN row sits beside a good row in the FFT pairing: its partner is intact.
        let clean = prepare_rows(&rows(120), &opts()).unwrap();
        let reference = clean.correlate_row(0).unwrap();
        assert!((c.values[1] - reference.values[1]).abs() < 1e-6);
        assert!(p.correlate_row(2).is_err());
    }

    #[test]
    fn options_are_validated_like_suma() {
        let data = rows(120);
        let no_tr = InstaCorrOptions {
            tr_seconds: None,
            ..opts()
        };
        assert!(prepare_rows(&data, &no_tr).is_err());
        let backwards = InstaCorrOptions {
            band: Some((0.2, 0.1)),
            ..opts()
        };
        assert!(prepare_rows(&data, &backwards).is_err());
        let bad_polort = InstaCorrOptions {
            polort: 200,
            ..opts()
        };
        assert!(prepare_rows(&data, &bad_polort).is_err());
        assert!(prepare_rows(&rows(8), &opts()).is_err()); // fewer than 9 samples
        let short_ort = InstaCorrOptions {
            orts: vec![vec![1.0; 5]],
            ..opts()
        };
        assert!(prepare_rows(&data, &short_ort).is_err());
        assert!(prepare_rows(&[vec![1.0; 20], vec![1.0; 19]], &opts()).is_err());
        assert!(prepare_rows(&[], &opts()).is_err());
    }

    #[test]
    fn datasets_in_and_correlation_datasets_out() {
        use crate::column::{ColumnData, DataColumn};
        let n = 24;
        let data = rows(n);
        let columns: Vec<DataColumn> = (0..n)
            .map(|t| {
                DataColumn::new(
                    format!("t{t}"),
                    ColumnRole::TimePoint,
                    ColumnData::Float64(data.iter().map(|r| r[t]).collect()),
                )
                .unwrap()
            })
            .collect();
        let domain = Domain::Surface(SurfaceDomain::new(None, 6).unwrap());
        let ds = Dataset::dense(DatasetKind::TimeSeries, domain, columns)
            .unwrap()
            .with_time_step_seconds(Some(2.0))
            .unwrap();
        // The TR comes from the dataset.
        let o = InstaCorrOptions {
            tr_seconds: None,
            band: Some((0.02, 0.2)),
            polort: 1,
            ..opts()
        };
        let p = prepare_dataset(&ds, &o).unwrap();
        assert_eq!((p.row_count(), p.sample_count()), (6, n));
        let c = p.correlate_row(0).unwrap();
        let out = c.to_dataset(&ds, "seed 0").unwrap();
        assert_eq!(out.row_count(), 6);
        assert_eq!(out.columns()[0].stat().unwrap().kind, StatKind::Correl);
        // Unmarked columns work too: a file stores its time points as generic numbers.
        let generic: Vec<DataColumn> = (0..n)
            .map(|t| {
                DataColumn::new(
                    format!("g{t}"),
                    ColumnRole::Generic,
                    ColumnData::Float64(data.iter().map(|r| r[t]).collect()),
                )
                .unwrap()
            })
            .collect();
        let from_file = Dataset::dense(
            DatasetKind::TimeSeries,
            Domain::Surface(SurfaceDomain::new(None, 6).unwrap()),
            generic,
        )
        .unwrap()
        .with_time_step_seconds(Some(2.0))
        .unwrap();
        let q = prepare_dataset(&from_file, &o).unwrap();
        assert_eq!(q.correlate_row(0).unwrap().values, c.values);
        // A scalar dataset is refused.
        let scalar = Dataset::dense(
            DatasetKind::Scalar,
            Domain::Surface(SurfaceDomain::new(None, 2).unwrap()),
            vec![DataColumn::new(
                "x",
                ColumnRole::Generic,
                ColumnData::Float64(vec![1.0, 2.0]),
            )
            .unwrap()],
        )
        .unwrap();
        assert!(prepare_dataset(&scalar, &o).is_err());
    }
}
