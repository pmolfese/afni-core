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
// False discovery rate (FDR) support, in two halves:
//
//  1. USING a stored curve. AFNI stores, next to each statistical sub-brick,
//     an "FDR curve" mapping a statistic threshold to z(q) (the z-score
//     equivalent of the FDR q-value) and an optional "MDF curve" giving the
//     missed-detection fraction versus log10(p). Here: threshold -> q,
//     q -> threshold, and p -> MDF.
//  2. BUILDING a curve from data, a port of AFNI's `mri_fdrize` and
//     `mri_fdr_curve`: convert statistics to p, sort, step-up to q, estimate the
//     number of true positives, optionally rescale, and tabulate 101 points.
//  Also a plain Benjamini-Hochberg helper for callers with raw p-values.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `curve.rs` holds the validated table and AFNI's exact interpolation;
//   `stats.rs` supplies p-values; `special.rs` the normal quantile.
// * `column.rs` stores a column's FDR/MDF curves; this file computes with them.
// * Phase 5 (thresholding) will call `threshold_for_q` so a viewer can show
//   "q <= 0.05" as a statistic threshold.
//
// P-VALUES AND Q-VALUES ARE DIFFERENT THINGS and are kept apart by type:
// `Probability` (stats.rs) is a p-value; `QValue` (here) is a q-value.
//
// FIDELITY TO AFNI
//
// AFNI computes in single precision and has quirks. This port reproduces the
// arithmetic where it affects results (stored arrays are `f32`, p-values are
// floored at 1e-15 and rounded to `f32`, the m1 histogram is binned in `f32`),
// and documents each quirk where it is copied. Differences from the C code are
// listed in the roadmap discovery log.
// ---------------------------------------------------------------------------

//! FDR curves, q-values, and missed-detection fractions.

use crate::curve::ThresholdCurve;
use crate::error::{Error, Result};
use crate::numeric::ensure_finite;
use crate::special::{normal_ln_tails, normal_upper_quantile_ln};
use crate::stat::{StatKind, StatSpec};
use crate::stats::{PValueEvaluator, Tail};

/// An FDR q-value in `[0, 1]`: the smallest false discovery rate at which a
/// result would be declared significant. Distinct from a p-value
/// ([`Probability`](crate::stats::Probability)).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct QValue(f64);

impl QValue {
    /// A q-value; must be finite and in `[0, 1]`.
    pub fn new(q: f64) -> Result<Self> {
        ensure_finite("q-value", q)?;
        if (0.0..=1.0).contains(&q) {
            Ok(Self(q))
        } else {
            Err(Error::InvalidParameter {
                name: "q-value".into(),
                reason: format!("{q} is outside [0, 1]"),
            })
        }
    }

    /// The value.
    pub fn get(self) -> f64 {
        self.0
    }
}

/// Largest z(q) AFNI writes into an FDR curve or z-score map (`ZTOP`).
pub const Z_TOP: f64 = 9.0;
/// Smallest q AFNI distinguishes (`QBOT`); smaller q map to [`Z_TOP`].
pub const Q_BOTTOM: f64 = 2.25718e-19;
/// p-values below this are raised to it before FDR processing (`PBOT`).
pub const P_BOTTOM: f64 = 1.0e-15;
/// p-values at or above this are not processed (`PMAX`).
const P_MAX: f32 = 0.9999;
/// q level below which the missed-detection curve is attempted (`QSTHRESH`).
const Q_SMALL: f64 = 0.15;
/// Number of samples in a generated FDR curve (`NCURV`).
const CURVE_POINTS: usize = 101;

/// The z-score whose two-sided normal tail is `q`: AFNI's `QTOZ(q) = qginv(q/2)`.
/// `q` must be in `(0, 1]`.
pub fn z_for_q(q: QValue) -> Result<f64> {
    let q = q.get();
    if q <= 0.0 {
        return Err(Error::NoSolution(
            "q = 0 corresponds to an infinite z-score".into(),
        ));
    }
    if q >= 1.0 {
        return Ok(0.0);
    }
    normal_upper_quantile_ln((0.5 * q).ln())
        .ok_or_else(|| Error::NoSolution("normal quantile did not converge".into()))
}

/// The two-sided q for a z-score: AFNI's `ZTOQ(z) = 2 qg(z)`, for `z >= 0`.
pub fn q_for_z(z: f64) -> Result<QValue> {
    ensure_finite("z", z)?;
    if z <= 0.0 {
        return QValue::new(1.0);
    }
    let (_, ln_upper) = normal_ln_tails(z)
        .ok_or_else(|| Error::NoSolution("normal tail did not converge".into()))?;
    QValue::new((2.0 * ln_upper.exp()).min(1.0))
}

// ---------------------------------------------------------------------------
// Using a stored FDR curve
// ---------------------------------------------------------------------------

/// The q-value at statistic `threshold` on an FDR `curve`
/// (`thd_fdrcurve.c: THD_fdrcurve_zval`, then `fdrval`).
///
/// The curve gives `z(q)` for a threshold on the *absolute* statistic, so pass
/// `|statistic|`. If the interpolated `z` is not positive (a threshold below the
/// curve's range) the q-value is 1, as in `fdrval`.
pub fn q_value_for_threshold(curve: &ThresholdCurve, threshold: f64) -> Result<QValue> {
    ensure_finite("threshold", threshold)?;
    q_for_z(curve.interpolate(threshold))
}

/// The statistic threshold at which an FDR `curve` reaches q-value `q`
/// (`thd_fdrcurve.c: THD_fdrcurve_zqtot`, as used by `fdrval -qinput`).
///
/// `max_abs` is the largest absolute value of the sub-brick's data. AFNI uses it
/// when `q` is smaller than the curve can express: the threshold is then at
/// least `max_abs * 1.000002`, i.e. just above every sample (so nothing passes).
/// Pass `None` if unknown.
///
/// Behaviour at the edges (AFNI's, kept):
///
/// * `q` too small for the curve: `x0 + dx * len` (one step past the last
///   sample), or `max_abs * 1.000002` if that is larger;
/// * `q` too large for the curve: `0`;
/// * otherwise [`ThresholdCurve::inverse_interpolate`].
///
/// `q` must be in `(0, 1]` (AFNI's `fdrval` silently clamps it to
/// `[1e-9, 0.99999]`; here `q = 0` is an error and `q = 1` gives threshold 0).
pub fn threshold_for_q(curve: &ThresholdCurve, q: QValue, max_abs: Option<f64>) -> Result<f64> {
    let z = z_for_q(q)?;
    let samples = curve.samples();
    let last = *samples.last().expect("curves have samples");
    Ok(if z > last {
        let mut threshold = curve.x0() + curve.dx() * samples.len() as f64;
        if let Some(mab) = max_abs {
            if mab >= threshold {
                threshold = mab * 1.000002;
            }
        }
        threshold
    } else if z < samples[0] {
        0.0
    } else {
        curve.inverse_interpolate(z)
    })
}

/// The smallest q an FDR `curve` can express (`DSET_BRICK_FDRMIN`): the q of
/// its largest stored z. A large value means few true detections.
pub fn minimum_q(curve: &ThresholdCurve) -> Result<QValue> {
    q_for_z(*curve.samples().last().expect("curves have samples"))
}

/// The missed-detection fraction at p-value `p` on an MDF `curve`
/// (`thd_fdrcurve.c: THD_mdfcurve_mval`). The curve's axis is `log10(p)`, NOT a
/// statistic like an FDR curve's.
///
/// `p <= 0` gives 0.999 and `p >= 1` gives 0 (AFNI's edge values).
pub fn missed_detection_fraction(curve: &ThresholdCurve, p: f64) -> Result<f64> {
    ensure_finite("p-value", p)?;
    Ok(if p <= 0.0 {
        0.999
    } else if p >= 1.0 {
        0.0
    } else {
        curve.interpolate(p.log10())
    })
}

// ---------------------------------------------------------------------------
// Building FDR curves (port of mri_fdrize / mri_fdr_curve)
// ---------------------------------------------------------------------------

/// What `fdrize` writes for each sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FdrOutput {
    /// z(q), so that larger means more significant (AFNI's default; what
    /// `3dFDR` writes).
    #[default]
    ZScore,
    /// The q-value itself.
    QValue,
}

/// Options for [`fdrize`] and [`fdr_curves`].
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct FdrOptions {
    /// Tail used to turn each statistic into a p-value. `None` uses AFNI's rule
    /// for the classic statistics (two-sided for t, z and signed correlation;
    /// upper tail for the rest) and is an error for any other kind.
    pub tail: Option<Tail>,
    /// Correct for arbitrary dependence between tests by inflating the test
    /// count by `ln(n) + 0.5772157 + 0.5/n` (AFNI `flags & 2`, the environment
    /// variable `AFNI_NON_INDEPENDENT_FDR`). Not needed for typical fMRI maps.
    pub arbitrary_dependence: bool,
    /// Skip AFNI's rescaling of q by the estimated fraction of true nulls
    /// (the inverse of leaving `AFNI_DONT_ADJUST_FDR` unset). Default `false`
    /// (AFNI's default: the adjustment is applied).
    pub skip_true_positive_adjustment: bool,
    /// z(q) or q for the per-sample output.
    pub output: FdrOutput,
    /// How samples with exactly equal z(q) are ordered when building a curve
    /// (see [`fdr_curves`]). `false` (the default) reproduces AFNI's unstable
    /// quicksort exactly, so curves match AFNI's point for point. `true` orders
    /// ties by increasing statistic, which is reproducible without reference to
    /// any particular sort and starts the curve at the smallest tied statistic.
    pub deterministic_ties: bool,
}

/// Result of [`fdrize`].
#[derive(Debug, Clone, PartialEq)]
pub struct Fdrized {
    /// One value per input sample: z(q) or q. Samples that were ignored
    /// (masked, zero, non-finite, or p >= 0.9999) get 0 for z(q) and 1 for q.
    pub values: Vec<f32>,
    /// How many samples took part in the FDR calculation.
    pub used: usize,
    /// AFNI's estimate of the number of truly active samples (`m1`), if it could
    /// be made (needs at least 233 samples with 160 in the p = 0.1..0.95 range).
    pub true_positives: Option<f64>,
    /// Missed-detection fraction versus log10(p), when it could be built.
    pub missed_detection: Option<ThresholdCurve>,
}

/// AFNI's tail rule (`THD_stat_is_2sided`) for the classic statistics.
pub fn afni_tail(spec: &StatSpec) -> Option<Tail> {
    match spec.kind {
        StatKind::Ttest | StatKind::Zscore => Some(Tail::TwoSided),
        // Signed correlation is two-sided; a multiple correlation (nfit > 1)
        // is |R| >= 0 and one-sided.
        StatKind::Correl if spec.params.get(1) == Some(&1.0) => Some(Tail::TwoSided),
        StatKind::Correl
        | StatKind::Ftest
        | StatKind::Chisq
        | StatKind::Beta
        | StatKind::Binom
        | StatKind::Gamma
        | StatKind::Poisson => Some(Tail::Upper),
        _ => None,
    }
}

/// Estimate of the number of truly active samples, from the p-values
/// (`estimate_m1`): histogram the p-values in 16 bins of width 0.05 starting at
/// 0.15, take the flat part of the sorted histogram as the null density, and
/// scale up. `None` when there are too few samples or too few p in range.
///
/// AFNI quirk, kept: the bin index is `(int)((p - 0.15) * 20)`, which truncates
/// toward zero, so p-values in `(0.10, 0.15)` land in bin 0 instead of being
/// excluded.
fn estimate_true_positives(p_sorted: &[f32]) -> Option<f32> {
    let nq = p_sorted.len();
    if nq < 233 {
        return None;
    }
    let mut hist = [0_i32; 16];
    let mut in_range = 0;
    for &p in p_sorted {
        let bin = ((p - 0.15_f32) * 20.0_f32) as i32;
        if (0..=15).contains(&bin) {
            hist[bin as usize] += 1;
            in_range += 1;
        }
    }
    if in_range < 160 {
        return None;
    }
    hist.sort_unstable();
    let h = |i: usize| hist[i] as f32;
    let n = nq as f32;
    // Two estimates of the null count from the middle bins; keep the smaller m1.
    let ma = (n - 20.0_f32 * (h(6) + 2.0 * h(7) + 2.0 * h(8) + h(9)) / 6.0_f32) as i32;
    let mb = (n - 20.0_f32 * (h(5) + 2.0 * h(6) + 2.0 * h(7) + 2.0 * h(8) + 2.0 * h(9) + h(10))
        / 10.0_f32) as i32;
    Some(ma.min(mb) as f32)
}

/// FDR-ize a set of p-values (AFNI `mri_fdrize` with p-value input).
///
/// `pvalues` outside `[0, 1]`, NaN, or masked out (`mask[i] == false`) are
/// ignored. Needs at least 20 usable values.
pub fn fdrize_p_values(
    pvalues: &[f32],
    mask: Option<&[bool]>,
    options: &FdrOptions,
) -> Result<Fdrized> {
    check_mask(pvalues.len(), mask)?;
    let p: Vec<f32> = pvalues
        .iter()
        .enumerate()
        .map(|(i, &v)| {
            let masked = mask.is_some_and(|m| !m[i]);
            // NaN fails both comparisons and is dropped as "not reasonable".
            if masked || !(0.0..=1.0).contains(&v) {
                1.0
            } else {
                v
            }
        })
        .collect();
    fdrize_core(&p, options)
}

/// FDR-ize a set of statistics (AFNI `mri_fdrize` with statistic input).
///
/// Each value is converted to a p-value with `spec` under `options.tail` (AFNI's
/// rule by default). Zero statistics, non-finite statistics, and masked samples
/// count as p = 1 and are ignored, as in AFNI. Statistics are single precision
/// because AFNI's are; convert deliberately with
/// [`narrow_to_f32`](crate::numeric::narrow_to_f32).
pub fn fdrize(
    spec: &StatSpec,
    statistics: &[f32],
    mask: Option<&[bool]>,
    options: &FdrOptions,
) -> Result<Fdrized> {
    check_mask(statistics.len(), mask)?;
    let evaluator = evaluator_for(spec, options)?;
    let p: Vec<f32> = statistics
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            let masked = mask.is_some_and(|m| !m[i]);
            if x == 0.0 || !x.is_finite() || masked {
                return Ok(1.0_f32);
            }
            // AFNI converts |statistic| and stores the p-value as a float.
            Ok(evaluator.evaluate(f64::from(x.abs()))?.p() as f32)
        })
        .collect::<Result<_>>()?;
    fdrize_core(&p, options)
}

fn check_mask(len: usize, mask: Option<&[bool]>) -> Result<()> {
    match mask {
        Some(m) if m.len() != len => Err(Error::LengthMismatch {
            what: "FDR mask".into(),
            expected: len,
            found: m.len(),
        }),
        _ => Ok(()),
    }
}

fn evaluator_for(spec: &StatSpec, options: &FdrOptions) -> Result<PValueEvaluator> {
    let tail = match options.tail.or_else(|| afni_tail(spec)) {
        Some(t) => t,
        None => {
            return Err(Error::InvalidParameter {
                name: "tail".into(),
                reason: format!(
                    "AFNI defines FDR only for the classic statistics (codes 2-10); for {} \
                     choose a tail explicitly",
                    spec.kind.name()
                ),
            })
        }
    };
    PValueEvaluator::new(spec, tail)
}

/// The shared algorithm. `p` holds one p-value per sample, with 1.0 meaning
/// "ignore". Mirrors `mri_fdrize` step by step (see the C source).
fn fdrize_core(p: &[f32], options: &FdrOptions) -> Result<Fdrized> {
    let want_z = options.output == FdrOutput::ZScore;
    let ignored_value: f32 = if want_z { 0.0 } else { 1.0 };

    // Reasonable p-values only (p < PMAX), raised to the floor PBOT.
    let mut sorted: Vec<(f32, usize)> = Vec::new();
    for (i, &v) in p.iter().enumerate() {
        if (0.0..P_MAX).contains(&v) {
            sorted.push((v.max(P_BOTTOM as f32), i));
        }
    }
    let nq = sorted.len();
    if nq <= 19 {
        // AFNI warns "will not process only N values (min=20)" and returns 0.
        return Err(Error::InvalidParameter {
            name: "samples".into(),
            reason: format!("FDR needs at least 20 usable p-values, found {nq}"),
        });
    }
    // Stable sort by p; ties give equal q, so their order does not matter.
    sorted.sort_by(|a, b| a.0.total_cmp(&b.0));

    let mut values = vec![ignored_value; p.len()];

    // Step-up: scan from the largest p down, q_j = min(q_{j+1}, n p_j / (j + 1)).
    let mut nthr = nq as f32;
    if options.arbitrary_dependence && nthr > 1.0 {
        nthr *= nthr.ln() + 0.577_215_7_f32 + 0.5_f32 / nthr;
    }
    let mut q_min = 1.0_f64;
    let mut first_small: usize = 0; // index of the first q <= 0.15 met scanning down (0 = none)
    let mut q_raw = vec![0.0_f32; nq];
    for j in (0..nq).rev() {
        // The product is a single-precision multiply in AFNI, then divided in double.
        let mut q = f64::from(nthr * sorted[j].0) / (j as f64 + 1.0);
        if q > q_min {
            q = q_min;
        } else {
            q_min = q;
        }
        if first_small == 0 && q <= Q_SMALL {
            first_small = j;
        }
        q_raw[j] = q as f32;
    }

    // Estimate the number of true positives (m1) and rescale q by it.
    let p_only: Vec<f32> = sorted.iter().map(|s| s.0).collect();
    let mut m1: f32 = 0.0;
    if first_small != 0 && sorted[0].0 > 0.0 {
        if let Some(m) = estimate_true_positives(&p_only) {
            m1 = m;
        }
    }
    let m1_estimate = (m1 != 0.0).then_some(f64::from(m1));
    if m1 > 0.0 && !options.skip_true_positive_adjustment {
        let mut factor = (nq as f32 - m1) / nq as f32;
        if factor < 0.5 {
            factor = 0.25 + factor * factor;
        }
        for q in &mut q_raw {
            *q *= factor;
        }
    }

    // Write q or z(q) back to each sample's position.
    for (j, &(_, index)) in sorted.iter().enumerate() {
        values[index] = if want_z {
            let q = f64::from(q_raw[j]);
            let z = if q < Q_BOTTOM {
                Z_TOP
            } else if q >= 1.0 {
                0.0
            } else {
                z_for_q(QValue::new(q)?)?
            };
            z as f32
        } else {
            q_raw[j]
        };
    }

    // Missed-detection fraction, only when m1 is comfortably positive.
    let missed_detection = if m1 > 8.0 {
        missed_detection_curve(&sorted, nthr, m1)
    } else {
        None
    };

    Ok(Fdrized {
        values,
        used: nq,
        true_positives: m1_estimate,
        missed_detection,
    })
}

/// Missed-detection fraction versus `log10(p)` (the second half of `mri_fdrize`).
///
/// Reasoning (AFNI's): at a threshold passing `j + 1` samples about
/// `q (j + 1)` are false positives, so `(1 - q)(j + 1)` are true detections and
/// `1 - (1 - q)(j + 1) / m1` is the fraction of the `m1` true positives missed.
/// It depends entirely on the quality of the `m1` estimate.
fn missed_detection_curve(sorted: &[(f32, usize)], nthr: f32, m1: f32) -> Option<ThresholdCurve> {
    let nq = sorted.len();
    let qq = |j: usize| sorted[j].0;
    let mut mdf = vec![0.0_f32; nq];
    let mut q_min = 1.0_f64;
    for j in (0..nq).rev() {
        // Not adjusted by the true-positive factor, as in AFNI.
        let mut q = f64::from(nthr * qq(j)) / (j as f64 + 1.0);
        if q > q_min {
            q = q_min;
        } else {
            q_min = q;
        }
        if q < Q_BOTTOM {
            q = Q_BOTTOM;
        }
        let m = 1.0 - (1.0 - q) * (j as f64 + 1.0) / f64::from(m1);
        mdf[j] = (m as f32).clamp(0.0, 1.0);
    }
    // Missed detections can only decrease as p increases.
    for j in 1..nq {
        if mdf[j] > mdf[j - 1] {
            mdf[j] = mdf[j - 1];
        }
    }
    // Stretch so the curve reaches 0 at p -> 1 ("cheapo trick" in AFNI).
    let last = mdf[nq - 1];
    if last > 0.0 {
        let alpha = 1.0_f32 / (1.0_f32 - last);
        let beta = alpha * last;
        for m in &mut mdf {
            *m = alpha * *m - beta;
        }
        mdf[nq - 1] = 0.0;
    }
    // Last nonzero entry.
    let mut j = nq - 2;
    while j > 0 && mdf[j] == 0.0 {
        j -= 1;
    }
    if j <= 1 || qq(j + 1) <= qq(0) {
        return None; // all zero, or p constant: no usable curve
    }
    let jtop = j + 1; // mdf[jtop] is zero
                      // First entry below 99.9%.
    let mut jj = 1;
    while jj < jtop && mdf[jj] >= 0.999 {
        jj += 1;
    }
    let jbot = if jj + 9 < jtop { jj - 1 } else { 0 };

    let span = (f64::from(qq(jtop) / qq(jbot)).log10()) as f32; // positive
    let npp = (0.99_f32 + 5.0_f32 * span) as i32;
    if npp < 3 {
        return None;
    }
    let npp = npp as usize;
    let dpl = span / (npp - 1) as f32;
    let x0 = f64::from(qq(jbot)).log10() as f32;
    let mut samples = vec![0.0_f32; npp];
    samples[0] = mdf[jbot];
    let mut jj = jbot;
    for (k, slot) in samples.iter_mut().enumerate().take(npp - 1).skip(1) {
        let pl = x0 + k as f32 * dpl; // log10(p) of this grid point
        let pv = 10.0_f32.powf(pl);
        while jj < jtop && qq(jj) < pv {
            jj += 1;
        }
        // Linear interpolation of mdf in log10(p) between samples jj-1 and jj.
        let p1 = f64::from(qq(jj - 1)).log10() as f32;
        let p2 = f64::from(qq(jj)).log10() as f32;
        let pf = (pl - p1) / (p2 - p1);
        *slot = pf * mdf[jj] + (1.0_f32 - pf) * mdf[jj - 1];
    }
    samples[npp - 1] = 0.0;
    ThresholdCurve::new(
        f64::from(x0),
        f64::from(dpl),
        samples.into_iter().map(f64::from).collect(),
    )
    .ok()
}

/// An FDR curve, its optional MDF curve, and how many samples produced them.
#[derive(Debug, Clone, PartialEq)]
pub struct FdrCurves {
    /// z(q) versus the absolute statistic, 101 points (AFNI `FDRCURVE_*`).
    pub fdr: ThresholdCurve,
    /// Missed-detection fraction versus log10(p) (AFNI `MDFCURVE_*`), if built.
    pub mdf: Option<ThresholdCurve>,
    /// Samples that took part.
    pub used: usize,
    /// AFNI's estimate of the number of truly active samples, if made.
    pub true_positives: Option<f64>,
}

/// Build the FDR (and MDF) curves for one statistical sub-brick: a port of
/// `mri_fdr_curve`.
///
/// The statistics are FDR-ized to z(q); samples with a meaningful (positive)
/// z(q) are sorted by z, and the curve is tabulated at 101 evenly spaced
/// absolute-statistic values from the smallest to the point where z(q) first
/// reaches [`Z_TOP`], interpolating z linearly in the statistic.
///
/// `options.output` is ignored (z(q) is always used). Fails if too few samples
/// (fewer than 20 usable, or fewer than 9 with positive z) or if no sample
/// reaches a z below the cap, as AFNI returns no curve then.
pub fn fdr_curves(
    spec: &StatSpec,
    statistics: &[f32],
    mask: Option<&[bool]>,
    options: &FdrOptions,
) -> Result<FdrCurves> {
    let z_options = FdrOptions {
        output: FdrOutput::ZScore,
        ..*options
    };
    let fdrized = fdrize(spec, statistics, mask, &z_options)?;
    if fdrized.used < 9 {
        return Err(Error::InvalidParameter {
            name: "samples".into(),
            reason: "too few samples to build an FDR curve".into(),
        });
    }
    // (z(q), |statistic|) for samples with a meaningful z.
    let mut pairs: Vec<(f32, f32)> = fdrized
        .values
        .iter()
        .zip(statistics)
        .filter(|(&z, _)| z > 0.0)
        .map(|(&z, &t)| (z, t.abs()))
        .collect();
    let nq = pairs.len();
    if nq < 9 {
        return Err(Error::InvalidParameter {
            name: "samples".into(),
            reason: format!("only {nq} samples have a positive z(q); need 9"),
        });
    }
    let (z, t) = if options.deterministic_ties {
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        (
            pairs.iter().map(|p| p.0).collect::<Vec<f32>>(),
            pairs.iter().map(|p| p.1).collect::<Vec<f32>>(),
        )
    } else {
        // AFNI's qsort_floatfloat, keyed on z with the statistic as payload.
        let mut z: Vec<f32> = pairs.iter().map(|p| p.0).collect();
        let mut t: Vec<f32> = pairs.iter().map(|p| p.1).collect();
        afni_qsort_float_float(&mut z, &mut t);
        (z, t)
    };

    // Largest z(q) not beyond the cap: the curve ends one past the last capped value.
    let mut klast = nq - 1;
    while klast > 0 && f64::from(z[klast]) >= Z_TOP {
        klast -= 1;
    }
    if klast == 0 {
        return Err(Error::NoSolution(
            "every sample is at the z(q) cap; no curve can be built".into(),
        ));
    }
    if klast < nq - 1 {
        klast += 1;
    }

    let (t_bottom, t_top) = (t[0], t[klast]);
    let dt = (t_top - t_bottom) / (CURVE_POINTS - 1) as f32;
    let mut samples = vec![0.0_f32; CURVE_POINTS];
    samples[0] = z[0];
    let mut j = 1;
    for (i, slot) in samples
        .iter_mut()
        .enumerate()
        .take(CURVE_POINTS - 1)
        .skip(1)
    {
        let tt = t_bottom + i as f32 * dt;
        while j < nq && t[j] < tt {
            j += 1;
        }
        let (t1, t2) = (t[j - 1], t[j]);
        // Linear interpolation of z in the statistic. AFNI divides by t2 - t1,
        // which is 0 only for a degenerate curve (all statistics equal); guard it.
        let frac = if t2 > t1 { (tt - t1) / (t2 - t1) } else { 0.0 };
        *slot = frac * z[j] + (1.0_f32 - frac) * z[j - 1];
    }
    samples[CURVE_POINTS - 1] = z[klast];
    let fdr = ThresholdCurve::new(
        f64::from(t_bottom),
        f64::from(dt),
        samples.into_iter().map(f64::from).collect(),
    )?;
    Ok(FdrCurves {
        fdr,
        mdf: fdrized.missed_detection,
        used: fdrized.used,
        true_positives: fdrized.true_positives,
    })
}

// ---------------------------------------------------------------------------
// AFNI's sort (cs_sort_ff.c), ported so curves match when z(q) values tie
// ---------------------------------------------------------------------------
//
// Many samples can share exactly the same z(q) (for example every sample on the
// flat top of the step-up procedure). A quicksort does not keep tied entries in
// input order, so WHICH statistic ends up first among the tied values depends on
// the algorithm. To reproduce AFNI's curve exactly this is its algorithm, line
// for line: median-of-three quicksort with an explicit stack, partitions of 10
// or fewer left for a final insertion sort.

/// Swap element `i` and `j` in both arrays.
fn swap_pair(a: &mut [f32], ia: &mut [f32], i: usize, j: usize) {
    a.swap(i, j);
    ia.swap(i, j);
}

/// Sort `a` ascending, carrying `ia` along, exactly as AFNI's
/// `qsort_floatfloat` does.
fn afni_qsort_float_float(a: &mut [f32], ia: &mut [f32]) {
    debug_assert_eq!(a.len(), ia.len());
    const CUTOFF: isize = 10;
    let n = a.len() as isize;
    if n >= CUTOFF.max(3) {
        let mut stack: Vec<isize> = vec![0, n - 1];
        while stack.len() >= 2 {
            let right = stack.pop().expect("stack has right");
            let left = stack.pop().expect("stack has left");
            let (l, r) = (left as usize, right as usize);
            let mid = ((left + right) / 2) as usize;
            // Order the left, middle and right entries; the middle is the pivot.
            if a[l] > a[mid] {
                swap_pair(a, ia, l, mid);
            }
            if a[l] > a[r] {
                swap_pair(a, ia, l, r);
            }
            if a[mid] > a[r] {
                swap_pair(a, ia, r, mid);
            }
            let pivot = a[mid];
            a[mid] = a[r];
            let pivot_payload = ia[mid];
            ia[mid] = ia[r];
            // Partition: scan in from both ends, swapping out-of-place pairs.
            let (mut i, mut j) = (left, right);
            loop {
                loop {
                    i += 1;
                    if a[i as usize] >= pivot {
                        break;
                    }
                }
                loop {
                    j -= 1;
                    if a[j as usize] <= pivot {
                        break;
                    }
                }
                if j <= i {
                    break;
                }
                swap_pair(a, ia, i as usize, j as usize);
            }
            // Restore the pivot.
            a[r] = a[i as usize];
            a[i as usize] = pivot;
            ia[r] = ia[i as usize];
            ia[i as usize] = pivot_payload;
            // Push the big-enough sub-ranges, shorter one to be handled first.
            let mut pushed = 0;
            if i - left > CUTOFF {
                stack.push(left);
                stack.push(i - 1);
                pushed += 1;
            }
            if right - i > CUTOFF {
                stack.push(i + 1);
                stack.push(right);
                pushed += 1;
            }
            if pushed == 2 {
                let m = stack.len();
                if stack[m - 3] - stack[m - 4] > stack[m - 1] - stack[m - 2] {
                    stack.swap(m - 4, m - 2);
                    stack.swap(m - 3, m - 1);
                }
            }
        }
    }
    // Insertion sort finishes the nearly sorted array (and handles short ones).
    for j in 1..a.len() {
        if a[j] < a[j - 1] {
            let (key, payload) = (a[j], ia[j]);
            let mut p = j;
            loop {
                a[p] = a[p - 1];
                ia[p] = ia[p - 1];
                p -= 1;
                if !(p > 0 && key < a[p - 1]) {
                    break;
                }
            }
            a[p] = key;
            ia[p] = payload;
        }
    }
}

// ---------------------------------------------------------------------------
// Benjamini-Hochberg (plain, for raw p-values)
// ---------------------------------------------------------------------------

/// How the number of tests is inflated for dependence in
/// [`benjamini_hochberg`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Dependence {
    /// Independent (or positively dependent) tests: Benjamini-Hochberg (1995).
    #[default]
    Independent,
    /// Arbitrary dependence: Benjamini-Yekutieli (2001), multiplying by
    /// `1 + 1/2 + ... + 1/n`.
    Arbitrary,
}

/// Benjamini-Hochberg q-values for a set of raw p-values.
///
/// For sorted p-values `p_(1) <= ... <= p_(n)` the q-value of the `i`-th is
/// `min over j >= i of  n p_(j) / j`, capped at 1. A NaN p-value is ignored (it
/// does not count toward `n`) and yields `None` in the output; a p-value outside
/// `[0, 1]` is an error.
///
/// This is the textbook procedure and is deliberately separate from AFNI's
/// stored-curve algorithm ([`fdrize`]), which floors p at 1e-15, ignores
/// p >= 0.9999, estimates the number of true positives, and works in single
/// precision. Results differ from AFNI's by design.
pub fn benjamini_hochberg(pvalues: &[f64], dependence: Dependence) -> Result<Vec<Option<QValue>>> {
    for &p in pvalues {
        if !p.is_nan() && !(0.0..=1.0).contains(&p) {
            return Err(Error::InvalidParameter {
                name: "p-value".into(),
                reason: format!("{p} is outside [0, 1]"),
            });
        }
    }
    let mut order: Vec<usize> = (0..pvalues.len())
        .filter(|&i| !pvalues[i].is_nan())
        .collect();
    order.sort_by(|&a, &b| pvalues[a].total_cmp(&pvalues[b]));
    let n = order.len();
    let inflate = match dependence {
        Dependence::Independent => 1.0,
        Dependence::Arbitrary => (1..=n).map(|i| 1.0 / i as f64).sum(),
    };
    let mut out = vec![None; pvalues.len()];
    let mut running = 1.0_f64;
    for (rank, &index) in order.iter().enumerate().rev() {
        let q = (n as f64 * inflate * pvalues[index] / (rank as f64 + 1.0)).min(1.0);
        running = running.min(q);
        out[index] = Some(QValue(running));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::critical_value;

    /// Small deterministic generator (xorshift) so tests need no dependencies.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> f64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64 // uniform in (0, 1)
        }
    }

    fn t_spec() -> StatSpec {
        StatSpec::new(StatKind::Ttest, &[23.0], 0.0)
    }

    /// t statistics: mostly null (uniform p), some signal (small p), random sign.
    fn t_data(n: usize, signal_every: usize) -> Vec<f32> {
        let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
        (0..n)
            .map(|i| {
                let u = rng.next();
                let p = if i % signal_every == 0 { u.powi(6) } else { u };
                let t = critical_value(&t_spec(), p.max(1e-12), Tail::TwoSided).unwrap();
                (if rng.next() < 0.5 { -t } else { t }) as f32
            })
            .collect()
    }

    #[test]
    fn afni_sort_sorts_and_keeps_pairs_together() {
        let mut rng = Rng(12345);
        for n in [0usize, 1, 2, 9, 10, 11, 50, 1000] {
            // Few distinct keys to force many ties.
            let keys: Vec<f32> = (0..n).map(|_| (rng.next() * 7.0).floor() as f32).collect();
            let payload: Vec<f32> = (0..n).map(|i| i as f32).collect();
            let (mut a, mut ia) = (keys.clone(), payload.clone());
            afni_qsort_float_float(&mut a, &mut ia);
            assert!(a.windows(2).all(|w| w[0] <= w[1]), "n={n}: not sorted");
            // Every payload still sits with its original key.
            for (k, p) in a.iter().zip(&ia) {
                assert_eq!(
                    *k, keys[*p as usize],
                    "n={n}: payload detached from its key"
                );
            }
            let mut seen = ia.clone();
            seen.sort_by(f32::total_cmp);
            assert_eq!(seen, payload, "n={n}: payloads must be a permutation");
        }
    }

    #[test]
    fn deterministic_ties_start_the_curve_at_the_smallest_tied_statistic() {
        let stats = t_data(6000, 6);
        let afni = fdr_curves(&t_spec(), &stats, None, &FdrOptions::default()).unwrap();
        let det = fdr_curves(
            &t_spec(),
            &stats,
            None,
            &FdrOptions {
                deterministic_ties: true,
                ..FdrOptions::default()
            },
        )
        .unwrap();
        // Same everywhere except where tied z values were ordered differently:
        // the deterministic curve starts no later than AFNI's.
        assert!(det.fdr.x0() <= afni.fdr.x0());
        assert_eq!(det.fdr.len(), afni.fdr.len());
    }

    #[test]
    fn q_value_wrapper_validates() {
        assert!(QValue::new(0.0).is_ok() && QValue::new(1.0).is_ok());
        assert!(
            QValue::new(-0.1).is_err()
                && QValue::new(1.1).is_err()
                && QValue::new(f64::NAN).is_err()
        );
    }

    #[test]
    fn z_and_q_are_inverse_and_edges_follow_afni() {
        for q in [1e-12, 1e-4, 0.01, 0.05, 0.3, 0.9] {
            let z = z_for_q(QValue::new(q).unwrap()).unwrap();
            let back = q_for_z(z).unwrap().get();
            assert!((back - q).abs() < 1e-12 * q.max(1e-3), "{q}: {z} -> {back}");
        }
        assert_eq!(z_for_q(QValue::new(1.0).unwrap()).unwrap(), 0.0);
        assert!(z_for_q(QValue::new(0.0).unwrap()).is_err());
        // z <= 0 means q = 1 (fdrval).
        assert_eq!(q_for_z(-2.0).unwrap().get(), 1.0);
        // Known: two-sided q = 0.05 at z = 1.959964.
        assert!(
            (z_for_q(QValue::new(0.05).unwrap()).unwrap() - 1.959_963_984_540_054).abs() < 1e-12
        );
    }

    #[test]
    fn benjamini_hochberg_matches_the_textbook_example() {
        let p = [0.01, 0.04, 0.03, 0.005];
        let q: Vec<f64> = benjamini_hochberg(&p, Dependence::Independent)
            .unwrap()
            .into_iter()
            .map(|q| q.unwrap().get())
            .collect();
        // sorted p: .005 .01 .03 .04 -> n p / rank = .02 .02 .04 .04 (already monotone).
        assert_eq!(q, vec![0.02, 0.04, 0.04, 0.02]);
        // Step-up monotonicity: a larger p never has a smaller q.
        let q2 = benjamini_hochberg(&[0.001, 0.2, 0.2, 0.9], Dependence::Independent).unwrap();
        let v: Vec<f64> = q2.iter().map(|q| q.unwrap().get()).collect();
        assert!(v[0] <= v[1] && v[1] == v[2] && v[2] <= v[3]);
        // Arbitrary dependence multiplies by the harmonic number (capped at 1).
        let by = benjamini_hochberg(&[0.001, 0.002, 0.5], Dependence::Arbitrary).unwrap();
        let h3: f64 = 1.0 + 0.5 + 1.0 / 3.0;
        assert!((by[0].unwrap().get() - (3.0 * h3 * 0.001).min(1.0)).abs() < 1e-15);
        // NaN is skipped and does not count toward n; out-of-range is an error.
        let with_nan =
            benjamini_hochberg(&[0.01, f64::NAN, 0.04], Dependence::Independent).unwrap();
        assert!(with_nan[1].is_none());
        assert!((with_nan[0].unwrap().get() - 0.02).abs() < 1e-15); // n = 2, not 3
        assert!(benjamini_hochberg(&[0.5, 1.5], Dependence::Independent).is_err());
        assert!(benjamini_hochberg(&[], Dependence::Independent)
            .unwrap()
            .is_empty());
    }

    #[test]
    fn fdrize_ignores_what_afni_ignores_and_needs_twenty_samples() {
        let mut stats = t_data(400, 10);
        stats[0] = 0.0; // zero statistic: p = 1, ignored
        stats[1] = f32::NAN; // non-finite: ignored
        let mut mask = vec![true; 400];
        mask[2] = false; // masked out
        let r = fdrize(&t_spec(), &stats, Some(&mask), &FdrOptions::default()).unwrap();
        for ignored in [0, 1, 2] {
            assert_eq!(r.values[ignored], 0.0, "z-score of an ignored sample is 0");
        }
        let q = fdrize(
            &t_spec(),
            &stats,
            Some(&mask),
            &FdrOptions {
                output: FdrOutput::QValue,
                ..FdrOptions::default()
            },
        )
        .unwrap();
        for ignored in [0, 1, 2] {
            assert_eq!(q.values[ignored], 1.0, "q-value of an ignored sample is 1");
        }
        assert!(r.used <= 397);
        assert!(fdrize(&t_spec(), &stats[..15], None, &FdrOptions::default()).is_err());
        assert!(fdrize(&t_spec(), &stats, Some(&mask[..3]), &FdrOptions::default()).is_err());
    }

    #[test]
    fn q_values_decrease_with_significance_and_stay_in_range() {
        let stats = t_data(2000, 8);
        let opts = FdrOptions {
            output: FdrOutput::QValue,
            ..FdrOptions::default()
        };
        let r = fdrize(&t_spec(), &stats, None, &opts).unwrap();
        let mut pairs: Vec<(f32, f32)> = stats
            .iter()
            .map(|t| t.abs())
            .zip(r.values.iter().copied())
            .collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        for w in pairs.windows(2) {
            assert!(
                w[1].1 <= w[0].1 + 1e-6,
                "larger |t| must not have larger q: {w:?}"
            );
        }
        assert!(r.values.iter().all(|q| (0.0..=1.0).contains(q)));
        // Strong signal exists, so some q is small; and the adjustment can be turned off.
        assert!(r.values.iter().any(|&q| q < 0.01));
        let raw = fdrize(
            &t_spec(),
            &stats,
            None,
            &FdrOptions {
                skip_true_positive_adjustment: true,
                ..opts
            },
        )
        .unwrap();
        // The adjustment (scaling by the null fraction) only ever lowers q.
        assert!(raw
            .values
            .iter()
            .zip(&r.values)
            .all(|(raw, adj)| adj <= &(raw + 1e-6)));
        assert!(r.true_positives.is_some());
    }

    #[test]
    fn arbitrary_dependence_raises_q() {
        let stats = t_data(1000, 8);
        let base = FdrOptions {
            output: FdrOutput::QValue,
            ..FdrOptions::default()
        };
        let a = fdrize(&t_spec(), &stats, None, &base).unwrap();
        let b = fdrize(
            &t_spec(),
            &stats,
            None,
            &FdrOptions {
                arbitrary_dependence: true,
                ..base
            },
        )
        .unwrap();
        assert!(a.values.iter().zip(&b.values).all(|(a, b)| b >= a));
        assert!(a.values.iter().zip(&b.values).any(|(a, b)| b > a));
    }

    #[test]
    fn p_value_input_matches_statistic_input() {
        let stats = t_data(600, 6);
        let evaluator = PValueEvaluator::new(&t_spec(), Tail::TwoSided).unwrap();
        let p: Vec<f32> = stats
            .iter()
            .map(|&t| evaluator.evaluate(f64::from(t.abs())).unwrap().p() as f32)
            .collect();
        let opts = FdrOptions::default();
        let from_p = fdrize_p_values(&p, None, &opts).unwrap();
        let from_t = fdrize(&t_spec(), &stats, None, &opts).unwrap();
        assert_eq!(from_p.values, from_t.values);
        // Bad p-values (negative, above 1, NaN) are ignored, not errors.
        let mut bad = p.clone();
        bad[0] = -0.5;
        bad[1] = 1.5;
        bad[2] = f32::NAN;
        let r = fdrize_p_values(&bad, None, &opts).unwrap();
        assert_eq!(&r.values[..3], &[0.0, 0.0, 0.0]);
    }

    #[test]
    fn built_curves_have_the_afni_shape_and_agree_with_per_sample_q() {
        let stats = t_data(6000, 6);
        let built = fdr_curves(&t_spec(), &stats, None, &FdrOptions::default()).unwrap();
        let c = &built.fdr;
        assert_eq!(c.len(), 101);
        // z(q) rises with the statistic (never decreases) and stays within the cap.
        for w in c.samples().windows(2) {
            assert!(w[1] >= w[0] - 1e-6, "z(q) curve must not decrease: {w:?}");
        }
        assert!(c.samples().iter().all(|z| (0.0..=Z_TOP + 1e-6).contains(z)));
        // Looking a sample's threshold up on the curve reproduces (roughly) its own q.
        let q_opts = FdrOptions {
            output: FdrOutput::QValue,
            ..FdrOptions::default()
        };
        let per_sample = fdrize(&t_spec(), &stats, None, &q_opts).unwrap();
        let mut checked = 0;
        for (t, q) in stats.iter().zip(&per_sample.values) {
            if (0.001..0.2).contains(q) {
                let from_curve = q_value_for_threshold(c, f64::from(t.abs())).unwrap().get();
                assert!(
                    (from_curve - f64::from(*q)).abs() < 0.25 * f64::from(*q) + 1e-3,
                    "|t|={} per-sample q={q} curve q={from_curve}",
                    t.abs()
                );
                checked += 1;
            }
        }
        assert!(checked > 50, "only {checked} samples compared");
        // And the inverse: the threshold for q reproduces q (approximately).
        for q in [0.2, 0.05, 0.01] {
            let thr = threshold_for_q(c, QValue::new(q).unwrap(), None).unwrap();
            let back = q_value_for_threshold(c, thr).unwrap().get();
            assert!((back - q).abs() < 0.1 * q, "q={q}: thr={thr} -> {back}");
        }
        assert!(minimum_q(c).unwrap().get() < 0.001);
    }

    #[test]
    fn threshold_for_q_edges() {
        let stats = t_data(4000, 6);
        let c = fdr_curves(&t_spec(), &stats, None, &FdrOptions::default())
            .unwrap()
            .fdr;
        // q = 1 -> z = 0 -> below the curve -> threshold 0.
        assert_eq!(
            threshold_for_q(&c, QValue::new(1.0).unwrap(), None).unwrap(),
            0.0
        );
        // q beyond what the curve can express: past the last sample, or just above max |stat|.
        let tiny = QValue::new(1e-15).unwrap();
        let end = c.x0() + c.dx() * c.len() as f64;
        assert_eq!(threshold_for_q(&c, tiny, None).unwrap(), end);
        assert_eq!(
            threshold_for_q(&c, tiny, Some(end * 3.0)).unwrap(),
            end * 3.0 * 1.000002
        );
        assert_eq!(threshold_for_q(&c, tiny, Some(end * 0.5)).unwrap(), end);
        assert!(threshold_for_q(&c, QValue::new(0.0).unwrap(), None).is_err());
    }

    #[test]
    fn missed_detection_curve_and_lookup() {
        let stats = t_data(8000, 4); // many true positives
        let built = fdr_curves(&t_spec(), &stats, None, &FdrOptions::default()).unwrap();
        let mdf = built
            .mdf
            .expect("an MDF curve is built when there are many true positives");
        assert!(built.true_positives.unwrap() > 8.0);
        // Axis is log10(p); values are fractions that fall toward 0 as p grows.
        assert!(mdf.x0() < 0.0 && mdf.dx() > 0.0);
        assert!(mdf.samples().iter().all(|m| (0.0..=1.0).contains(m)));
        assert_eq!(*mdf.samples().last().unwrap(), 0.0);
        assert_eq!(missed_detection_fraction(&mdf, 0.0).unwrap(), 0.999);
        assert_eq!(missed_detection_fraction(&mdf, 1.0).unwrap(), 0.0);
        assert!(missed_detection_fraction(&mdf, f64::NAN).is_err());
        let strict = missed_detection_fraction(&mdf, 1e-8).unwrap();
        let loose = missed_detection_fraction(&mdf, 0.05).unwrap();
        assert!(
            strict >= loose,
            "a stricter p misses more: {strict} vs {loose}"
        );
    }

    #[test]
    fn classic_statistics_use_afni_tails_and_others_need_one() {
        let f = StatSpec::new(StatKind::Ftest, &[3.0, 40.0], 0.0);
        let stats: Vec<f32> = (1..400).map(|i| i as f32 * 0.02).collect();
        assert!(fdrize(&f, &stats, None, &FdrOptions::default()).is_ok());
        assert_eq!(afni_tail(&f), Some(Tail::Upper));
        assert_eq!(afni_tail(&t_spec()), Some(Tail::TwoSided));
        let normal = StatSpec::new(StatKind::Normal, &[0.0, 1.0], 0.0);
        assert!(fdrize(&normal, &stats, None, &FdrOptions::default()).is_err());
        let explicit = FdrOptions {
            tail: Some(Tail::Upper),
            ..FdrOptions::default()
        };
        assert!(fdrize(&normal, &stats, None, &explicit).is_ok());
        // Direct-probability kinds cannot be treated as statistics.
        let direct = StatSpec::new(StatKind::Pval, &[], 0.0);
        assert!(fdrize(&direct, &stats, None, &explicit).is_err());
    }

    #[test]
    fn degenerate_inputs_do_not_panic() {
        // All identical statistics: no curve, but also no crash or NaN.
        let same = vec![2.5_f32; 500];
        let r = fdr_curves(&t_spec(), &same, None, &FdrOptions::default());
        assert!(r.is_err() || r.unwrap().fdr.samples().iter().all(|z| z.is_finite()));
        // Nothing significant at all.
        let nulls: Vec<f32> = (1..=300).map(|i| 0.0005 * i as f32).collect();
        let _ = fdr_curves(&t_spec(), &nulls, None, &FdrOptions::default());
    }
}
