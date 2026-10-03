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
// Thresholds: which samples "pass" a display threshold, how failing samples may
// be faded instead of hidden (AFNI's "transparent thresholding"), and how a
// threshold chosen for one statistic is carried to another by matching p-values.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `overlay.rs` applies a `Threshold` and a `FadeModel` to a column of values
//   when it builds a colored overlay.
// * `stats.rs` supplies `p_value` / `critical_value` for `transfer_threshold`.
// * Ported from sumaru's `Threshold`, checked against SUMA's
//   `SUMA_ScaleToMap_Interactive` and AFNI's `AFNI_newnewfunc_overlay` (see the
//   per-item notes below for what each reference does).
//
// THE BOUNDARY IS PART OF THE DEFINITION
//
// SUMA's threshold modes disagree about whether the boundary value passes, and
// users notice a voxel exactly at the threshold. Here every mode states it:
//
//   Above(t)            passes  value >= t          (boundary passes)
//   Below(t)            passes  value <= t          (boundary passes)
//   Between{lo, hi}     passes  lo <= value <= hi   (both ends pass)
//   Outside{lo, hi}     passes  value < lo || value > hi   (both ends FAIL)
//   AbsoluteAbove(t)    passes  |value| >= t        (boundary passes)
//
// `Outside` is SUMA's "hide values in [lo, hi]" (so the ends are hidden) while
// `AbsoluteAbove` is SUMA's "hide values with -t < v < t" (so the ends pass). They
// look like the same thing for lo = -t, hi = t but differ at exactly +-t; both are
// kept, honestly named, rather than one pretending to be the other.
//
// A non-finite threshold value never passes (unless the threshold is Off); see
// `MissingThreshold` for AFNI/SUMA's different behavior.
// ---------------------------------------------------------------------------

//! Display thresholds, transparent thresholding, and matched-p-value transfer.

use crate::error::{Error, Result};
use crate::numeric::ensure_finite;
use crate::stat::StatSpec;
use crate::stats::{critical_value_ln, p_value, Tail};

/// A display threshold on a threshold column. See the module docs for the exact
/// boundary behavior of each mode.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum Threshold {
    /// No thresholding: every sample passes.
    #[default]
    Off,
    /// `value >= t`.
    Above(f64),
    /// `value <= t`. (Not a SUMA mode; kept from sumaru.)
    Below(f64),
    /// `lo <= value <= hi`.
    Between {
        /// Lower end (passes).
        lo: f64,
        /// Upper end (passes).
        hi: f64,
    },
    /// `value < lo || value > hi`: SUMA's "hide values inside [lo, hi]".
    Outside {
        /// Lower end (fails).
        lo: f64,
        /// Upper end (fails).
        hi: f64,
    },
    /// `|value| >= t`: symmetric thresholding stated directly, so callers never
    /// construct a `[-t, t]` pair by hand. `t` must be non-negative.
    AbsoluteAbove(f64),
}

impl Threshold {
    /// Check the numbers: all finite, `lo <= hi`, and `AbsoluteAbove` non-negative.
    pub fn validate(&self) -> Result<()> {
        match *self {
            Threshold::Off => Ok(()),
            Threshold::Above(t) | Threshold::Below(t) => ensure_finite("threshold", t).map(|_| ()),
            Threshold::AbsoluteAbove(t) => {
                ensure_finite("threshold", t)?;
                if t < 0.0 {
                    return Err(Error::InvalidParameter {
                        name: "threshold".into(),
                        reason: format!("an absolute threshold must be non-negative, got {t}"),
                    });
                }
                Ok(())
            }
            Threshold::Between { lo, hi } | Threshold::Outside { lo, hi } => {
                ensure_finite("threshold low end", lo)?;
                ensure_finite("threshold high end", hi)?;
                if lo > hi {
                    return Err(Error::InvalidParameter {
                        name: "threshold range".into(),
                        reason: format!("low end {lo} is above high end {hi}"),
                    });
                }
                Ok(())
            }
        }
    }

    /// Whether `value` passes. NaN and infinities never pass (except `Off`).
    pub fn passes(&self, value: f64) -> bool {
        if matches!(self, Threshold::Off) {
            return true;
        }
        if !value.is_finite() {
            return false;
        }
        match *self {
            Threshold::Off => true,
            Threshold::Above(t) => value >= t,
            Threshold::Below(t) => value <= t,
            Threshold::Between { lo, hi } => value >= lo && value <= hi,
            Threshold::Outside { lo, hi } => value < lo || value > hi,
            Threshold::AbsoluteAbove(t) => value.abs() >= t,
        }
    }

    /// The threshold magnitude used by the AFNI and SUMA fade models, which are
    /// defined only for "above" style thresholds.
    fn fade_magnitude(&self) -> Option<f64> {
        match *self {
            Threshold::Above(t) | Threshold::AbsoluteAbove(t) => Some(t),
            _ => None,
        }
    }
}

/// What to do with a sample whose THRESHOLD value is missing (NaN or infinite).
///
/// AFNI and SUMA compare with `<`/`>`, which are false for NaN, so a NaN threshold
/// value is never hidden and the sample is drawn. That is almost certainly not
/// what anyone wants, so the default here is to hide it; `Show` reproduces the
/// reference programs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MissingThreshold {
    /// A missing threshold value fails the threshold (the default).
    #[default]
    Hide,
    /// A missing threshold value passes, as in AFNI and SUMA.
    Show,
}

/// The shape of the opacity ramp for samples below threshold.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum FadeCurve {
    /// Opacity proportional to `|value| / threshold` (AFNI `alcode` 1).
    Linear,
    /// The square of the linear ramp (AFNI `alcode` 2; the default).
    #[default]
    Quadratic,
    /// Cube. Not an AFNI option; a sumaru enhancement.
    Cubic,
    /// Fourth power. Not an AFNI option; a sumaru enhancement.
    Quartic,
}

impl FadeCurve {
    fn apply(self, ratio: f64) -> f64 {
        match self {
            FadeCurve::Linear => ratio,
            FadeCurve::Quadratic => ratio * ratio,
            FadeCurve::Cubic => ratio * ratio * ratio,
            FadeCurve::Quartic => (ratio * ratio) * (ratio * ratio),
        }
    }
}

/// AFNI's ceiling on the opacity of a sub-threshold voxel, 222/255, "to make sure
/// there is some distinction between above-threshold voxels and below threshold
/// rebel scum" (`ALFABYTE` in `afni_func.c`).
pub const AFNI_SUBTHRESHOLD_MAX_ALPHA: f32 = 222.0 / 255.0;

/// How samples that FAIL the threshold are made less prominent, instead of being
/// hidden outright. Passing samples always have factor 1.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FadeModel {
    /// AFNI's volumetric transparent thresholding (`AFNI_newnewfunc_overlay`).
    ///
    /// For a failing sample with threshold value `v` and threshold `t`
    /// (`Above(t)`, or `AbsoluteAbove(t)`): `x = |v| (1 - floor) / t`, opacity
    /// `255 x + 255 floor` (linear) or `255 x^2 + 255 floor` (quadratic), rounded to
    /// the nearest byte (ties to even) and clamped to 0..=222, then divided by 255.
    /// A threshold value of exactly 0 is transparent. `floor` is AFNI's
    /// `thr_alpha_floor` (0 in practice).
    ///
    /// For `Above(t)` a negative value cannot pass, so it is transparent here;
    /// whether AFNI's one-sided GUI path draws it is not pinned down (see roadmap).
    Afni {
        /// Linear or quadratic ramp. Cubic and quartic are not AFNI curves and are
        /// rejected by `validate`.
        curve: FadeCurve,
        /// The opacity floor in `[0, 1)`.
        floor: f32,
    },
    /// SUMA's surface "alpha opacity falloff" (`alphaOpacitiesForOverlay`):
    /// opacity `min(1, |v| / t)`, squared if quadratic; 1 if `t == 0`. No byte
    /// quantization and no 222 ceiling.
    Suma {
        /// Linear or quadratic ramp.
        curve: FadeCurve,
    },
    /// Fade with distance from the nearest edge of the passing region, over
    /// `width` data units (or the magnitude of that edge). Works for every mode,
    /// including `Between`, `Outside` and `Below`. This is sumaru's generalization,
    /// not an AFNI behavior. For `Above(t)`/`AbsoluteAbove(t)` with
    /// `FadeWidth::BoundaryMagnitude` it equals the linear/quadratic ramps above
    /// (without the byte quantization or ceiling).
    Boundary {
        /// Ramp exponent.
        curve: FadeCurve,
        /// Distance over which opacity falls from 1 to 0.
        width: FadeWidth,
    },
}

/// The distance over which a [`FadeModel::Boundary`] fade falls to zero.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FadeWidth {
    /// The magnitude of the nearest threshold edge (equivalent to AFNI's
    /// fade-to-zero for `Above(t)` and symmetric thresholds).
    BoundaryMagnitude,
    /// An explicit distance, in the units of the threshold column.
    Absolute(f64),
}

/// Round half to even, like C's `rintf` in the default rounding mode.
fn rint(x: f32) -> f32 {
    let r = x.round(); // rounds halves away from zero
    if (x - x.trunc()).abs() == 0.5 && r % 2.0 != 0.0 {
        r - x.signum() // a tie that landed on an odd integer: step back to even
    } else {
        r
    }
}

impl FadeModel {
    /// Reject combinations that are not defined (for example an AFNI fade with
    /// `Between`, or a cubic curve under an AFNI model).
    pub fn validate(&self, threshold: &Threshold) -> Result<()> {
        match *self {
            FadeModel::Afni { curve, floor } => {
                if !matches!(curve, FadeCurve::Linear | FadeCurve::Quadratic) {
                    return Err(Error::InvalidParameter {
                        name: "fade curve".into(),
                        reason: "AFNI fades are linear or quadratic only".into(),
                    });
                }
                ensure_finite("alpha floor", f64::from(floor))?;
                if !(0.0..1.0).contains(&floor) {
                    return Err(Error::InvalidParameter {
                        name: "alpha floor".into(),
                        reason: format!("{floor} is outside [0, 1)"),
                    });
                }
                Self::require_above(threshold, "AFNI")
            }
            FadeModel::Suma { curve } => {
                if !matches!(curve, FadeCurve::Linear | FadeCurve::Quadratic) {
                    return Err(Error::InvalidParameter {
                        name: "fade curve".into(),
                        reason: "SUMA fades are linear or quadratic only".into(),
                    });
                }
                Self::require_above(threshold, "SUMA")
            }
            FadeModel::Boundary { width, .. } => {
                if let FadeWidth::Absolute(w) = width {
                    ensure_finite("fade width", w)?;
                    if w <= 0.0 {
                        return Err(Error::InvalidParameter {
                            name: "fade width".into(),
                            reason: format!("must be positive, got {w}"),
                        });
                    }
                }
                Ok(())
            }
        }
    }

    fn require_above(threshold: &Threshold, who: &str) -> Result<()> {
        if matches!(threshold, Threshold::Off) || threshold.fade_magnitude().is_some() {
            Ok(())
        } else {
            Err(Error::Unsupported(format!(
                "the {who} fade is defined only for Above and AbsoluteAbove thresholds; \
                 use FadeModel::Boundary for {threshold:?}"
            )))
        }
    }

    /// The opacity multiplier for a sample whose threshold value is `value`: 1 if
    /// it passes, otherwise the fade ramp in `[0, 1]`. A non-finite value is 0.
    /// The model must have been validated against `threshold`.
    pub fn factor(&self, threshold: &Threshold, value: f64) -> f32 {
        if matches!(threshold, Threshold::Off) {
            return 1.0;
        }
        if !value.is_finite() {
            return 0.0;
        }
        if threshold.passes(value) {
            return 1.0;
        }
        match *self {
            FadeModel::Afni { curve, floor } => {
                let Some(t) = threshold.fade_magnitude() else {
                    return 0.0;
                };
                // AFNI: a threshold value of exactly 0 is rejected outright, and for
                // a one-sided threshold a negative value never fades in.
                let v = match threshold {
                    Threshold::Above(_) => value,
                    _ => value.abs(),
                };
                if v <= 0.0 || t <= 0.0 {
                    return 0.0;
                }
                let scale = (1.0 - floor) / t as f32; // `ft` in afni_func.c
                let x = v as f32 * scale;
                let ramp = match curve {
                    FadeCurve::Linear => x,
                    _ => x * x,
                };
                // 255 * ramp + 255 * floor, as a byte capped at 222.
                let byte = rint(255.0 * ramp + 255.0 * floor).clamp(0.0, 222.0);
                byte / 255.0
            }
            FadeModel::Suma { curve } => {
                let Some(t) = threshold.fade_magnitude() else {
                    return 0.0;
                };
                let denom = t.max(0.0);
                let o = if denom == 0.0 {
                    1.0
                } else {
                    (value.abs() / denom).min(1.0)
                };
                curve.apply(o) as f32
            }
            FadeModel::Boundary { curve, width } => {
                let (distance, edge) = match *threshold {
                    Threshold::Off => return 1.0,
                    Threshold::Above(t) => (t - value, t),
                    Threshold::Below(t) => (value - t, t),
                    Threshold::AbsoluteAbove(t) => (t - value.abs(), t),
                    Threshold::Between { lo, hi } => {
                        if value < lo {
                            (lo - value, lo)
                        } else {
                            (value - hi, hi)
                        }
                    }
                    Threshold::Outside { lo, hi } => {
                        let (to_lo, to_hi) = (value - lo, hi - value);
                        if to_lo <= to_hi {
                            (to_lo, lo)
                        } else {
                            (to_hi, hi)
                        }
                    }
                };
                let span = match width {
                    FadeWidth::BoundaryMagnitude => edge.abs(),
                    FadeWidth::Absolute(w) => w,
                };
                // `span` is NaN-checked by the comparison form below on purpose.
                if !distance.is_finite() || distance < 0.0 || span.is_nan() || span <= 0.0 {
                    return 0.0;
                }
                curve.apply((1.0 - distance / span).clamp(0.0, 1.0)) as f32
            }
        }
    }
}

/// Carry a threshold from one statistic to another by matching its p-value: the
/// probability of `threshold` under `source` (for the chosen `tail`) is looked up
/// in `destination` for the same `tail`.
///
/// The tail is preserved on purpose: a two-sided t threshold of 2.0 means "p = 0.07
/// two-sided", and the equivalent z threshold is the one with the same two-sided
/// probability, not the one-sided probability. Fails if either statistic does not
/// support the tail, if `threshold` is outside the source's domain, or if the
/// probability has no finite counterpart in `destination`.
///
/// ```
/// use afni_core::stat::{StatKind, StatSpec};
/// use afni_core::stats::Tail;
/// use afni_core::threshold::transfer_threshold;
///
/// let t = StatSpec::new(StatKind::Ttest, &[1000.0], 0.0);
/// let z = StatSpec::new(StatKind::Zscore, &[], 0.0);
/// // With 1000 degrees of freedom, t is nearly normal: the thresholds nearly agree.
/// let zt = transfer_threshold(&t, &z, Tail::TwoSided, 1.96)?;
/// assert!((zt - 1.96).abs() < 0.01);
/// # Ok::<(), afni_core::Error>(())
/// ```
pub fn transfer_threshold(
    source: &StatSpec,
    destination: &StatSpec,
    tail: Tail,
    threshold: f64,
) -> Result<f64> {
    let probability = p_value(source, threshold, tail)?;
    if probability.ln_p() == 0.0 {
        // p = 1: the threshold sits at the centre of the source distribution, where
        // every value is "at least as central", so the destination has none either.
        return Err(Error::NoSolution(
            "a probability of 1 has no finite threshold in the destination".into(),
        ));
    }
    // The log form keeps probabilities far below 1e-308 invertible.
    critical_value_ln(destination, probability.ln_p(), tail)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stat::StatKind;

    #[test]
    fn boundary_behavior_is_exact_for_every_mode() {
        let t = 2.0;
        assert!(Threshold::Above(t).passes(2.0) && !Threshold::Above(t).passes(1.999));
        assert!(Threshold::Below(t).passes(2.0) && !Threshold::Below(t).passes(2.001));
        let between = Threshold::Between { lo: -1.0, hi: 1.0 };
        assert!(between.passes(-1.0) && between.passes(1.0) && !between.passes(1.0001));
        let outside = Threshold::Outside { lo: -1.0, hi: 1.0 };
        assert!(
            !outside.passes(-1.0) && !outside.passes(1.0),
            "SUMA hides the ends"
        );
        assert!(outside.passes(-1.0001) && outside.passes(1.0001) && !outside.passes(0.0));
        let abs = Threshold::AbsoluteAbove(1.0);
        assert!(
            abs.passes(-1.0) && abs.passes(1.0),
            "|v| >= t passes the ends"
        );
        assert!(!abs.passes(0.999) && !abs.passes(-0.999));
        // Same-looking thresholds differ only at the boundary, by design.
        for v in [-1.0, 1.0] {
            assert_ne!(outside.passes(v), abs.passes(v));
        }
        for v in [-3.0, -1.5, 0.0, 0.5, 2.0] {
            assert_eq!(outside.passes(v), abs.passes(v), "{v}");
        }
        assert!(Threshold::Off.passes(f64::NAN) && Threshold::Off.passes(-5.0));
    }

    #[test]
    fn non_finite_values_never_pass_a_real_threshold() {
        for th in [
            Threshold::Above(0.0),
            Threshold::Below(0.0),
            Threshold::Between { lo: -1.0, hi: 1.0 },
            Threshold::Outside { lo: -1.0, hi: 1.0 },
            Threshold::AbsoluteAbove(0.0),
        ] {
            for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                assert!(!th.passes(v), "{th:?} {v}");
            }
        }
    }

    #[test]
    fn validation() {
        assert!(Threshold::Above(f64::NAN).validate().is_err());
        assert!(Threshold::AbsoluteAbove(-1.0).validate().is_err());
        assert!(Threshold::Between { lo: 2.0, hi: 1.0 }.validate().is_err());
        assert!(Threshold::Outside {
            lo: 0.0,
            hi: f64::INFINITY
        }
        .validate()
        .is_err());
        assert!(Threshold::Between { lo: 1.0, hi: 1.0 }.validate().is_ok());
        assert!(Threshold::Off.validate().is_ok());
    }

    #[test]
    fn afni_fade_follows_the_c_formula_with_byte_rounding_and_the_222_cap() {
        let m = FadeModel::Afni {
            curve: FadeCurve::Linear,
            floor: 0.0,
        };
        let th = Threshold::AbsoluteAbove(4.0);
        m.validate(&th).unwrap();
        // v = 2 -> x = 0.5 -> 127.5 -> rint (ties to even) = 128 -> 128/255.
        assert_eq!(m.factor(&th, 2.0), 128.0 / 255.0);
        assert_eq!(m.factor(&th, -2.0), 128.0 / 255.0, "negatives use |v|");
        // v = 1 -> 63.75 -> 64.
        assert_eq!(m.factor(&th, 1.0), 64.0 / 255.0);
        // Just below threshold: ramp would be ~255, the cap holds it at 222.
        assert_eq!(m.factor(&th, 3.99), 222.0 / 255.0);
        // At and above the threshold: passes, factor exactly 1 (opaque, uncapped).
        assert_eq!(m.factor(&th, 4.0), 1.0);
        assert_eq!(m.factor(&th, -9.0), 1.0);
        // Exactly zero is rejected outright.
        assert_eq!(m.factor(&th, 0.0), 0.0);
        // Quadratic: 255 x^2.
        let q = FadeModel::Afni {
            curve: FadeCurve::Quadratic,
            floor: 0.0,
        };
        assert_eq!(q.factor(&th, 2.0), 64.0 / 255.0, "0.25 * 255 = 63.75 -> 64");
        // Opacity floor: 255 * ((1 - f) v / t)^2 + 255 f.
        let fl = FadeModel::Afni {
            curve: FadeCurve::Quadratic,
            floor: 0.2,
        };
        let x = 2.0_f32 * (1.0 - 0.2) / 4.0;
        assert_eq!(
            fl.factor(&th, 2.0),
            rint(255.0 * x * x + 255.0 * 0.2) / 255.0
        );
    }

    #[test]
    fn rint_rounds_ties_to_even() {
        assert_eq!(rint(0.5), 0.0);
        assert_eq!(rint(1.5), 2.0);
        assert_eq!(rint(2.5), 2.0);
        assert_eq!(rint(3.5), 4.0);
        assert_eq!(rint(127.5), 128.0);
        assert_eq!(rint(63.75), 64.0);
        assert_eq!(rint(-1.5), -2.0);
        assert_eq!(rint(2.4), 2.0);
    }

    #[test]
    fn suma_fade_has_no_cap_or_quantization() {
        let th = Threshold::AbsoluteAbove(4.0);
        let lin = FadeModel::Suma {
            curve: FadeCurve::Linear,
        };
        let quad = FadeModel::Suma {
            curve: FadeCurve::Quadratic,
        };
        lin.validate(&th).unwrap();
        assert_eq!(lin.factor(&th, 3.99), (3.99_f64 / 4.0) as f32);
        assert_eq!(quad.factor(&th, 2.0), 0.25);
        assert_eq!(lin.factor(&th, 5.0), 1.0);
        // A zero threshold means "no fade".
        let zero = Threshold::Above(0.0);
        assert_eq!(lin.factor(&zero, -1.0), 1.0);
    }

    #[test]
    fn one_sided_afni_fade_treats_negative_values_as_transparent() {
        let m = FadeModel::Afni {
            curve: FadeCurve::Linear,
            floor: 0.0,
        };
        let th = Threshold::Above(4.0);
        assert!(m.factor(&th, 2.0) > 0.0);
        assert_eq!(m.factor(&th, -2.0), 0.0);
    }

    #[test]
    fn boundary_fade_covers_every_mode_and_matches_the_simple_ramps() {
        let b = FadeModel::Boundary {
            curve: FadeCurve::Linear,
            width: FadeWidth::BoundaryMagnitude,
        };
        let th = Threshold::AbsoluteAbove(4.0);
        // Same ramp as SUMA's linear fade (no cap, no quantization).
        for v in [0.5, 1.0, 2.0, 3.5, -3.0] {
            let suma = FadeModel::Suma {
                curve: FadeCurve::Linear,
            }
            .factor(&th, v);
            assert!((b.factor(&th, v) - suma).abs() < 1e-6, "{v}");
        }
        let outside = Threshold::Outside { lo: -2.0, hi: 2.0 };
        assert_eq!(b.factor(&outside, 0.0), 0.0, "deep inside the hidden band");
        assert!(b.factor(&outside, 1.9) > b.factor(&outside, 1.0));
        let between = Threshold::Between { lo: 1.0, hi: 3.0 };
        assert!(b.factor(&between, 0.9) > 0.0 && b.factor(&between, 2.0) == 1.0);
        let wide = FadeModel::Boundary {
            curve: FadeCurve::Linear,
            width: FadeWidth::Absolute(10.0),
        };
        assert!((wide.factor(&Threshold::Above(5.0), 0.0) - 0.5).abs() < 1e-6);
        assert_eq!(b.factor(&th, f64::NAN), 0.0);
    }

    #[test]
    fn invalid_fade_combinations_are_rejected() {
        let afni = FadeModel::Afni {
            curve: FadeCurve::Quadratic,
            floor: 0.0,
        };
        assert!(afni
            .validate(&Threshold::Between { lo: 0.0, hi: 1.0 })
            .is_err());
        assert!(afni
            .validate(&Threshold::Outside { lo: 0.0, hi: 1.0 })
            .is_err());
        assert!(afni.validate(&Threshold::Below(1.0)).is_err());
        assert!(afni.validate(&Threshold::AbsoluteAbove(1.0)).is_ok());
        assert!(FadeModel::Afni {
            curve: FadeCurve::Cubic,
            floor: 0.0
        }
        .validate(&Threshold::Above(1.0))
        .is_err());
        assert!(FadeModel::Afni {
            curve: FadeCurve::Linear,
            floor: 1.0
        }
        .validate(&Threshold::Above(1.0))
        .is_err());
        assert!(FadeModel::Suma {
            curve: FadeCurve::Quartic
        }
        .validate(&Threshold::Above(1.0))
        .is_err());
        assert!(FadeModel::Boundary {
            curve: FadeCurve::Cubic,
            width: FadeWidth::Absolute(0.0)
        }
        .validate(&Threshold::Off)
        .is_err());
    }

    #[test]
    fn off_threshold_never_fades() {
        for m in [
            FadeModel::Afni {
                curve: FadeCurve::Linear,
                floor: 0.0,
            },
            FadeModel::Suma {
                curve: FadeCurve::Quadratic,
            },
            FadeModel::Boundary {
                curve: FadeCurve::Linear,
                width: FadeWidth::BoundaryMagnitude,
            },
        ] {
            assert_eq!(m.factor(&Threshold::Off, 0.0), 1.0);
            assert_eq!(m.factor(&Threshold::Off, f64::NAN), 1.0);
        }
    }

    #[test]
    fn matched_p_value_transfer_preserves_the_tail() {
        let t = StatSpec::new(StatKind::Ttest, &[10.0], 0.0);
        let z = StatSpec::new(StatKind::Zscore, &[], 0.0);
        // Two-sided p of t = 2.2281 with 10 dof is 0.05, whose z threshold is 1.96.
        let two = transfer_threshold(&t, &z, Tail::TwoSided, 2.228_138_851_986_274).unwrap();
        assert!((two - 1.959_963_984_540_054).abs() < 1e-9, "{two}");
        // One-sided p = 0.025 gives the same z; keeping the tail is what matters.
        let one = transfer_threshold(&t, &z, Tail::Upper, 2.228_138_851_986_274).unwrap();
        assert!((one - 1.959_963_984_540_054).abs() < 1e-9, "{one}");
        // Round trip.
        let back = transfer_threshold(&z, &t, Tail::TwoSided, two).unwrap();
        assert!((back - 2.228_138_851_986_274).abs() < 1e-8);
        // F has no two-sided tail, so the transfer refuses rather than guessing.
        let f = StatSpec::new(StatKind::Ftest, &[3.0, 20.0], 0.0);
        assert!(transfer_threshold(&t, &f, Tail::TwoSided, 2.0).is_err());
        // A threshold at the centre has p = 1 and no finite counterpart.
        assert!(transfer_threshold(&t, &z, Tail::TwoSided, 0.0).is_err());
        // Very small probabilities survive in log space.
        let far = transfer_threshold(&z, &t, Tail::Upper, 12.0).unwrap();
        assert!(far > 12.0);
    }
}
