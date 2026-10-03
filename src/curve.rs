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
// `ThresholdCurve`: a validated, evenly spaced table of numbers that AFNI
// stores next to a statistic to describe false-discovery-rate (FDR) or
// missed-detection (MDF) behaviour. Sample `i` sits at `x0 + i * dx`.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `column.rs` attaches optional FDR and MDF curves to a data column, so
//   consumers get typed curves instead of re-parsing attribute strings.
// * Phase 3 (FDR/MDF) will add the lookups and AFNI's exact interpolation on
//   top of this type. Here we only guarantee the table is well-formed.
// * `afni-io` has a raw `ThresholdCurve` that mirrors the file layout without
//   validation (so odd files still round-trip). Its adapter converts into this
//   validated type and reports an error for a malformed curve.
//
// The meaning of the x axis differs by curve (statistic value for FDR,
// log10(p) for MDF); that is Phase 3's job to document and enforce. Here `x` is
// just "the axis".
// ---------------------------------------------------------------------------

//! Validated evenly spaced threshold curves.

use crate::error::{Error, Result};
use crate::numeric::ensure_finite;

/// An evenly spaced curve: `values[i]` is the value at `x0 + i * dx`.
#[derive(Debug, Clone, PartialEq)]
pub struct ThresholdCurve {
    x0: f64,
    dx: f64,
    samples: Vec<f64>,
}

impl ThresholdCurve {
    /// Build a curve. Requires finite `x0`, finite non-zero `dx`, at least two
    /// samples (AFNI cannot interpolate fewer), and all samples finite.
    pub fn new(x0: f64, dx: f64, samples: Vec<f64>) -> Result<Self> {
        ensure_finite("curve x0", x0)?;
        ensure_finite("curve dx", dx)?;
        if dx == 0.0 {
            return Err(Error::InvalidParameter {
                name: "dx".into(),
                reason: "must be non-zero".into(),
            });
        }
        if samples.len() < 2 {
            return Err(Error::InvalidParameter {
                name: "samples".into(),
                reason: format!("need at least 2 samples, got {}", samples.len()),
            });
        }
        for &s in &samples {
            ensure_finite("curve sample", s)?;
        }
        Ok(Self { x0, dx, samples })
    }

    /// Axis position of the first sample.
    pub fn x0(&self) -> f64 {
        self.x0
    }

    /// Axis spacing between samples (never zero; may be negative).
    pub fn dx(&self) -> f64 {
        self.dx
    }

    /// The sample values.
    pub fn samples(&self) -> &[f64] {
        &self.samples
    }

    /// Number of samples (at least 2).
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Always false: a valid curve has at least two samples. Provided because
    /// Rust convention pairs `len` with `is_empty`.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Axis position of sample `i`, or `None` if `i` is out of range.
    pub fn x_at(&self, i: usize) -> Option<f64> {
        (i < self.samples.len()).then_some(self.x0 + self.dx * i as f64)
    }
}

// ---------------------------------------------------------------------------
// AFNI's interpolation (mri_floatvec.c: interp_floatvec, interp_inverse_floatvec)
// ---------------------------------------------------------------------------

/// The four cubic weights of AFNI's `interp_floatvec` at fractional offset `x`
/// in `[0, 1)`, for the points at offsets -1, 0, +1, +2. Note AFNI writes the
/// `1/6` factors as the decimal `0.1666667`, not `1/6`; matching it exactly
/// requires the same constant (the difference is ~2e-8, but it is the
/// definition of "AFNI's" curve).
fn cubic_weights(x: f64) -> [f64; 4] {
    const SIXTH: f64 = 0.1666667;
    [
        x * (1.0 - x) * (x - 2.0) * SIXTH,
        (x + 1.0) * (x - 1.0) * (x - 2.0) * 0.5,
        x * (x + 1.0) * (2.0 - x) * 0.5,
        x * (x + 1.0) * (x - 1.0) * SIXTH,
    ]
}

impl ThresholdCurve {
    /// The curve's value at `x`: AFNI's clamped four-point cubic interpolation,
    /// ported exactly from `interp_floatvec` (`mri_floatvec.c`).
    ///
    /// * Outside the sampled range the nearest end sample is returned.
    /// * Between samples `i` and `i + 1` the cubic through samples
    ///   `i - 1 ..= i + 2` is used (end samples repeated at the edges) and the
    ///   result is then clamped to lie between samples `i` and `i + 1`, so it
    ///   can never overshoot a monotone curve. This is *not* a Catmull-Rom
    ///   spline (sumaru's former choice); the weights differ.
    /// * AFNI quirk, kept for parity: with only **two** samples (`itop <= 1`)
    ///   AFNI returns the first sample for every `x`.
    /// * A NaN `x` gives NaN.
    pub fn interpolate(&self, x: f64) -> f64 {
        let ar = &self.samples;
        let itop = ar.len() - 1;
        if itop <= 1 {
            return ar[0];
        }
        if x.is_nan() {
            return f64::NAN;
        }
        let fx = (x - self.x0) / self.dx;
        if fx <= 0.0 {
            return ar[0];
        }
        if fx >= itop as f64 {
            return ar[itop];
        }
        // x lies between sample `ix` and `ix + 1`, at fractional offset `fx`.
        let ix = fx as usize;
        let fx = fx - ix as f64;
        let im1 = ix.saturating_sub(1);
        let ip1 = ix + 1;
        let ip2 = (ip1 + 1).min(itop);
        let w = cubic_weights(fx);
        let value = w[0] * ar[im1] + w[1] * ar[ix] + w[2] * ar[ip1] + w[3] * ar[ip2];
        // Keep the result inside the local range of values.
        let (lo, hi) = if ar[ix] > ar[ip1] {
            (ar[ip1], ar[ix])
        } else {
            (ar[ix], ar[ip1])
        };
        value.clamp(lo, hi)
    }

    /// The `x` at which a (roughly monotone) curve reaches `y`: AFNI's
    /// `interp_inverse_floatvec`, ported exactly.
    ///
    /// * `y` off either end gives the corresponding end `x` (the first or the
    ///   last sample position).
    /// * Otherwise the first bracketing pair of samples is found, linear
    ///   interpolation gives a first estimate, and two regula-falsi steps
    ///   (about `0.05 * dx` either side) refine it against the *cubic* forward
    ///   interpolation; the candidate with the smallest residual wins. Because
    ///   forward and inverse use different methods the round trip is only
    ///   approximate (AFNI documents this).
    /// * AFNI quirk, kept: with only two samples the first `x` is returned.
    pub fn inverse_interpolate(&self, y: f64) -> f64 {
        let ar = &self.samples;
        let itop = ar.len() - 1;
        if itop <= 1 {
            return self.x0;
        }
        let first = ar[0];
        let last = ar[itop];
        let increasing = first < last;
        let decreasing = first > last;
        if (increasing && y <= first) || (decreasing && y >= first) {
            return self.x0;
        }
        let x_last = self.x0 + self.dx * itop as f64;
        if (increasing && y >= last) || (decreasing && y <= last) {
            return x_last;
        }
        // One regula-falsi step from x0 toward x1, as in AFNI.
        let regula_falsi = |x0: f64, x1: f64| -> f64 {
            let (y0, y1) = (self.interpolate(x0), self.interpolate(x1));
            let dy = y1 - y0;
            if dy == 0.0 || dy.abs() < 0.00666 * ((y - y0).abs() + (y - y1).abs()) {
                x0
            } else {
                x0 + (x1 - x0) / dy * (y - y0)
            }
        };
        for ip in 1..=itop {
            let (ym, yp) = (ar[ip - 1], ar[ip]);
            if (y - ym) * (y - yp) <= 0.0 {
                // `y` is bracketed by samples ip-1 and ip. (A flat segment that
                // equals `y` gives AFNI a 0/0 here; the left end is used.)
                let frac = if yp == ym { 0.0 } else { (y - ym) / (yp - ym) };
                let x0 = self.x0 + self.dx * (ip as f64 - 1.0 + frac);
                let xp = regula_falsi(x0, x0 + 0.05 * self.dx);
                let xm = regula_falsi(x0, x0 - 0.05 * self.dx);
                let residual = |x: f64| (self.interpolate(x) - y).abs();
                // Smallest residual wins; ties keep the earlier candidate.
                let mut best = (x0, residual(x0));
                for candidate in [xm, xp] {
                    let r = residual(candidate);
                    if r < best.1 {
                        best = (candidate, r);
                    }
                }
                return best.0;
            }
        }
        // Not reached for a curve whose end values bracket `y`.
        self.x0 + self.dx * 0.5 * itop as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_curve_exposes_its_axis() {
        let c = ThresholdCurve::new(1.0, 0.5, vec![3.0, 2.0, 1.0]).unwrap();
        assert_eq!(c.len(), 3);
        assert_eq!(c.x_at(2), Some(2.0));
        assert_eq!(c.x_at(3), None);
        assert!(!c.is_empty());
    }

    #[test]
    fn invalid_curves_are_rejected() {
        assert!(ThresholdCurve::new(f64::NAN, 1.0, vec![1.0, 2.0]).is_err());
        assert!(ThresholdCurve::new(0.0, 0.0, vec![1.0, 2.0]).is_err());
        assert!(ThresholdCurve::new(0.0, 1.0, vec![1.0]).is_err());
        assert!(ThresholdCurve::new(0.0, 1.0, vec![1.0, f64::INFINITY]).is_err());
        // A negative spacing is legitimate (a curve over decreasing log10(p)).
        assert!(ThresholdCurve::new(0.0, -1.0, vec![1.0, 2.0]).is_ok());
    }

    fn curve(x0: f64, dx: f64, v: &[f64]) -> ThresholdCurve {
        ThresholdCurve::new(x0, dx, v.to_vec()).unwrap()
    }

    #[test]
    fn interpolation_hits_the_samples_and_clamps_the_ends() {
        let c = curve(1.0, 0.5, &[0.0, 1.0, 4.0, 9.0, 16.0]);
        for i in 0..5 {
            assert_eq!(
                c.interpolate(1.0 + 0.5 * i as f64),
                c.samples()[i],
                "sample {i}"
            );
        }
        assert_eq!(c.interpolate(-100.0), 0.0);
        assert_eq!(c.interpolate(100.0), 16.0);
        assert!(c.interpolate(f64::NAN).is_nan());
    }

    #[test]
    fn interpolation_is_not_catmull_rom_and_never_overshoots() {
        // AFNI's weights at the midpoint of a unit interval (0.1666667 for 1/6).
        let w = cubic_weights(0.5);
        assert!((w[0] - 0.5 * 0.5 * (-1.5) * 0.1666667).abs() < 1e-15);
        assert!(
            (w.iter().sum::<f64>() - 1.0).abs() < 2e-7,
            "weights nearly sum to 1"
        );
        // At the midpoint AFNI's cubic and Catmull-Rom happen to coincide
        // (-1/16, 9/16, 9/16, -1/16), but elsewhere they do not: at an offset of
        // 0.25 Catmull-Rom's first weight is -0.0703125 and AFNI's is -0.0546875.
        let quarter = cubic_weights(0.25);
        assert!((quarter[0] + 0.054_687_5).abs() < 1e-6, "{}", quarter[0]);
        assert!((quarter[0] + 0.070_312_5).abs() > 1e-3);
        // A steep monotone step: the cubic would overshoot, the clamp prevents it.
        let c = curve(0.0, 1.0, &[0.0, 0.0, 0.0, 10.0, 10.0, 10.0]);
        for k in 0..=500 {
            let v = c.interpolate(k as f64 * 0.01);
            assert!((0.0..=10.0).contains(&v), "overshoot {v}");
        }
        // Between two equal neighbours the value is pinned to them.
        assert_eq!(c.interpolate(0.5), 0.0);
    }

    #[test]
    fn a_two_sample_curve_is_constant_like_afni() {
        let c = curve(0.0, 1.0, &[3.0, 8.0]);
        assert_eq!(c.interpolate(0.5), 3.0);
        assert_eq!(c.interpolate(10.0), 3.0);
        assert_eq!(c.inverse_interpolate(5.0), 0.0);
    }

    #[test]
    fn inverse_interpolation_edges_and_round_trip() {
        let c = curve(2.0, 0.5, &[0.0, 0.3, 1.0, 2.2, 4.0, 7.0]);
        // Off the ends.
        assert_eq!(c.inverse_interpolate(-1.0), 2.0);
        assert_eq!(c.inverse_interpolate(0.0), 2.0);
        assert_eq!(c.inverse_interpolate(7.0), 2.0 + 0.5 * 5.0);
        assert_eq!(c.inverse_interpolate(99.0), 4.5);
        // Round trip is close (not exact: forward is cubic, inverse is refined linear).
        for y in [0.1, 0.5, 1.5, 3.0, 5.5] {
            let x = c.inverse_interpolate(y);
            let back = c.interpolate(x);
            assert!((back - y).abs() < 0.02 * (1.0 + y), "{y}: x={x} -> {back}");
        }
        // Decreasing curves work too.
        let d = curve(0.0, 1.0, &[9.0, 6.0, 4.0, 1.0, 0.0]);
        let x = d.inverse_interpolate(5.0);
        assert!((d.interpolate(x) - 5.0).abs() < 1e-3);
        assert_eq!(d.inverse_interpolate(10.0), 0.0);
    }
}
