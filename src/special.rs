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
// The special functions behind every p-value: log-gamma, log-beta, the
// regularized incomplete beta and gamma functions, and the normal tail. Each
// tail function returns the LOGARITHM of both tails (lower and upper), so a
// probability far below the smallest representable `f64` (about 1e-308) is
// still reported accurately instead of becoming zero.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `stats.rs` builds t, F, chi-square, beta, binomial, Poisson, gamma, normal,
//   and correlation p-values out of these functions.
// * Nothing here knows about statistics kinds, tails, or AFNI; these are plain
//   mathematical building blocks, so they are tested against closed forms.
//
// WHY WRITTEN HERE
//
// Rust's standard library has no gamma/beta/erf functions, and the crate keeps
// a tiny dependency list (roadmap: "small dependency surface"). The algorithms
// are the classic ones (Stirling series, Lentz continued fractions, power
// series), documented below, with accuracy checked against closed forms and
// against AFNI's own `nifticdf` library in `tests/stats_conformance.rs`.
//
// NOTATION
//
//   "lower" tail  P = P(X <= x)        "upper" tail  Q = P(X > x) = 1 - P
//   `ln_p`, `ln_q` are natural logs of P and Q.
// ---------------------------------------------------------------------------

//! Special functions with log-space tails.

use std::f64::consts::{LN_2, PI};

/// Relative convergence target for series and continued fractions.
const EPS: f64 = 1.0e-16;
/// Stand-in for zero inside Lentz's continued-fraction algorithm.
const TINY: f64 = 1.0e-300;
/// Iteration cap. Convergence takes roughly `sqrt(parameter)` steps, so this is
/// generous; hitting it returns `None` rather than a wrong number.
const MAX_ITER: usize = 200_000;

/// `ln(1 - exp(l))` for `l <= 0`, accurate for both tiny and near-one `exp(l)`.
///
/// This is how a log-probability of one tail is turned into the other tail
/// without losing precision: `ln_q = ln_1m_exp(ln_p)`.
///
/// Two formulas are needed because `1 - exp(l)` cancels badly when `l` is close
/// to zero, while `ln(1 - e)` wastes precision when `e = exp(l)` is tiny. The
/// switch point `-ln 2` is the standard choice (Maechler, "Accurately Computing
/// log(1 - exp(-|a|))").
pub fn ln_1m_exp(l: f64) -> f64 {
    if l > 0.0 || l.is_nan() {
        f64::NAN
    } else if l == 0.0 {
        f64::NEG_INFINITY
    } else if l > -LN_2 {
        // exp(l) is near 1: use expm1 to keep the small difference exact.
        (-l.exp_m1()).ln()
    } else {
        // exp(l) is small: ln_1p is exact for small arguments.
        (-l.exp()).ln_1p()
    }
}

/// Correction term of Stirling's series for `ln Gamma(x)`:
/// `1/(12x) - 1/(360x^3) + 1/(1260x^5) - ...`. Accurate to about 1e-17 for
/// `x >= 10`.
fn stirling_correction(x: f64) -> f64 {
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    // Horner evaluation in 1/x^2 of the Bernoulli-number series.
    inv * (1.0 / 12.0
        + inv2
            * (-1.0 / 360.0
                + inv2
                    * (1.0 / 1260.0
                        + inv2
                            * (-1.0 / 1680.0
                                + inv2
                                    * (1.0 / 1188.0
                                        + inv2 * (-691.0 / 360_360.0 + inv2 / 156.0))))))
}

/// Natural log of the gamma function for `x > 0` (NaN otherwise).
///
/// Small arguments are shifted up with the recurrence
/// `Gamma(x) = Gamma(x + 1) / x` until `x >= 10`, where Stirling's series is
/// accurate to double precision.
pub fn ln_gamma(x: f64) -> f64 {
    if x.is_nan() || x <= 0.0 {
        return f64::NAN;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }
    // Accumulate the product x (x+1) (x+2) ... so only one log is taken.
    let mut shift_product = 1.0;
    let mut z = x;
    while z < 10.0 {
        shift_product *= z;
        z += 1.0;
    }
    // Stirling: ln Gamma(z) = (z - 1/2) ln z - z + (1/2) ln(2 pi) + correction.
    let stirling = (z - 0.5) * z.ln() - z + 0.5 * (2.0 * PI).ln() + stirling_correction(z);
    stirling - shift_product.ln()
}

/// `ln Gamma(a) - ln Gamma(a + b)` for `a >= 10`, `b >= 0`, computed from
/// Stirling's series *analytically* so the huge cancelling terms never appear.
fn ln_gamma_ratio_large(a: f64, b: f64) -> f64 {
    // (a - 1/2) ln a - a - [(a + b - 1/2) ln(a + b) - (a + b)] + corrections
    //   = -(a - 1/2) ln(1 + b/a) - b ln(a + b) + b + corrections
    -(a - 0.5) * (b / a).ln_1p() - b * (a + b).ln() + b + stirling_correction(a)
        - stirling_correction(a + b)
}

/// `ln B(a, b) = ln Gamma(a) + ln Gamma(b) - ln Gamma(a + b)` for `a, b > 0`.
pub fn ln_beta(a: f64, b: f64) -> f64 {
    let (small, large) = if a < b { (a, b) } else { (b, a) };
    if large >= 10.0 {
        // ln Gamma(large) - ln Gamma(large + small) without cancellation.
        ln_gamma(small) + ln_gamma_ratio_large(large, small)
    } else {
        ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
    }
}

/// `ln(x)` that stays accurate when `x` is near 1, given `y = 1 - x`.
fn ln_with_complement(x: f64, y: f64) -> f64 {
    if x > 0.5 {
        (-y).ln_1p()
    } else {
        x.ln()
    }
}

/// Continued fraction for the incomplete beta function (Lentz's method; see
/// Numerical Recipes `betacf`). Returns `None` if it fails to converge.
fn beta_continued_fraction(a: f64, b: f64, x: f64) -> Option<f64> {
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=MAX_ITER {
        let m = m as f64;
        let m2 = 2.0 * m;
        // Even step of the recurrence.
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        h *= d * c;
        // Odd step.
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let delta = d * c;
        h *= delta;
        if (delta - 1.0).abs() < EPS {
            return Some(h);
        }
    }
    None
}

/// Log of both tails of the regularized incomplete beta function
/// `I_x(a, b)`, returned as `(ln I, ln (1 - I))`.
///
/// `y` must equal `1 - x`; passing it separately lets callers who know `y`
/// exactly (for example `t^2 / (nu + t^2)`) avoid the rounding in `1 - x`.
/// Requires `a > 0`, `b > 0`, `x` and `y` in `[0, 1]`. Returns `None` if the
/// continued fraction does not converge.
///
/// The continued fraction converges fastest when `x < (a + 1) / (a + b + 2)`;
/// otherwise the symmetry `I_x(a, b) = 1 - I_y(b, a)` is used. Either way the
/// *small* tail is computed directly and the other is obtained with
/// [`ln_1m_exp`], so both stay accurate.
pub fn inc_beta_ln(x: f64, y: f64, a: f64, b: f64) -> Option<(f64, f64)> {
    if x <= 0.0 {
        return Some((f64::NEG_INFINITY, 0.0));
    }
    if y <= 0.0 {
        return Some((0.0, f64::NEG_INFINITY));
    }
    let ln_x = ln_with_complement(x, y);
    let ln_y = ln_with_complement(y, x);
    let ln_front = a * ln_x + b * ln_y - ln_beta(a, b);
    if x < (a + 1.0) / (a + b + 2.0) {
        let ln_p = ln_front - a.ln() + beta_continued_fraction(a, b, x)?.ln();
        Some((ln_p, ln_1m_exp(ln_p)))
    } else {
        let ln_q = ln_front - b.ln() + beta_continued_fraction(b, a, y)?.ln();
        Some((ln_1m_exp(ln_q), ln_q))
    }
}

/// Log of both tails of the regularized incomplete gamma function:
/// `(ln P(a, x), ln Q(a, x))` with `P = gamma(a, x) / Gamma(a)` and `Q = 1 - P`.
///
/// Requires `a > 0` and `x >= 0`. Uses the power series for `x < a + 1` and a
/// continued fraction otherwise (Numerical Recipes `gser`/`gcf`), each in log
/// form. Returns `None` on non-convergence.
pub fn inc_gamma_ln(a: f64, x: f64) -> Option<(f64, f64)> {
    if x <= 0.0 {
        return Some((f64::NEG_INFINITY, 0.0));
    }
    if x.is_infinite() {
        return Some((0.0, f64::NEG_INFINITY));
    }
    // The factor common to both expansions: x^a e^-x / Gamma(a).
    let ln_front = a * x.ln() - x - ln_gamma(a);
    if x < a + 1.0 {
        // Series: P = front * sum_{n>=0} x^n / (a (a+1) ... (a+n)).
        let mut term = 1.0 / a;
        let mut sum = term;
        let mut ap = a;
        for _ in 0..MAX_ITER {
            ap += 1.0;
            term *= x / ap;
            sum += term;
            if term.abs() < sum.abs() * EPS {
                let ln_p = ln_front + sum.ln();
                return Some((ln_p, ln_1m_exp(ln_p.min(0.0))));
            }
        }
        None
    } else {
        // Continued fraction for Q (modified Lentz).
        let mut b = x + 1.0 - a;
        let mut c = 1.0 / TINY;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..=MAX_ITER {
            let i = i as f64;
            let an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if d.abs() < TINY {
                d = TINY;
            }
            c = b + an / c;
            if c.abs() < TINY {
                c = TINY;
            }
            d = 1.0 / d;
            let delta = d * c;
            h *= delta;
            if (delta - 1.0).abs() < EPS {
                let ln_q = ln_front + h.ln();
                return Some((ln_1m_exp(ln_q.min(0.0)), ln_q));
            }
        }
        None
    }
}

/// Log of both tails of the standard normal distribution:
/// `(ln Phi(z), ln (1 - Phi(z)))`.
///
/// Uses the identity `Q(z) = (1/2) Q_gamma(1/2, z^2 / 2)` for `z >= 0` and
/// symmetry for `z < 0`, so the result is accurate far into the tail
/// (`z = 30` gives `ln Q ~ -454`, where `Q` itself underflows).
pub fn normal_ln_tails(z: f64) -> Option<(f64, f64)> {
    if z.is_nan() {
        return None;
    }
    if z.is_infinite() {
        return Some(if z > 0.0 {
            (0.0, f64::NEG_INFINITY)
        } else {
            (f64::NEG_INFINITY, 0.0)
        });
    }
    let (_, ln_q_gamma) = inc_gamma_ln(0.5, 0.5 * z * z)?;
    // Half of the two-sided gamma tail is the one-sided normal tail.
    let ln_small = ln_q_gamma - LN_2;
    let ln_big = ln_1m_exp(ln_small.min(0.0));
    // For z >= 0 the upper tail is the small one; for z < 0 it is the big one.
    Some(if z >= 0.0 {
        (ln_big, ln_small)
    } else {
        (ln_small, ln_big)
    })
}

/// The `x` for which the standard normal upper tail `P(Z > x)` equals
/// `exp(ln_p)`, for `ln_p` in `(-inf, 0)`; `None` if `ln_p` is outside that
/// range or the iteration fails.
///
/// This is AFNI's `qginv` strategy (Abramowitz & Stegun 26.2.23 as a starting
/// guess, then Newton refinement), but refined on the *logarithm* of the tail
/// so it stays accurate for probabilities far below `1e-300`, and run to
/// convergence rather than for a fixed three steps. It is fast enough to call
/// once per voxel.
pub fn normal_upper_quantile_ln(ln_p: f64) -> Option<f64> {
    if ln_p.is_nan() || ln_p >= 0.0 || ln_p == f64::NEG_INFINITY {
        return None;
    }
    // Work with the smaller tail (<= 1/2) and use symmetry for the other.
    let (negate, ln_dp) = if ln_p > -LN_2 {
        (true, ln_1m_exp(ln_p))
    } else {
        (false, ln_p)
    };
    // A&S 26.2.23: |error in x| < 4.5e-4, good enough to start Newton.
    let t = (-2.0 * ln_dp).sqrt();
    let mut x = t
        - ((0.010328 * t + 0.802853) * t + 2.515517)
            / (((0.001308 * t + 0.189269) * t + 1.432788) * t + 1.0);
    // Newton on f(x) = ln Q(x) - ln dp, with f'(x) = -hazard(x) = -phi(x)/Q(x).
    let ln_sqrt_2pi = 0.5 * (2.0 * PI).ln();
    for _ in 0..12 {
        let (_, ln_q) = normal_ln_tails(x)?;
        let hazard = (-0.5 * x * x - ln_sqrt_2pi - ln_q).exp();
        let step = (ln_q - ln_dp) / hazard;
        x += step;
        if step.abs() <= 1e-15 * x.abs().max(1.0) {
            break;
        }
    }
    Some(if negate { -x } else { x })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a - b).abs() <= rel * b.abs().max(1e-300)
    }

    #[test]
    fn ln_1m_exp_matches_direct_formula_and_edges() {
        for l in [-0.1_f64, -0.69, -0.7, -3.0] {
            let direct = (1.0 - l.exp()).ln();
            assert!(close(ln_1m_exp(l), direct, 1e-9), "{l}");
        }
        // Deep negative: ln(1 - e^-50) = -e^-50 to full precision (the naive form gives 0).
        assert!(close(ln_1m_exp(-50.0), -1.928_749_847_963_917_8e-22, 1e-14));
        assert_eq!(ln_1m_exp(0.0), f64::NEG_INFINITY);
        assert_eq!(ln_1m_exp(f64::NEG_INFINITY), 0.0);
        assert!(ln_1m_exp(0.5).is_nan());
        // Tiny argument: 1 - e^-1e-20 = 1e-20 exactly, which the naive form loses.
        assert!(close(ln_1m_exp(-1e-20), (1e-20_f64).ln(), 1e-12));
    }

    #[test]
    fn ln_gamma_matches_known_values() {
        assert!(close(ln_gamma(0.5), 0.5 * PI.ln(), 1e-14)); // ln sqrt(pi)
        assert!(ln_gamma(1.0).abs() < 1e-14);
        assert!(ln_gamma(2.0).abs() < 1e-14);
        assert!(close(ln_gamma(10.0), 362_880_f64.ln(), 1e-14));
        assert!(close(ln_gamma(5.5), 3.957_813_967_618_717, 1e-13));
        assert!(close(ln_gamma(100.0), 359.134_205_369_575_4, 1e-14));
        // Recurrence Gamma(x+1) = x Gamma(x) across the Stirling switch at 10.
        for x in [0.3, 1.7, 8.9, 9.99, 10.0, 25.5, 1e4] {
            assert!(
                (ln_gamma(x + 1.0) - ln_gamma(x) - x.ln()).abs()
                    < 1e-12 * ln_gamma(x + 1.0).abs().max(1.0),
                "{x}"
            );
        }
        assert!(ln_gamma(0.0).is_nan() && ln_gamma(-1.0).is_nan());
    }

    #[test]
    fn ln_beta_is_symmetric_and_continuous_across_methods() {
        assert!(ln_beta(1.0, 1.0).abs() < 1e-14); // B(1,1) = 1
        assert!(close(ln_beta(2.0, 3.0), (1.0_f64 / 12.0).ln(), 1e-13));
        for (a, b) in [(0.5, 20.0), (3.0, 1e5), (9.9, 10.1), (50.0, 0.25)] {
            assert!(close(ln_beta(a, b), ln_beta(b, a), 1e-14));
            let naive = ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b);
            assert!(close(ln_beta(a, b), naive, 1e-8), "({a},{b})");
        }
    }

    #[test]
    fn inc_beta_matches_closed_forms() {
        // I_x(1, 1) = x.
        for x in [0.01, 0.3, 0.5, 0.9] {
            let (p, q) = inc_beta_ln(x, 1.0 - x, 1.0, 1.0).unwrap();
            assert!(
                close(p.exp(), x, 1e-13) && close(q.exp(), 1.0 - x, 1e-13),
                "{x}"
            );
        }
        // I_x(1, b) = 1 - (1-x)^b ; upper tail is exactly (1-x)^b.
        let (x, b) = (0.2, 7.0);
        let (_, q) = inc_beta_ln(x, 1.0 - x, 1.0, b).unwrap();
        assert!(close(q, b * (1.0 - x).ln(), 1e-13));
        // I_x(a, 1) = x^a.
        let (p, _) = inc_beta_ln(0.3, 0.7, 4.0, 1.0).unwrap();
        assert!(close(p, 4.0 * 0.3_f64.ln(), 1e-13));
        // Edges.
        assert_eq!(
            inc_beta_ln(0.0, 1.0, 2.0, 3.0).unwrap(),
            (f64::NEG_INFINITY, 0.0)
        );
        assert_eq!(
            inc_beta_ln(1.0, 0.0, 2.0, 3.0).unwrap(),
            (0.0, f64::NEG_INFINITY)
        );
        // Symmetry I_x(a,b) = 1 - I_{1-x}(b,a).
        let (p1, q1) = inc_beta_ln(0.37, 0.63, 2.5, 6.0).unwrap();
        let (p2, q2) = inc_beta_ln(0.63, 0.37, 6.0, 2.5).unwrap();
        assert!(close(p1, q2, 1e-12) && close(q1, p2, 1e-12));
    }

    #[test]
    fn inc_gamma_matches_closed_forms() {
        // a = 1: P = 1 - e^-x, Q = e^-x.
        for x in [0.1, 1.0, 5.0, 40.0] {
            let (p, q) = inc_gamma_ln(1.0, x).unwrap();
            assert!(close(q, -x, 1e-13), "{x}");
            if x < 30.0 {
                assert!(close(p.exp(), 1.0 - (-x).exp(), 1e-12), "{x}");
            }
        }
        // a = 2: Q = (1 + x) e^-x.
        let x = 3.0;
        let (_, q) = inc_gamma_ln(2.0, x).unwrap();
        assert!(close(q.exp(), (1.0 + x) * (-x).exp(), 1e-13));
        assert_eq!(inc_gamma_ln(2.0, 0.0).unwrap(), (f64::NEG_INFINITY, 0.0));
    }

    #[test]
    fn normal_tails_match_known_quantiles_and_the_deep_tail() {
        let (p0, q0) = normal_ln_tails(0.0).unwrap();
        assert!(close(p0, -LN_2, 1e-14) && close(q0, -LN_2, 1e-14));
        // Phi-bar(1.959963984540054) = 0.025.
        let (_, q) = normal_ln_tails(1.959_963_984_540_054).unwrap();
        assert!(close(q.exp(), 0.025, 1e-12));
        // Symmetry.
        let (pm, qm) = normal_ln_tails(-1.3).unwrap();
        let (pp, qp) = normal_ln_tails(1.3).unwrap();
        assert!(close(pm, qp, 1e-13) && close(qm, pp, 1e-13));
        // Deep tail against the Mills-ratio asymptotic series, where Q underflows.
        let z = 40.0_f64;
        let (_, q) = normal_ln_tails(z).unwrap();
        let asymptotic = -0.5 * z * z - z.ln() - 0.5 * (2.0 * PI).ln()
            + (1.0 - 1.0 / z.powi(2) + 3.0 / z.powi(4) - 15.0 / z.powi(6)).ln();
        assert!((q - asymptotic).abs() < 1e-7, "{q} vs {asymptotic}");
        assert_eq!(
            q.exp(),
            0.0,
            "the plain probability underflows; the log does not"
        );
        assert!(normal_ln_tails(f64::NAN).is_none());
    }

    #[test]
    fn normal_quantile_inverts_the_tail_across_many_magnitudes() {
        // Known: upper tail 0.025 at 1.959963984540054; 0.5 at 0.
        assert!(
            (normal_upper_quantile_ln(0.025_f64.ln()).unwrap() - 1.959_963_984_540_054).abs()
                < 1e-12
        );
        assert!(normal_upper_quantile_ln(0.5_f64.ln()).unwrap().abs() < 1e-12);
        // Round trip from far in the tail (where p underflows) to the far other side.
        for ln_p in [-1.0e-8, -0.1, -0.5, -2.0, -20.0, -200.0, -2000.0, -50_000.0] {
            let x = normal_upper_quantile_ln(ln_p).unwrap();
            let (_, back) = normal_ln_tails(x).unwrap();
            assert!(
                (back - ln_p).abs() < 1e-9 * ln_p.abs().max(1.0),
                "{ln_p}: x={x} -> {back}"
            );
        }
        // Symmetry: p and 1 - p give opposite signs.
        let a = normal_upper_quantile_ln(0.1_f64.ln()).unwrap();
        let b = normal_upper_quantile_ln(0.9_f64.ln()).unwrap();
        assert!((a + b).abs() < 1e-12);
        assert!(
            normal_upper_quantile_ln(0.0).is_none() && normal_upper_quantile_ln(f64::NAN).is_none()
        );
    }
}
