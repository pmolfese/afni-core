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
// Statistical interpretation: given a `StatSpec` (for example `Ttest(23)`) and
// a statistic value, compute the probability that value implies (`p_value`),
// and the reverse (`critical_value`).
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `stat.rs` says WHAT a statistic is (kind + parameters); this file says what
//   it MEANS. `special.rs` supplies the incomplete beta/gamma/normal functions.
// * `error.rs` carries the typed errors. Phase 3 (FDR) and Phase 5 (thresholds
//   and overlays) call into this file to turn thresholds into p-values.
// * Ported from, and checked against, AFNI's `thd_statpval.c`, `mri_stats.c`
//   and the NIfTI `nifticdf.c` library (see tests/stats_conformance.rs).
//
// DESIGN DECISIONS (each corrects something in older AFNI/sumaru behaviour)
//
// 1. The TAIL IS ALWAYS EXPLICIT (`Tail`). AFNI silently picks two-sided for
//    t/z/correlation and upper-tail for F/chi-square/etc. There is no global
//    "two-sided p-value" function here.
// 2. PROBABILITIES KEEP THEIR LOGARITHM (`Probability`), so p = 1e-400 is
//    representable and displayable even though no `f64` can hold it.
// 3. NO SENTINELS. AFNI returns 0, -1, 99.99 or 1.0 on error. Here every
//    failure is an `Error`: bad parameters, non-finite input, statistics outside
//    the distribution's domain, unsupported distributions.
// 4. NO INDISCRIMINATE `abs()`. A two-sided probability is only defined for
//    distributions symmetric about a centre (t, z, signed correlation, normal,
//    logistic, Laplace, uniform); for F, chi-square, gamma, ... asking for
//    `TwoSided` is an error rather than a guess.
// 5. DISCRETE DISTRIBUTIONS (binomial, Poisson) use the true step function:
//    a non-integer statistic `x` counts as `floor(x)`. AFNI/CDFLIB instead
//    evaluates a continuous extension, so they differ at non-integer inputs.
// 6. The gamma distribution's second parameter is a RATE (density
//    ~ x^(shape-1) e^(-rate x)), even though AFNI calls it "scale". This
//    matches CDFLIB and the NIfTI definition.
// 7. NONCENTRAL t/F/chi-square are computed exactly (Poisson mixtures of the
//    central forms; Lenth's AS 243 for t), never approximated by the central
//    distribution, and checked against AFNI's nifticdf library.
//
// TAIL CONVENTIONS
//
//   Lower     P(X <= x)
//   Upper     P(X >  x)
//   TwoSided  P(|X - c| >= |x - c|) = min(1, 2 * min(Lower, Upper)),
//             where c is the distribution's centre of symmetry.
//
// Direct probability statistics (`Pval`, `LogPval`, `Log10Pval`) store a
// probability already; the tail was fixed by whoever produced the data and
// cannot be reinterpreted, so the `tail` argument is ignored for them.
// ---------------------------------------------------------------------------

//! Tail-aware p-values and critical values for AFNI statistics.

use std::f64::consts::LN_2;

use crate::error::{Error, Result};
use crate::numeric::ensure_finite;
use crate::special::{inc_beta_ln, inc_gamma_ln, ln_1m_exp, ln_gamma, normal_ln_tails};
use crate::stat::{StatKind, StatSpec};

/// Which tail of the distribution a probability refers to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Tail {
    /// `P(X <= x)`.
    Lower,
    /// `P(X > x)`.
    Upper,
    /// `P(|X - c| >= |x - c|)` for a distribution symmetric about `c`.
    /// An error for asymmetric distributions.
    TwoSided,
}

/// A probability that remembers its natural logarithm.
///
/// `p()` is a best-effort ordinary `f64`; it underflows to `0.0` below about
/// `1e-308`. `ln_p()` does not, so very small probabilities can still be
/// compared, displayed (`log10_p`), and inverted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Probability {
    ln_p: f64,
}

impl Probability {
    /// From a probability in `[0, 1]`.
    pub fn from_p(p: f64) -> Result<Self> {
        ensure_finite("probability", p)?;
        if !(0.0..=1.0).contains(&p) {
            return Err(Error::InvalidParameter {
                name: "probability".into(),
                reason: format!("{p} is outside [0, 1]"),
            });
        }
        Ok(Self { ln_p: p.ln() })
    }

    /// From the natural log of a probability (`<= 0`, `-inf` allowed).
    pub fn from_ln_p(ln_p: f64) -> Result<Self> {
        if ln_p.is_nan() || ln_p > 0.0 {
            return Err(Error::InvalidParameter {
                name: "ln(probability)".into(),
                reason: format!("{ln_p} is not in (-inf, 0]"),
            });
        }
        Ok(Self { ln_p })
    }

    /// The probability as a plain float (may underflow to 0).
    pub fn p(self) -> f64 {
        self.ln_p.exp()
    }

    /// The natural log of the probability.
    pub fn ln_p(self) -> f64 {
        self.ln_p
    }

    /// `log10` of the probability.
    pub fn log10_p(self) -> f64 {
        self.ln_p / std::f64::consts::LN_10
    }

    /// Whether this probability is at most `threshold`: the test a viewer applies
    /// to decide if a sample passes a "p <= threshold" display setting.
    ///
    /// Compares logarithms, so it stays correct below `1e-308`.
    ///
    /// Two different questions share the letter p, and this type keeps them
    /// apart by how it is used:
    ///
    /// * the exact probability of ONE sample's statistic ("this voxel has
    ///   p = 0.003") is what [`p_value`] returns, and
    /// * a THRESHOLD probability ("show voxels with p <= 0.01") is what a user
    ///   chooses; a sample passes when its probability `is_at_most` the
    ///   threshold, and [`critical_value`] converts the threshold into the
    ///   statistic value that a viewer compares against.
    pub fn is_at_most(self, threshold: Probability) -> bool {
        self.ln_p <= threshold.ln_p
    }

    /// Whether `p()` underflowed to zero although the probability is not
    /// exactly zero.
    pub fn underflows(self) -> bool {
        self.ln_p.is_finite() && self.p() == 0.0
    }
}

// ---------------------------------------------------------------------------
// The validated distribution
// ---------------------------------------------------------------------------

/// A `StatSpec` whose parameters have been checked and decoded into the form
/// the formulas use. Building one is where "validate before calculating"
/// happens, so the math below can assume sane parameters.
#[derive(Debug, Clone, Copy)]
enum Model {
    /// Signed correlation (`nfit == 1`): `r` in `[-1, 1]`, symmetric.
    /// `dof` is `samples - nfit - nort`.
    Correl {
        dof: f64,
    },
    /// Multiple correlation (`nfit > 1`): `R` in `[0, 1]`, one-sided.
    MultiCorrel {
        dof: f64,
        nfit: f64,
    },
    T {
        dof: f64,
    },
    F {
        num: f64,
        den: f64,
    },
    Z,
    ChiSq {
        dof: f64,
    },
    Beta {
        a: f64,
        b: f64,
    },
    Binomial {
        n: f64,
        p: f64,
    },
    Gamma {
        shape: f64,
        rate: f64,
    },
    Poisson {
        lambda: f64,
    },
    Normal {
        mean: f64,
        sd: f64,
    },
    Logistic {
        loc: f64,
        scale: f64,
    },
    Laplace {
        loc: f64,
        scale: f64,
    },
    Uniform {
        start: f64,
        end: f64,
    },
    Weibull {
        loc: f64,
        scale: f64,
        power: f64,
    },
    Chi {
        dof: f64,
    },
    InvGauss {
        mu: f64,
        lambda: f64,
    },
    Extval {
        loc: f64,
        scale: f64,
    },
    /// Noncentral chi-square: `dof` degrees of freedom, noncentrality `lambda`.
    NoncChi {
        dof: f64,
        lambda: f64,
    },
    /// Noncentral F with noncentrality `lambda`.
    NoncF {
        num: f64,
        den: f64,
        lambda: f64,
    },
    /// Noncentral t with noncentrality `delta` (the mean shift; either sign).
    NoncT {
        dof: f64,
        delta: f64,
    },
    /// The stored value is itself a probability (`Pval`, `LogPval`,
    /// `Log10Pval`).
    Direct(StatKind),
}

/// Require `value > 0` and finite.
fn positive(name: &str, value: f64) -> Result<f64> {
    ensure_finite(name, value)?;
    if value > 0.0 {
        Ok(value)
    } else {
        Err(Error::InvalidParameter {
            name: name.into(),
            reason: format!("must be positive, got {value}"),
        })
    }
}

/// Require `value >= 0` and finite.
fn non_negative(name: &str, value: f64) -> Result<f64> {
    ensure_finite(name, value)?;
    if value >= 0.0 {
        Ok(value)
    } else {
        Err(Error::InvalidParameter {
            name: name.into(),
            reason: format!("must be non-negative, got {value}"),
        })
    }
}

impl Model {
    /// Validate `spec` and decode it.
    fn from_spec(spec: &StatSpec) -> Result<Self> {
        let k = spec.kind;
        let q = &spec.params;
        if q.len() != k.num_params() {
            return Err(Error::InvalidParameter {
                name: format!("{} parameters", k.name()),
                reason: format!("expected {}, found {}", k.num_params(), q.len()),
            });
        }
        Ok(match k {
            StatKind::Correl => {
                let (samples, nfit, nort) = (q[0], q[1], q[2]);
                ensure_finite("correlation samples", samples)?;
                ensure_finite("correlation nfit", nfit)?;
                ensure_finite("correlation nort", nort)?;
                if nfit < 1.0 || nort < 0.0 {
                    return Err(Error::InvalidParameter {
                        name: "Correl".into(),
                        reason: format!(
                            "need nfit >= 1 and nort >= 0, got nfit={nfit}, nort={nort}"
                        ),
                    });
                }
                let dof = samples - nfit - nort;
                if dof <= 0.0 {
                    return Err(Error::InvalidParameter {
                        name: "Correl".into(),
                        reason: format!(
                            "samples ({samples}) must exceed nfit + nort ({})",
                            nfit + nort
                        ),
                    });
                }
                if nfit == 1.0 {
                    Model::Correl { dof }
                } else {
                    Model::MultiCorrel { dof, nfit }
                }
            }
            StatKind::Ttest => Model::T {
                dof: positive("t dof", q[0])?,
            },
            StatKind::Ftest => Model::F {
                num: positive("F numerator dof", q[0])?,
                den: positive("F denominator dof", q[1])?,
            },
            StatKind::Zscore => Model::Z,
            StatKind::Chisq => Model::ChiSq {
                dof: positive("chi-square dof", q[0])?,
            },
            StatKind::Beta => Model::Beta {
                a: positive("beta a", q[0])?,
                b: positive("beta b", q[1])?,
            },
            StatKind::Binom => {
                let n = positive("binomial trials", q[0])?;
                let p = ensure_finite("binomial probability", q[1])?;
                if n.fract() != 0.0 || n > 9.0e15 {
                    return Err(Error::InvalidParameter {
                        name: "binomial trials".into(),
                        reason: format!("must be a whole number up to 9e15, got {n}"),
                    });
                }
                if !(0.0..=1.0).contains(&p) {
                    return Err(Error::InvalidParameter {
                        name: "binomial probability".into(),
                        reason: format!("{p} is outside [0, 1]"),
                    });
                }
                Model::Binomial { n, p }
            }
            StatKind::Gamma => Model::Gamma {
                shape: positive("gamma shape", q[0])?,
                rate: positive("gamma rate", q[1])?,
            },
            StatKind::Poisson => {
                let lambda = ensure_finite("Poisson mean", q[0])?;
                if lambda < 0.0 {
                    return Err(Error::InvalidParameter {
                        name: "Poisson mean".into(),
                        reason: format!("must be non-negative, got {lambda}"),
                    });
                }
                Model::Poisson { lambda }
            }
            StatKind::Normal => Model::Normal {
                mean: ensure_finite("normal mean", q[0])?,
                sd: positive("normal sd", q[1])?,
            },
            StatKind::Logistic => Model::Logistic {
                loc: ensure_finite("logistic location", q[0])?,
                scale: positive("logistic scale", q[1])?,
            },
            StatKind::Laplace => Model::Laplace {
                loc: ensure_finite("Laplace location", q[0])?,
                scale: positive("Laplace scale", q[1])?,
            },
            StatKind::Uniform => {
                let (start, end) = (
                    ensure_finite("uniform start", q[0])?,
                    ensure_finite("uniform end", q[1])?,
                );
                if end <= start {
                    return Err(Error::InvalidParameter {
                        name: "uniform range".into(),
                        reason: format!("end ({end}) must exceed start ({start})"),
                    });
                }
                Model::Uniform { start, end }
            }
            StatKind::Weibull => Model::Weibull {
                loc: ensure_finite("Weibull location", q[0])?,
                scale: positive("Weibull scale", q[1])?,
                power: positive("Weibull power", q[2])?,
            },
            StatKind::Chi => Model::Chi {
                dof: positive("chi dof", q[0])?,
            },
            StatKind::Invgauss => Model::InvGauss {
                mu: positive("inverse Gaussian mu", q[0])?,
                lambda: positive("inverse Gaussian lambda", q[1])?,
            },
            StatKind::Extval => Model::Extval {
                loc: ensure_finite("extreme-value location", q[0])?,
                scale: positive("extreme-value scale", q[1])?,
            },
            StatKind::ChisqNonc => Model::NoncChi {
                dof: positive("noncentral chi-square dof", q[0])?,
                lambda: non_negative("noncentral chi-square noncentrality", q[1])?,
            },
            StatKind::FtestNonc => Model::NoncF {
                num: positive("noncentral F numerator dof", q[0])?,
                den: positive("noncentral F denominator dof", q[1])?,
                lambda: non_negative("noncentral F noncentrality", q[2])?,
            },
            StatKind::TtestNonc => Model::NoncT {
                dof: positive("noncentral t dof", q[0])?,
                delta: ensure_finite("noncentral t noncentrality", q[1])?,
            },
            StatKind::Pval | StatKind::LogPval | StatKind::Log10Pval => Model::Direct(k),
        })
    }

    /// Centre of symmetry, if the distribution has one (needed for `TwoSided`).
    fn centre(&self) -> Option<f64> {
        match *self {
            Model::Correl { .. } | Model::T { .. } | Model::Z => Some(0.0),
            Model::Normal { mean, .. } => Some(mean),
            Model::Logistic { loc, .. } | Model::Laplace { loc, .. } => Some(loc),
            Model::Uniform { start, end } => Some(0.5 * (start + end)),
            // Noncentral t is symmetric only when it is not shifted.
            Model::NoncT { delta, .. } => (delta == 0.0).then_some(0.0),
            _ => None,
        }
    }

    /// The closed interval on which the distribution has mass.
    fn support(&self) -> (f64, f64) {
        let inf = f64::INFINITY;
        match *self {
            Model::Correl { .. } => (-1.0, 1.0),
            Model::MultiCorrel { .. } => (0.0, 1.0),
            Model::T { .. } | Model::Z | Model::Normal { .. } | Model::Logistic { .. } => {
                (-inf, inf)
            }
            Model::Laplace { .. } | Model::Extval { .. } => (-inf, inf),
            Model::F { .. } | Model::ChiSq { .. } | Model::Gamma { .. } | Model::Chi { .. } => {
                (0.0, inf)
            }
            Model::InvGauss { .. } => (0.0, inf),
            Model::Beta { .. } => (0.0, 1.0),
            Model::Binomial { n, .. } => (0.0, n),
            Model::Poisson { .. } => (0.0, inf),
            Model::Uniform { start, end } => (start, end),
            Model::Weibull { loc, .. } => (loc, inf),
            Model::NoncChi { .. } | Model::NoncF { .. } => (0.0, inf),
            Model::NoncT { .. } => (-inf, inf),
            Model::Direct(_) => (f64::NAN, f64::NAN),
        }
    }

    /// A starting point and scale for bracketing a root.
    fn guess(&self) -> (f64, f64) {
        match *self {
            Model::Correl { .. } | Model::MultiCorrel { .. } => (0.5, 0.25),
            Model::Beta { .. } => (0.5, 0.25),
            Model::Binomial { n, p } => (n * p, (n * p * (1.0 - p)).sqrt().max(1.0)),
            Model::Poisson { lambda } => (lambda, lambda.sqrt().max(1.0)),
            Model::Gamma { shape, rate } => (shape / rate, (shape.sqrt() / rate).max(1e-3)),
            Model::Normal { mean, sd } => (mean, sd),
            Model::Logistic { loc, scale } | Model::Laplace { loc, scale } => (loc, scale),
            Model::Extval { loc, scale } => (loc, scale),
            Model::Uniform { start, end } => (0.5 * (start + end), 0.25 * (end - start)),
            Model::Weibull { loc, scale, .. } => (loc + scale, scale),
            Model::InvGauss { mu, .. } => (mu, mu),
            Model::F { .. } | Model::ChiSq { .. } | Model::Chi { .. } => (1.0, 1.0),
            Model::NoncChi { dof, lambda } => (dof + lambda, (2.0 * (dof + 2.0 * lambda)).sqrt()),
            Model::NoncF { .. } => (1.0, 1.0),
            Model::NoncT { delta, .. } => (delta, 1.0 + delta.abs()),
            _ => (0.0, 1.0),
        }
    }

    fn is_discrete(&self) -> bool {
        matches!(self, Model::Binomial { .. } | Model::Poisson { .. })
    }

    /// `(ln P(X <= x), ln P(X > x))` for a finite `x`.
    fn ln_tails(&self, x: f64) -> Result<(f64, f64)> {
        let fail = || Error::NoSolution(format!("tail probability did not converge at {x}"));
        let neg_inf = f64::NEG_INFINITY;
        Ok(match *self {
            Model::Z => normal_ln_tails(x).ok_or_else(fail)?,
            Model::Normal { mean, sd } => normal_ln_tails((x - mean) / sd).ok_or_else(fail)?,
            Model::T { dof } => symmetric_beta_tails(x, dof, 0.5, None).ok_or_else(fail)?,
            Model::Correl { dof } => {
                if x.abs() > 1.0 {
                    return Err(Error::InvalidParameter {
                        name: "correlation".into(),
                        reason: format!("{x} is outside [-1, 1]"),
                    });
                }
                // Two-sided p = I_{1-r^2}(dof/2, 1/2); see module docs.
                let r2 = x * x;
                let (xb, yb) = ((1.0 - x.abs()) * (1.0 + x.abs()), r2);
                signed_half_tails(x, inc_beta_ln(xb, yb, 0.5 * dof, 0.5).ok_or_else(fail)?)
            }
            Model::MultiCorrel { dof, nfit } => {
                if x > 1.0 {
                    return Err(Error::InvalidParameter {
                        name: "multiple correlation".into(),
                        reason: format!("{x} is above 1"),
                    });
                }
                if x < 0.0 {
                    (neg_inf, 0.0)
                } else {
                    let (xb, yb) = ((1.0 - x) * (1.0 + x), x * x);
                    let (ln_i, ln_ni) =
                        inc_beta_ln(xb, yb, 0.5 * dof, 0.5 * nfit).ok_or_else(fail)?;
                    // The beta tail IS the upper tail of R.
                    (ln_ni, ln_i)
                }
            }
            Model::F { num, den } => {
                if x < 0.0 {
                    (neg_inf, 0.0)
                } else {
                    let (xb, yb) = ratio_complement(den, num * x);
                    let (ln_i, ln_ni) =
                        inc_beta_ln(xb, yb, 0.5 * den, 0.5 * num).ok_or_else(fail)?;
                    (ln_ni, ln_i)
                }
            }
            Model::ChiSq { dof } => {
                if x <= 0.0 {
                    (neg_inf, 0.0)
                } else {
                    inc_gamma_ln(0.5 * dof, 0.5 * x).ok_or_else(fail)?
                }
            }
            Model::Chi { dof } => {
                if x <= 0.0 {
                    (neg_inf, 0.0)
                } else {
                    inc_gamma_ln(0.5 * dof, 0.5 * x * x).ok_or_else(fail)?
                }
            }
            Model::Gamma { shape, rate } => {
                if x <= 0.0 {
                    (neg_inf, 0.0)
                } else {
                    inc_gamma_ln(shape, rate * x).ok_or_else(fail)?
                }
            }
            Model::Beta { a, b } => {
                if x <= 0.0 {
                    (neg_inf, 0.0)
                } else if x >= 1.0 {
                    (0.0, neg_inf)
                } else {
                    inc_beta_ln(x, 1.0 - x, a, b).ok_or_else(fail)?
                }
            }
            Model::Binomial { n, p } => {
                let k = x.floor();
                if k < 0.0 {
                    (neg_inf, 0.0)
                } else if k >= n || p == 0.0 {
                    (0.0, neg_inf)
                } else if p == 1.0 {
                    (neg_inf, 0.0)
                } else {
                    // P(X > k) = I_p(k + 1, n - k).
                    let (ln_i, ln_ni) = inc_beta_ln(p, 1.0 - p, k + 1.0, n - k).ok_or_else(fail)?;
                    (ln_ni, ln_i)
                }
            }
            Model::Poisson { lambda } => {
                let k = x.floor();
                if k < 0.0 {
                    (neg_inf, 0.0)
                } else if lambda == 0.0 {
                    (0.0, neg_inf)
                } else {
                    // P(X <= k) = Q(k + 1, lambda): the regularized upper gamma.
                    let (ln_p_gamma, ln_q_gamma) =
                        inc_gamma_ln(k + 1.0, lambda).ok_or_else(fail)?;
                    (ln_q_gamma, ln_p_gamma)
                }
            }
            Model::Logistic { loc, scale } => {
                let z = (x - loc) / scale;
                // ln Q = -softplus(z), ln P = -softplus(-z).
                let soft = |v: f64| v.max(0.0) + (-v.abs()).exp().ln_1p();
                (-soft(-z), -soft(z))
            }
            Model::Laplace { loc, scale } => {
                let z = (x - loc) / scale;
                let ln_small = -LN_2 - z.abs();
                let ln_big = ln_1m_exp(ln_small);
                if z >= 0.0 {
                    (ln_big, ln_small)
                } else {
                    (ln_small, ln_big)
                }
            }
            Model::Uniform { start, end } => {
                if x <= start {
                    (neg_inf, 0.0)
                } else if x >= end {
                    (0.0, neg_inf)
                } else {
                    (
                        ((x - start) / (end - start)).ln(),
                        ((end - x) / (end - start)).ln(),
                    )
                }
            }
            Model::Weibull { loc, scale, power } => {
                let z = (x - loc) / scale;
                if z <= 0.0 {
                    (neg_inf, 0.0)
                } else {
                    let w = z.powf(power);
                    (ln_1m_exp(-w), -w)
                }
            }
            Model::Extval { loc, scale } => {
                // Type I (Gumbel, maximum): cdf = exp(-exp(-z)).
                let z = (x - loc) / scale;
                let y = (-z).exp();
                let ln_upper = if y < 1.0e-4 {
                    // 1 - exp(-y) = y (1 - y/2 + y^2/6 - y^3/24 ...); keeps
                    // the far tail (-z) even when y itself underflows.
                    -z + (-y / 2.0 + y * y / 6.0 - y * y * y / 24.0).ln_1p()
                } else {
                    ln_1m_exp(-y)
                };
                (-y, ln_upper)
            }
            Model::InvGauss { mu, lambda } => {
                let z = x / mu;
                let c = lambda / mu;
                if z <= 0.0 {
                    (neg_inf, 0.0)
                } else {
                    // cdf = Phi(v1) + e^(2c) Phi(v2). The upper tail is found
                    // by subtraction, so it loses accuracy below ~1e-12.
                    let s = (c / z).sqrt();
                    let ln_a = normal_ln_tails(s * (z - 1.0)).ok_or_else(fail)?.0;
                    let ln_b = normal_ln_tails(-s * (z + 1.0)).ok_or_else(fail)?.0;
                    let ln_lower = ln_add_exp(ln_a, 2.0 * c + ln_b).min(0.0);
                    (ln_lower, ln_1m_exp(ln_lower))
                }
            }
            Model::NoncChi { dof, lambda } => {
                if x <= 0.0 {
                    (neg_inf, 0.0)
                } else {
                    // Poisson(lambda/2) mixture of central chi-squares.
                    poisson_mixture(0.5 * lambda, |j| {
                        inc_gamma_ln(0.5 * dof + j as f64, 0.5 * x)
                    })
                    .ok_or_else(fail)?
                }
            }
            Model::NoncF { num, den, lambda } => {
                if x <= 0.0 {
                    (neg_inf, 0.0)
                } else {
                    // Poisson mixture of central F's, written through the beta
                    // variable u = num*x / (num*x + den).
                    let (yb, xb) = ratio_complement(den, num * x);
                    poisson_mixture(0.5 * lambda, |j| {
                        inc_beta_ln(xb, yb, 0.5 * num + j as f64, 0.5 * den)
                    })
                    .ok_or_else(fail)?
                }
            }
            Model::NoncT { dof, delta } => noncentral_t_tails(x, dof, delta).ok_or_else(fail)?,
            Model::Direct(_) => {
                return Err(Error::Unsupported(
                    "direct probability statistics have no tail structure".into(),
                ))
            }
        })
    }
}

/// `ln(exp(a) + exp(b))` without overflow.
fn ln_add_exp(a: f64, b: f64) -> f64 {
    let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
    if lo == f64::NEG_INFINITY {
        hi
    } else {
        hi + (lo - hi).exp().ln_1p()
    }
}

/// `(d/(d + s), s/(d + s))` without overflow for huge `s`: the two
/// complementary beta arguments used by the t and F distributions.
fn ratio_complement(d: f64, s: f64) -> (f64, f64) {
    if s <= d {
        (d / (d + s), s / (d + s))
    } else {
        let r = d / s; // small
        (r / (1.0 + r), 1.0 / (1.0 + r))
    }
}

/// Turn the two-sided beta result `(ln I, ln(1 - I))` for `|x|` into signed
/// lower/upper tails: the upper tail of a positive `x` is half the two-sided
/// probability, and symmetry handles negative `x`.
fn signed_half_tails(x: f64, (ln_i, _): (f64, f64)) -> (f64, f64) {
    let ln_small = ln_i - LN_2; // P(X > |x|)
    let ln_big = ln_1m_exp(ln_small.min(0.0));
    if x >= 0.0 {
        (ln_big, ln_small)
    } else {
        (ln_small, ln_big)
    }
}

/// Tails of Student's t with `dof` degrees of freedom at `t`.
/// (`_unused` keeps the signature parallel to future non-central variants.)
fn symmetric_beta_tails(t: f64, dof: f64, half: f64, _unused: Option<f64>) -> Option<(f64, f64)> {
    if t == 0.0 {
        return Some((-LN_2, -LN_2));
    }
    let (xb, yb) = ratio_complement(dof, t * t);
    let both = inc_beta_ln(xb, yb, 0.5 * dof, half)?;
    Some(signed_half_tails(t, both))
}

// ---------------------------------------------------------------------------
// Noncentral distributions
// ---------------------------------------------------------------------------
//
// The noncentral chi-square and F distributions are Poisson mixtures of their
// central forms (Johnson, Kotz & Balakrishnan):
//
//     P(X <= x) = sum_j Poisson(j; lambda/2) * P_central(parameter + j; x)
//
// Every term is positive, so lower and upper tails are each a sum of positive
// numbers, accumulated in log space. That keeps both tails accurate (no
// "1 - lower" subtraction) and lets p below 1e-308 survive. The noncentral t
// follows Lenth's Algorithm AS 243, which is a mixture of incomplete betas.

/// `ln` of the Poisson probability mass at `j` for mean `h`.
fn ln_poisson_pmf(h: f64, j: f64) -> f64 {
    if h == 0.0 {
        return if j == 0.0 { 0.0 } else { f64::NEG_INFINITY };
    }
    -h + j * h.ln() - ln_gamma(j + 1.0)
}

/// Most terms any mixture will add in each direction before giving up.
const MAX_MIXTURE_TERMS: u64 = 5_000_000;

/// The smaller of the finite values among `a` and `b` (or `-inf` if neither).
fn min_finite(a: f64, b: f64) -> f64 {
    match (a.is_finite(), b.is_finite()) {
        (true, true) => a.min(b),
        (true, false) => a,
        (false, true) => b,
        _ => f64::NEG_INFINITY,
    }
}

/// `sum_j Poisson(j; h) * term(j)` for both tails, in log space.
///
/// `term(j)` returns `(ln lower_j, ln upper_j)` for the central distribution
/// with the `j`th shifted parameter. Summation starts at the Poisson mode and
/// moves outward in both directions, stopping once the Poisson weight is too
/// small to matter for the *smaller* of the two sums (each term is at most its
/// weight, since both tails are probabilities). Returns `None` if a term fails
/// or the term limit is reached.
fn poisson_mixture(h: f64, term: impl Fn(u64) -> Option<(f64, f64)>) -> Option<(f64, f64)> {
    let ln_eps = f64::EPSILON.ln() - 2.0; // two extra digits of margin
    let (mut ln_lo, mut ln_up) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let add = |j: u64, lo: &mut f64, up: &mut f64| -> Option<f64> {
        let ln_w = ln_poisson_pmf(h, j as f64);
        let (tl, tu) = term(j)?;
        *lo = ln_add_exp(*lo, ln_w + tl);
        *up = ln_add_exp(*up, ln_w + tu);
        Some(ln_w)
    };
    let mode = h.floor() as u64;
    add(mode, &mut ln_lo, &mut ln_up)?;
    // Upward from the mode.
    for step in 1..=MAX_MIXTURE_TERMS {
        let j = mode + step;
        if ln_poisson_pmf(h, j as f64) < ln_eps + min_finite(ln_lo, ln_up) {
            break;
        }
        add(j, &mut ln_lo, &mut ln_up)?;
        if step == MAX_MIXTURE_TERMS {
            return None;
        }
    }
    // Downward from the mode.
    for step in 1..=mode {
        let j = mode - step;
        if ln_poisson_pmf(h, j as f64) < ln_eps + min_finite(ln_lo, ln_up) {
            break;
        }
        add(j, &mut ln_lo, &mut ln_up)?;
    }
    Some((ln_lo.min(0.0), ln_up.min(0.0)))
}

/// Tails of the noncentral t distribution (Lenth, AS 243).
///
/// For `t > 0` and `h = delta^2 / 2`, with `x = t^2 / (t^2 + dof)`,
/// `p_j = e^-h h^j / j!` and `q_j = |delta| e^-h h^j / (sqrt(2) Gamma(j + 3/2))`:
///
/// ```text
/// lower = Phi(-delta) + (1/2) sum_j [ p_j I_x(j + 1/2, dof/2) + s q_j I_x(j + 1, dof/2) ]
/// upper =               (1/2) sum_j [ p_j (1 - I_x(...))      + s q_j (1 - I_x(...)) ]
/// ```
///
/// where `s = sign(delta)`. For `delta >= 0` both are sums of positive terms,
/// so both tails are accurate in log space. For `delta < 0` the `q` terms
/// subtract, which loses digits in the (small) far tail; if so few digits
/// survive that the answer is meaningless, this returns `None`. Negative `t`
/// is handled by the reflection `T -> -T, delta -> -delta`.
fn noncentral_t_tails(t: f64, dof: f64, delta: f64) -> Option<(f64, f64)> {
    let (ln_phi_neg, ln_phi_pos) = {
        let (lo, up) = normal_ln_tails(-delta)?;
        (lo, up) // Phi(-delta), 1 - Phi(-delta) = Phi(delta)
    };
    if t == 0.0 {
        return Some((ln_phi_neg, ln_phi_pos));
    }
    if t < 0.0 {
        // P(T <= t; delta) = P(T' >= -t; -delta) with T' = -T.
        let (lo, up) = noncentral_t_tails(-t, dof, -delta)?;
        return Some((up, lo));
    }
    let a = delta.abs();
    let h = 0.5 * a * a;
    let (xb, yb) = {
        let (d_over, s_over) = ratio_complement(dof, t * t);
        (s_over, d_over) // x = t^2/(t^2+dof), 1-x = dof/(t^2+dof)
    };
    let ln_a = if a > 0.0 { a.ln() } else { f64::NEG_INFINITY };
    // Four positive sums: p-terms and q-terms, for each beta tail.
    let (mut p_lo, mut q_lo) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let (mut p_up, mut q_up) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let ln_eps = f64::EPSILON.ln() - 2.0;
    let ln_half_ln2 = -0.5 * LN_2;
    let step_term = |j: u64, sums: &mut [f64; 4]| -> Option<f64> {
        let jf = j as f64;
        let ln_p = ln_poisson_pmf(h, jf);
        let ln_q = if a > 0.0 {
            ln_a - h + jf * h.ln() + ln_half_ln2 - ln_gamma(jf + 1.5)
        } else {
            f64::NEG_INFINITY
        };
        let (lo_p, up_p) = inc_beta_ln(xb, yb, jf + 0.5, 0.5 * dof)?;
        let (lo_q, up_q) = inc_beta_ln(xb, yb, jf + 1.0, 0.5 * dof)?;
        sums[0] = ln_add_exp(sums[0], ln_p + lo_p);
        sums[1] = ln_add_exp(sums[1], ln_q + lo_q);
        sums[2] = ln_add_exp(sums[2], ln_p + up_p);
        sums[3] = ln_add_exp(sums[3], ln_q + up_q);
        Some(ln_p.max(ln_q))
    };
    let mut sums = [p_lo, q_lo, p_up, q_up];
    let mode = h.floor() as u64;
    step_term(mode, &mut sums)?;
    let smallest = |s: &[f64; 4]| s.iter().fold(f64::NEG_INFINITY, |m, &v| min_finite(m, v));
    let weight = |j: u64| ln_poisson_pmf(h, j as f64);
    // Outward from the mode in both directions, as in `poisson_mixture`.
    for step in 1..=MAX_MIXTURE_TERMS {
        let j = mode + step;
        if weight(j) < ln_eps + smallest(&sums) {
            break;
        }
        step_term(j, &mut sums)?;
        if step == MAX_MIXTURE_TERMS {
            return None;
        }
    }
    for step in 1..=mode {
        let j = mode - step;
        if weight(j) < ln_eps + smallest(&sums) {
            break;
        }
        step_term(j, &mut sums)?;
    }
    [p_lo, q_lo, p_up, q_up] = sums;

    if delta >= 0.0 {
        let ln_lower = ln_add_exp(ln_phi_neg, -LN_2 + ln_add_exp(p_lo, q_lo));
        let ln_upper = -LN_2 + ln_add_exp(p_up, q_up);
        Some((ln_lower.min(0.0), ln_upper.min(0.0)))
    } else {
        // delta < 0: the q terms subtract.
        let lower = ln_phi_neg.exp() + 0.5 * (p_lo.exp() - q_lo.exp());
        if lower <= 0.0 {
            return None;
        }
        let ln_lower = lower.ln().min(0.0);
        // The upper tail is the DIFFERENCE of two positive sums. When they are
        // nearly equal most digits cancel (this is the far tail on the "wrong"
        // side of the shift), so fall back to integrating the definition,
        // whose integrand is positive and has no cancellation.
        let ln_diff = if q_up < p_up {
            -LN_2 + p_up + ln_1m_exp(q_up - p_up)
        } else {
            f64::NEG_INFINITY
        };
        let digits_lost = (-LN_2 + p_up) - ln_diff; // ln of the cancellation factor
        let ln_upper = if digits_lost.is_finite() && digits_lost < 1.0e6_f64.ln() {
            ln_diff
        } else {
            noncentral_t_upper_by_quadrature(t, dof, delta)?
        };
        Some((ln_lower, ln_upper.min(0.0)))
    }
}

/// `ln P(T > t)` for the noncentral t with `t > 0`, by integrating the definition
///
/// ```text
/// P(T > t) = E[ Q(t W - delta) ],   W = sqrt(V / dof),  V ~ chi-square(dof),
/// ```
///
/// where `Q` is the standard normal upper tail. Unlike the AS 243 series this
/// has a positive integrand for either sign of `delta`, so it stays accurate
/// when the series would cancel. Substituting `w = e^x` makes the integrand
/// smooth and log-concave in `x` for every `dof > 0`; it is then summed with the
/// trapezoid rule over the region within `e^-60` of its peak, which converges
/// exponentially for such integrands.
fn noncentral_t_upper_by_quadrature(t: f64, dof: f64, delta: f64) -> Option<f64> {
    let half = 0.5 * dof;
    // ln of the density of W = chi / sqrt(dof), plus the Jacobian of w = e^x.
    let ln_const = LN_2 + half * half.ln() - ln_gamma(half);
    let log_integrand = |x: f64| -> Option<f64> {
        let w = x.exp();
        let ln_q = normal_ln_tails(t * w - delta)?.1;
        // w^dof e^(-dof w^2 / 2): the density's w^(dof-1) times the Jacobian w.
        Some(ln_const + dof * x - half * w * w + ln_q)
    };
    // Locate the peak (the integrand is concave in x): golden-section search.
    let (mut lo, mut hi) = (-(60.0 / dof) - 5.0, (1.0 + 20.0 / dof.sqrt()).ln() + 1.0);
    lo = lo.max(-700.0);
    let golden = 0.618_033_988_749_894_9;
    let (mut a, mut b) = (hi - golden * (hi - lo), lo + golden * (hi - lo));
    let (mut fa, mut fb) = (log_integrand(a)?, log_integrand(b)?);
    for _ in 0..120 {
        if fa < fb {
            lo = a;
            a = b;
            fa = fb;
            b = lo + golden * (hi - lo);
            fb = log_integrand(b)?;
        } else {
            hi = b;
            b = a;
            fb = fa;
            a = hi - golden * (hi - lo);
            fa = log_integrand(a)?;
        }
    }
    let peak = 0.5 * (lo + hi);
    let l_peak = log_integrand(peak)?;
    // Walk outward until the integrand has fallen by 60 e-folds on each side.
    let edge = |direction: f64| -> Option<f64> {
        let mut step = 0.05_f64;
        let mut x = peak;
        for _ in 0..80 {
            x += direction * step;
            if x < -700.0 {
                return Some(-700.0);
            }
            if log_integrand(x)? < l_peak - 60.0 {
                return Some(x);
            }
            step *= 1.5;
        }
        None
    };
    let (left, right) = (edge(-1.0)?, edge(1.0)?);
    const POINTS: usize = 4001;
    let h = (right - left) / (POINTS - 1) as f64;
    let mut sum = 0.0;
    for i in 0..POINTS {
        let x = left + h * i as f64;
        let weight = if i == 0 || i == POINTS - 1 { 0.5 } else { 1.0 };
        sum += weight * (log_integrand(x)? - l_peak).exp();
    }
    Some(l_peak + (h * sum).ln())
}

// ---------------------------------------------------------------------------
// Forward: statistic -> probability
// ---------------------------------------------------------------------------

/// Probability implied by `statistic` under `spec`, for the requested `tail`.
///
/// # Errors
///
/// * [`Error::NonFinite`] if `statistic` is NaN or infinite.
/// * [`Error::InvalidParameter`] for bad distribution parameters or a statistic
///   outside the distribution's domain (for example a correlation above 1).
/// * [`Error::InvalidParameter`] if `tail` is `TwoSided` for a distribution
///   with no centre of symmetry (F, chi-square, gamma, beta, binomial, ...).
/// * [`Error::Unsupported`] for noncentral distributions.
///
/// For `Pval`, `LogPval` and `Log10Pval` the value is already a probability and
/// `tail` is ignored: `Pval` uses the value itself (must lie in `[0, 1]`),
/// `LogPval` uses `exp(-|value|)` and `Log10Pval` uses `10^(-|value|)`.
///
/// ```
/// use afni_core::stat::{StatKind, StatSpec};
/// use afni_core::stats::{p_value, Tail};
///
/// let t = StatSpec::new(StatKind::Ttest, &[10.0], 0.0);
/// // Two-sided and one-sided probabilities are different questions:
/// let two = p_value(&t, 2.0, Tail::TwoSided)?.p();
/// let one = p_value(&t, 2.0, Tail::Upper)?.p();
/// assert!((two - 2.0 * one).abs() < 1e-15);
/// assert!((two - 0.073388).abs() < 1e-6); // matches AFNI's `cdf -t2p fitt 2.0 10`
/// # Ok::<(), afni_core::Error>(())
/// ```
pub fn p_value(spec: &StatSpec, statistic: f64, tail: Tail) -> Result<Probability> {
    let model = Model::from_spec(spec)?;
    ensure_finite("statistic", statistic)?;
    if let Model::Direct(kind) = model {
        return direct_probability(kind, statistic);
    }
    let (ln_lower, ln_upper) = model.ln_tails(statistic)?;
    let ln_p = match tail {
        Tail::Lower => ln_lower,
        Tail::Upper => ln_upper,
        Tail::TwoSided => {
            require_centre(&model, spec)?;
            LN_2 + ln_lower.min(ln_upper)
        }
    };
    // Rounding can push a log-probability a hair above zero.
    Probability::from_ln_p(ln_p.min(0.0))
}

/// Interpret a stored direct probability.
fn direct_probability(kind: StatKind, value: f64) -> Result<Probability> {
    match kind {
        StatKind::Pval => Probability::from_p(value),
        // exp(-|v|) and 10^(-|v|): the sign is discarded by NIfTI's definition.
        StatKind::LogPval => Probability::from_ln_p(-value.abs()),
        StatKind::Log10Pval => Probability::from_ln_p(-value.abs() * std::f64::consts::LN_10),
        _ => unreachable!("only direct kinds reach here"),
    }
}

fn require_centre(model: &Model, spec: &StatSpec) -> Result<f64> {
    model.centre().ok_or_else(|| Error::InvalidParameter {
        name: "tail".into(),
        reason: format!(
            "{} has no centre of symmetry, so a two-sided probability is not defined; \
                 ask for Tail::Upper or Tail::Lower",
            spec.kind.name()
        ),
    })
}

/// A validated statistic and tail, ready to evaluate many values quickly.
///
/// [`p_value`] re-validates its parameters on every call. When converting a
/// whole dataset (millions of voxels), build one `PValueEvaluator` and call
/// [`evaluate`](Self::evaluate) per value instead.
#[derive(Debug, Clone, Copy)]
pub struct PValueEvaluator {
    model: Model,
    tail: Tail,
}

impl PValueEvaluator {
    /// Validate `spec` and `tail` once. Direct-probability kinds are rejected:
    /// their values are already probabilities.
    pub fn new(spec: &StatSpec, tail: Tail) -> Result<Self> {
        let model = Model::from_spec(spec)?;
        if let Model::Direct(_) = model {
            return Err(Error::Unsupported(format!(
                "{} values are already probabilities; there is nothing to evaluate",
                spec.kind.name()
            )));
        }
        if tail == Tail::TwoSided {
            require_centre(&model, spec)?;
        }
        Ok(Self { model, tail })
    }

    /// The probability for `statistic` (same rules as [`p_value`]).
    pub fn evaluate(&self, statistic: f64) -> Result<Probability> {
        ensure_finite("statistic", statistic)?;
        let (ln_lower, ln_upper) = self.model.ln_tails(statistic)?;
        let ln_p = match self.tail {
            Tail::Lower => ln_lower,
            Tail::Upper => ln_upper,
            Tail::TwoSided => LN_2 + ln_lower.min(ln_upper),
        };
        Probability::from_ln_p(ln_p.min(0.0))
    }
}

// ---------------------------------------------------------------------------
// Inverse: probability -> statistic
// ---------------------------------------------------------------------------

/// The statistic whose tail probability equals `p`.
///
/// Convenience for [`critical_value_ln`] with `ln(p)`.
///
/// Conventions:
///
/// * `p` must lie in `(0, 1)`. (`TwoSided` also accepts `p = 1`, which is the
///   centre.) The endpoints otherwise correspond to infinity or a support edge
///   and are reported as errors rather than invented numbers.
/// * `Upper`: the `x` with `P(X > x) = p`. `Lower`: the `x` with `P(X <= x) = p`.
/// * `TwoSided`: the *upper* critical value, `x >= c`, with
///   `P(|X - c| >= x - c) = p`. The matching lower value is `2c - x`.
/// * Discrete distributions (binomial, Poisson) have no exact inverse. `Upper`
///   returns the smallest integer `k` with `P(X > k) <= p` (the least value
///   that is significant at level `p`); `Lower` returns the largest integer `k`
///   with `P(X <= k) <= p`. AFNI instead returns a non-integer root of a
///   continuous extension.
/// * For direct probability kinds the "statistic" is the stored value: `p`,
///   `-ln p`, or `-log10 p`.
pub fn critical_value(spec: &StatSpec, p: f64, tail: Tail) -> Result<f64> {
    ensure_finite("probability", p)?;
    if !(0.0..=1.0).contains(&p) {
        return Err(Error::InvalidParameter {
            name: "probability".into(),
            reason: format!("{p} is outside [0, 1]"),
        });
    }
    critical_value_ln(spec, p.ln(), tail)
}

/// Like [`critical_value`], taking the natural log of the probability so that
/// probabilities below `1e-308` can be inverted.
pub fn critical_value_ln(spec: &StatSpec, ln_p: f64, tail: Tail) -> Result<f64> {
    let model = Model::from_spec(spec)?;
    if ln_p.is_nan() || ln_p > 0.0 {
        return Err(Error::InvalidParameter {
            name: "ln(probability)".into(),
            reason: format!("{ln_p} is not in (-inf, 0]"),
        });
    }
    if let Model::Direct(kind) = model {
        // The tail cannot change what the stored value means.
        return match kind {
            StatKind::Pval => Ok(ln_p.exp()),
            StatKind::LogPval => finite_value(-ln_p),
            _ => finite_value(-ln_p / std::f64::consts::LN_10),
        };
    }
    // Which one-sided problem to solve, and at what level.
    let (upper, ln_target) = match tail {
        Tail::Upper => (true, ln_p),
        Tail::Lower => (false, ln_p),
        Tail::TwoSided => {
            let centre = require_centre(&model, spec)?;
            if ln_p == 0.0 {
                return Ok(centre); // p = 1: every value is at least as central.
            }
            (true, ln_p - LN_2)
        }
    };
    if ln_target == f64::NEG_INFINITY || (ln_target == 0.0 && tail != Tail::TwoSided) {
        return Err(Error::NoSolution(
            "probabilities of exactly 0 or 1 correspond to infinity or the edge of the \
             support; use a value strictly between 0 and 1"
                .into(),
        ));
    }
    if model.is_discrete() {
        solve_discrete(&model, ln_target, upper)
    } else {
        solve_continuous(&model, ln_target, upper)
    }
}

fn finite_value(v: f64) -> Result<f64> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(Error::NoSolution(
            "the corresponding value is infinite".into(),
        ))
    }
}

/// Find `x` with `ln_tail(x) = ln_target` for a continuous distribution.
///
/// `g(x)` is built to be increasing in `x` whichever tail is used, so one
/// bracket-then-bisect routine serves both: first the interval is grown
/// geometrically from a sensible starting point until `g` changes sign (never
/// leaving the support), then it is bisected to the limit of `f64` precision.
fn solve_continuous(model: &Model, ln_target: f64, upper: bool) -> Result<f64> {
    let g = |x: f64| -> Result<f64> {
        let (ln_lower, ln_upper) = model.ln_tails(x)?;
        Ok(if upper {
            ln_target - ln_upper
        } else {
            ln_lower - ln_target
        })
    };
    let (support_lo, support_hi) = model.support();
    let (start, scale) = model.guess();
    let start = start.clamp(support_lo, support_hi);
    let unreachable = || Error::NoSolution("could not bracket the requested probability".into());

    // Find lo <= hi with g(lo) <= 0 <= g(hi).
    let (mut lo, mut hi) = (start, start);
    let mut step = scale.abs().max(1e-3);
    if g(start)? <= 0.0 {
        for _ in 0..2100 {
            hi = (lo + step).min(support_hi);
            if g(hi)? >= 0.0 {
                break;
            }
            if hi >= support_hi || !hi.is_finite() {
                return Err(unreachable());
            }
            lo = hi;
            step *= 2.0;
        }
    } else {
        for _ in 0..2100 {
            lo = (hi - step).max(support_lo);
            if g(lo)? <= 0.0 {
                break;
            }
            if lo <= support_lo || !lo.is_finite() {
                return Err(unreachable());
            }
            hi = lo;
            step *= 2.0;
        }
    }
    if !(g(lo)? <= 0.0 && g(hi)? >= 0.0) {
        return Err(unreachable());
    }

    // Bisection until the interval cannot shrink in f64.
    for _ in 0..400 {
        let mid = lo + 0.5 * (hi - lo);
        if mid <= lo || mid >= hi {
            break;
        }
        if g(mid)? <= 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

/// Inverse for binomial/Poisson: an integer quantile by binary search (see
/// [`critical_value`] for the convention).
fn solve_discrete(model: &Model, ln_target: f64, upper: bool) -> Result<f64> {
    let tails = |k: u64| model.ln_tails(k as f64);
    // The predicate that is false below the answer and true from it onward.
    // Upper: smallest k with P(X > k) <= p.  Lower: largest k with P(X <= k) <= p,
    // i.e. one below the smallest k with P(X <= k) > p.
    let holds = |k: u64| -> Result<bool> {
        let (lower, upper_tail) = tails(k)?;
        Ok(if upper {
            upper_tail <= ln_target
        } else {
            lower > ln_target
        })
    };
    // Largest candidate: the binomial's n, or a doubling search for Poisson.
    let top = match *model {
        Model::Binomial { n, .. } => n as u64,
        _ => {
            let mut k: u64 = 1;
            while !holds(k)? {
                k = k
                    .checked_mul(2)
                    .filter(|&v| v < (1 << 53))
                    .ok_or_else(|| Error::NoSolution("quantile is beyond 2^53".into()))?;
            }
            k
        }
    };
    if !holds(top)? {
        // Only possible for the binomial Lower tail when p is so large that
        // even P(X <= n) = 1 does not exceed it, i.e. p = 1 (rejected earlier).
        return Err(Error::NoSolution(
            "no integer satisfies the requested probability".into(),
        ));
    }
    // Smallest k in [0, top] with holds(k).
    let (mut lo, mut hi) = (0_u64, top);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if holds(mid)? {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    if upper {
        Ok(lo as f64)
    } else if lo == 0 {
        Err(Error::NoSolution(
            "even P(X <= 0) exceeds the requested lower-tail probability".into(),
        ))
    } else {
        Ok((lo - 1) as f64)
    }
}

// ---------------------------------------------------------------------------
// Convenience methods on StatSpec
// ---------------------------------------------------------------------------

impl StatSpec {
    /// Probability implied by `statistic`; see [`p_value`].
    pub fn p_value(&self, statistic: f64, tail: Tail) -> Result<Probability> {
        p_value(self, statistic, tail)
    }

    /// The statistic with tail probability `p`; see [`critical_value`].
    pub fn critical_value(&self, p: f64, tail: Tail) -> Result<f64> {
        critical_value(self, p, tail)
    }

    /// Whether the stored value is already a probability (`Pval`, `LogPval`,
    /// `Log10Pval`). A UI should not offer a tail choice for these.
    pub fn is_direct_probability(&self) -> bool {
        matches!(
            self.kind,
            StatKind::Pval | StatKind::LogPval | StatKind::Log10Pval
        )
    }

    /// Whether [`Tail::TwoSided`] is defined for this statistic (it has a
    /// centre of symmetry). `false` for invalid parameters.
    pub fn supports_two_sided(&self) -> bool {
        Model::from_spec(self)
            .ok()
            .and_then(|m| m.centre())
            .is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(kind: StatKind, params: &[f64]) -> StatSpec {
        StatSpec::new(kind, params, 0.0)
    }
    fn close(a: f64, b: f64, rel: f64) -> bool {
        (a - b).abs() <= rel * b.abs().max(1e-300)
    }

    #[test]
    fn tail_is_explicit_and_consistent() {
        let t = spec(StatKind::Ttest, &[10.0]);
        let lo = p_value(&t, 1.3, Tail::Lower).unwrap().p();
        let up = p_value(&t, 1.3, Tail::Upper).unwrap().p();
        let two = p_value(&t, 1.3, Tail::TwoSided).unwrap().p();
        assert!(close(lo + up, 1.0, 1e-14));
        assert!(close(two, 2.0 * up, 1e-14));
        // Negative statistics: lower tail is the small one; two-sided matches |t|.
        let two_neg = p_value(&t, -1.3, Tail::TwoSided).unwrap().p();
        assert!(close(two_neg, two, 1e-14));
        assert!(close(
            p_value(&t, -1.3, Tail::Lower).unwrap().p(),
            up,
            1e-14
        ));
        // At the centre everything is one half or one.
        assert!(close(
            p_value(&t, 0.0, Tail::Upper).unwrap().p(),
            0.5,
            1e-15
        ));
        assert!(close(
            p_value(&t, 0.0, Tail::TwoSided).unwrap().p(),
            1.0,
            1e-15
        ));
    }

    #[test]
    fn two_sided_is_refused_where_it_is_not_defined() {
        for (kind, params) in [
            (StatKind::Ftest, &[3.0, 20.0][..]),
            (StatKind::Chisq, &[4.0]),
            (StatKind::Gamma, &[2.0, 1.0]),
            (StatKind::Beta, &[2.0, 3.0]),
            (StatKind::Binom, &[10.0, 0.3]),
            (StatKind::Poisson, &[3.0]),
        ] {
            let s = spec(kind, params);
            assert!(!s.supports_two_sided(), "{kind:?}");
            assert!(p_value(&s, 1.0, Tail::TwoSided).is_err(), "{kind:?}");
            assert!(
                critical_value(&s, 0.05, Tail::TwoSided).is_err(),
                "{kind:?}"
            );
        }
        assert!(spec(StatKind::Zscore, &[]).supports_two_sided());
    }

    #[test]
    fn tiny_probabilities_survive_in_log_space() {
        let z = spec(StatKind::Zscore, &[]);
        let p = p_value(&z, 40.0, Tail::Upper).unwrap();
        assert!(p.underflows());
        assert!(p.ln_p() < -800.0);
        // And invert it back from the log.
        let back = critical_value_ln(&z, p.ln_p(), Tail::Upper).unwrap();
        assert!(close(back, 40.0, 1e-12), "{back}");
        // Very large t: p is astronomically small but still ordered.
        let t = spec(StatKind::Ttest, &[5.0]);
        let a = p_value(&t, 1e10, Tail::Upper).unwrap().ln_p();
        let b = p_value(&t, 1e12, Tail::Upper).unwrap().ln_p();
        assert!(b < a && a < -50.0);
    }

    #[test]
    fn errors_instead_of_sentinels() {
        let t = spec(StatKind::Ttest, &[10.0]);
        assert!(matches!(
            p_value(&t, f64::NAN, Tail::Upper),
            Err(Error::NonFinite { .. })
        ));
        assert!(p_value(&t, f64::INFINITY, Tail::Upper).is_err());
        assert!(p_value(&spec(StatKind::Ttest, &[0.0]), 1.0, Tail::Upper).is_err());
        assert!(p_value(&spec(StatKind::Ttest, &[-3.0]), 1.0, Tail::Upper).is_err());
        assert!(p_value(&spec(StatKind::Ftest, &[3.0, f64::NAN]), 1.0, Tail::Upper).is_err());
        // Wrong parameter count (the field is public, so it can happen).
        let bad = StatSpec {
            kind: StatKind::Ftest,
            params: vec![3.0],
        };
        assert!(p_value(&bad, 1.0, Tail::Upper).is_err());
        // Out-of-domain statistics.
        let c = spec(StatKind::Correl, &[30.0, 1.0, 1.0]);
        assert!(p_value(&c, 1.5, Tail::Upper).is_err());
        assert!(critical_value(&t, 1.5, Tail::Upper).is_err());
        assert!(critical_value(&t, f64::NAN, Tail::Upper).is_err());
        // Probabilities at the endpoints have no finite critical value.
        assert!(critical_value(&t, 0.0, Tail::Upper).is_err());
        assert!(critical_value(&t, 1.0, Tail::Upper).is_err());
        assert_eq!(critical_value(&t, 1.0, Tail::TwoSided).unwrap(), 0.0);
        // Noncentral distributions are computed (never approximated by the central
        // form), so they have no Unsupported path; invalid parameters still error.
        assert!(p_value(&spec(StatKind::TtestNonc, &[10.0, 1.0]), 1.0, Tail::Upper).is_ok());
    }

    #[test]
    fn correlation_forms_agree_with_the_t_transform() {
        // Simple correlation with n samples, 1 fit, 1 nort => dof = n - 2.
        let r = 0.45;
        let n = 25.0;
        let c = spec(StatKind::Correl, &[n, 1.0, 1.0]);
        let dof = n - 2.0;
        let t = r * (dof / (1.0 - r * r)).sqrt();
        let from_t = p_value(&spec(StatKind::Ttest, &[dof]), t, Tail::TwoSided)
            .unwrap()
            .p();
        let from_r = p_value(&c, r, Tail::TwoSided).unwrap().p();
        assert!(close(from_r, from_t, 1e-12), "{from_r} vs {from_t}");
        // NIfTI one-parameter form maps to the same distribution.
        let std = StatSpec::from_intent(2, [dof, 0.0, 0.0], crate::stat::IntentOrigin::Standard)
            .unwrap()
            .unwrap();
        assert!(close(
            p_value(&std, r, Tail::TwoSided).unwrap().p(),
            from_t,
            1e-12
        ));
        // r = +-1 is certain (p = 0), r = 0 is p = 1.
        assert_eq!(p_value(&c, 1.0, Tail::TwoSided).unwrap().p(), 0.0);
        assert!(close(
            p_value(&c, 0.0, Tail::TwoSided).unwrap().p(),
            1.0,
            1e-15
        ));
        // Parameter validation specific to correlation.
        assert!(p_value(&spec(StatKind::Correl, &[2.0, 1.0, 1.0]), 0.5, Tail::Upper).is_err());
        assert!(p_value(&spec(StatKind::Correl, &[30.0, 0.0, 1.0]), 0.5, Tail::Upper).is_err());
    }

    #[test]
    fn multiple_correlation_is_one_sided() {
        let c = spec(StatKind::Correl, &[40.0, 3.0, 1.0]);
        assert!(!c.supports_two_sided());
        let up = p_value(&c, 0.5, Tail::Upper).unwrap().p();
        let lo = p_value(&c, 0.5, Tail::Lower).unwrap().p();
        assert!(close(up + lo, 1.0, 1e-13));
        assert_eq!(p_value(&c, -0.1, Tail::Upper).unwrap().p(), 1.0);
    }

    #[test]
    fn direct_probability_kinds() {
        let p = spec(StatKind::Pval, &[]);
        assert!(close(
            p_value(&p, 0.03, Tail::Upper).unwrap().p(),
            0.03,
            1e-15
        ));
        // The tail argument cannot change a stored probability.
        assert_eq!(
            p_value(&p, 0.03, Tail::Lower).unwrap(),
            p_value(&p, 0.03, Tail::TwoSided).unwrap()
        );
        assert!(p_value(&p, 1.5, Tail::Upper).is_err());
        assert!(p_value(&p, -0.1, Tail::Upper).is_err());
        assert!(p.is_direct_probability());
        // Log forms use |value| and keep tiny probabilities exactly.
        let l = spec(StatKind::LogPval, &[]);
        assert!(close(
            p_value(&l, 3.0, Tail::Upper).unwrap().p(),
            (-3.0_f64).exp(),
            1e-15
        ));
        assert_eq!(
            p_value(&l, -3.0, Tail::Upper).unwrap(),
            p_value(&l, 3.0, Tail::Upper).unwrap()
        );
        assert_eq!(p_value(&l, 1000.0, Tail::Upper).unwrap().ln_p(), -1000.0);
        let l10 = spec(StatKind::Log10Pval, &[]);
        assert!(close(
            p_value(&l10, 5.0, Tail::Upper).unwrap().p(),
            1e-5,
            1e-14
        ));
        assert!(close(
            p_value(&l10, 400.0, Tail::Upper).unwrap().log10_p(),
            -400.0,
            1e-14
        ));
        // Inverses return the stored representation.
        assert!(close(
            critical_value(&p, 0.02, Tail::Upper).unwrap(),
            0.02,
            1e-15
        ));
        assert!(close(
            critical_value(&l, 0.01, Tail::Upper).unwrap(),
            -(0.01_f64).ln(),
            1e-14
        ));
        assert!(close(
            critical_value(&l10, 0.001, Tail::Upper).unwrap(),
            3.0,
            1e-13
        ));
    }

    #[test]
    fn forward_inverse_round_trips_for_continuous_kinds() {
        let cases: Vec<(StatSpec, Vec<Tail>)> = vec![
            (
                spec(StatKind::Zscore, &[]),
                vec![Tail::Lower, Tail::Upper, Tail::TwoSided],
            ),
            (
                spec(StatKind::Ttest, &[7.5]),
                vec![Tail::Lower, Tail::Upper, Tail::TwoSided],
            ),
            (
                spec(StatKind::Correl, &[30.0, 1.0, 1.0]),
                vec![Tail::Lower, Tail::Upper, Tail::TwoSided],
            ),
            (
                spec(StatKind::Ftest, &[3.0, 20.5]),
                vec![Tail::Lower, Tail::Upper],
            ),
            (
                spec(StatKind::Chisq, &[4.0]),
                vec![Tail::Lower, Tail::Upper],
            ),
            (
                spec(StatKind::Beta, &[2.0, 5.0]),
                vec![Tail::Lower, Tail::Upper],
            ),
            (
                spec(StatKind::Gamma, &[3.0, 2.0]),
                vec![Tail::Lower, Tail::Upper],
            ),
            (
                spec(StatKind::Normal, &[1.0, 2.0]),
                vec![Tail::Lower, Tail::Upper, Tail::TwoSided],
            ),
            (
                spec(StatKind::Logistic, &[0.5, 2.0]),
                vec![Tail::Lower, Tail::Upper, Tail::TwoSided],
            ),
            (
                spec(StatKind::Laplace, &[0.5, 2.0]),
                vec![Tail::Lower, Tail::Upper, Tail::TwoSided],
            ),
            (
                spec(StatKind::Uniform, &[-1.0, 3.0]),
                vec![Tail::Lower, Tail::Upper, Tail::TwoSided],
            ),
            (
                spec(StatKind::Weibull, &[0.5, 2.0, 1.5]),
                vec![Tail::Lower, Tail::Upper],
            ),
            (spec(StatKind::Chi, &[3.0]), vec![Tail::Lower, Tail::Upper]),
            (
                spec(StatKind::Extval, &[0.5, 2.0]),
                vec![Tail::Lower, Tail::Upper],
            ),
            (spec(StatKind::Invgauss, &[2.0, 3.0]), vec![Tail::Lower]),
        ];
        for (s, tails) in &cases {
            for &tail in tails {
                for &p in &[1e-10, 1e-4, 0.01, 0.05, 0.3, 0.5, 0.9] {
                    // A uniform's quantile for p = 1e-10 sits within 1e-9 of the
                    // support edge, where f64 cannot resolve p better than ~1e-7;
                    // that is conditioning of the question, not an error.
                    if s.kind == StatKind::Uniform && p < 1e-6 {
                        continue;
                    }
                    let x = critical_value(s, p, tail)
                        .unwrap_or_else(|e| panic!("{:?} {tail:?} p={p}: {e}", s.kind));
                    let back = p_value(s, x, tail).unwrap().p();
                    assert!(
                        close(back, p, 1e-8),
                        "{:?} {tail:?}: p={p} -> x={x} -> {back}",
                        s.kind
                    );
                }
            }
        }
    }

    #[test]
    fn probabilities_are_monotone_in_the_statistic() {
        let specs = [
            spec(StatKind::Ttest, &[12.0]),
            spec(StatKind::Ftest, &[4.0, 30.0]),
            spec(StatKind::Chisq, &[6.0]),
            spec(StatKind::Zscore, &[]),
            spec(StatKind::Gamma, &[2.5, 1.5]),
            spec(StatKind::Beta, &[2.0, 2.0]),
            spec(StatKind::Binom, &[20.0, 0.4]),
            spec(StatKind::Poisson, &[6.0]),
        ];
        for s in &specs {
            let mut prev_up = 2.0;
            let mut prev_lo = -1.0;
            for i in 0..200 {
                let x = -3.0 + 0.1 * f64::from(i);
                let up = p_value(s, x, Tail::Upper).unwrap().p();
                let lo = p_value(s, x, Tail::Lower).unwrap().p();
                assert!(
                    up <= prev_up + 1e-15 && lo >= prev_lo - 1e-15,
                    "{:?} at {x}",
                    s.kind
                );
                assert!((up + lo - 1.0).abs() < 1e-12, "{:?} at {x}", s.kind);
                prev_up = up;
                prev_lo = lo;
            }
        }
    }

    #[test]
    fn discrete_distributions_use_the_step_function_and_integer_quantiles() {
        let b = spec(StatKind::Binom, &[10.0, 0.3]);
        // P(X > 5) for Binomial(10, 0.3) = 0.047348987...
        assert!(close(
            p_value(&b, 5.0, Tail::Upper).unwrap().p(),
            0.047_348_987_4,
            1e-8
        ));
        // A non-integer statistic counts as its floor.
        assert_eq!(
            p_value(&b, 5.9, Tail::Upper).unwrap(),
            p_value(&b, 5.0, Tail::Upper).unwrap()
        );
        // Quantile convention: smallest k with P(X > k) <= p.
        let k = critical_value(&b, 0.05, Tail::Upper).unwrap();
        assert_eq!(k, 5.0);
        assert!(p_value(&b, k, Tail::Upper).unwrap().p() <= 0.05);
        assert!(p_value(&b, k - 1.0, Tail::Upper).unwrap().p() > 0.05);
        // Lower: largest k with P(X <= k) <= p.
        let k = critical_value(&b, 0.05, Tail::Lower).unwrap();
        assert!(p_value(&b, k, Tail::Lower).unwrap().p() <= 0.05);
        assert!(p_value(&b, k + 1.0, Tail::Lower).unwrap().p() > 0.05);
        // Bounds and rejections.
        assert_eq!(p_value(&b, 10.0, Tail::Upper).unwrap().p(), 0.0);
        assert_eq!(p_value(&b, -1.0, Tail::Lower).unwrap().p(), 0.0);
        assert!(p_value(&spec(StatKind::Binom, &[10.5, 0.3]), 1.0, Tail::Upper).is_err());
        assert!(p_value(&spec(StatKind::Binom, &[10.0, 1.3]), 1.0, Tail::Upper).is_err());
        // Poisson: P(X > 3) with mean 3 = 0.352768...
        let po = spec(StatKind::Poisson, &[3.0]);
        assert!(close(
            p_value(&po, 3.0, Tail::Upper).unwrap().p(),
            0.352_768_111,
            1e-7
        ));
        let k = critical_value(&po, 1e-6, Tail::Upper).unwrap();
        assert!(p_value(&po, k, Tail::Upper).unwrap().p() <= 1e-6);
        assert!(p_value(&po, k - 1.0, Tail::Upper).unwrap().p() > 1e-6);
        // A lower-tail probability smaller than P(X = 0) has no integer answer.
        assert!(critical_value(&po, 1e-6, Tail::Lower).is_err());
    }

    #[test]
    fn threshold_comparison_uses_logs() {
        let z = spec(StatKind::Zscore, &[]);
        let strong = p_value(&z, 40.0, Tail::Upper).unwrap(); // underflows as a plain f64
        let weak = p_value(&z, 39.0, Tail::Upper).unwrap();
        assert_eq!(strong.p(), 0.0);
        assert_eq!(weak.p(), 0.0); // both read as zero...
        assert!(strong.is_at_most(weak) && !weak.is_at_most(strong)); // ...but are ordered
        let threshold = Probability::from_p(0.05).unwrap();
        assert!(p_value(&z, 2.0, Tail::Upper).unwrap().is_at_most(threshold));
        assert!(!p_value(&z, 1.0, Tail::Upper).unwrap().is_at_most(threshold));
        // The boundary is inclusive: a sample exactly at the threshold passes.
        assert!(threshold.is_at_most(threshold));
    }

    #[test]
    fn noncentral_chi_and_f_match_independent_summations() {
        // Reference values from a separate pure-Python implementation of the same
        // Poisson mixtures (different code, double precision); AFNI's CDFLIB is
        // only accurate to ~1e-5 here, so it cannot arbitrate these.
        let chi = spec(StatKind::ChisqNonc, &[3.0, 2.0]);
        let p = p_value(&chi, 2.0, Tail::Lower).unwrap().p();
        assert!(close(p, 0.220_733_087_074_121_3, 1e-11), "{p}");
        let f = spec(StatKind::FtestNonc, &[3.0, 20.0, 2.0]);
        let p = p_value(&f, 3.0, Tail::Lower).unwrap().p();
        assert!(close(p, 0.823_755_494_809_820_7, 1e-11), "{p}");
        let f2 = spec(StatKind::FtestNonc, &[5.0, 12.5, 8.0]);
        let p = p_value(&f2, 10.0, Tail::Lower).unwrap().p();
        assert!(close(p, 0.983_563_361_612_018, 1e-11), "{p}");
    }

    #[test]
    fn zero_noncentrality_is_the_central_distribution() {
        let pairs = [
            (
                spec(StatKind::ChisqNonc, &[4.0, 0.0]),
                spec(StatKind::Chisq, &[4.0]),
                3.3,
            ),
            (
                spec(StatKind::FtestNonc, &[3.0, 12.0, 0.0]),
                spec(StatKind::Ftest, &[3.0, 12.0]),
                2.1,
            ),
            (
                spec(StatKind::TtestNonc, &[9.0, 0.0]),
                spec(StatKind::Ttest, &[9.0]),
                1.4,
            ),
            (
                spec(StatKind::TtestNonc, &[9.0, 0.0]),
                spec(StatKind::Ttest, &[9.0]),
                -2.2,
            ),
        ];
        for (nc, central, x) in &pairs {
            for tail in [Tail::Lower, Tail::Upper] {
                let a = p_value(nc, *x, tail).unwrap().p();
                let b = p_value(central, *x, tail).unwrap().p();
                assert!(
                    close(a, b, 1e-12),
                    "{:?} {tail:?} at {x}: {a} vs {b}",
                    nc.kind
                );
            }
        }
        // A central-equivalent noncentral t keeps its two-sided definition.
        assert!(spec(StatKind::TtestNonc, &[9.0, 0.0]).supports_two_sided());
        assert!(!spec(StatKind::TtestNonc, &[9.0, 1.0]).supports_two_sided());
        assert!(p_value(
            &spec(StatKind::FtestNonc, &[3.0, 4.0, 1.0]),
            1.0,
            Tail::TwoSided
        )
        .is_err());
    }

    #[test]
    fn noncentral_parameters_are_validated() {
        for bad in [
            spec(StatKind::ChisqNonc, &[0.0, 1.0]),
            spec(StatKind::ChisqNonc, &[3.0, -1.0]),
            spec(StatKind::FtestNonc, &[3.0, 0.0, 1.0]),
            spec(StatKind::FtestNonc, &[3.0, 4.0, f64::NAN]),
            spec(StatKind::TtestNonc, &[-2.0, 1.0]),
            spec(StatKind::TtestNonc, &[5.0, f64::INFINITY]),
        ] {
            assert!(
                p_value(&bad, 1.0, Tail::Upper).is_err(),
                "{:?} {:?}",
                bad.kind,
                bad.params
            );
        }
    }

    #[test]
    fn noncentral_t_series_and_quadrature_agree() {
        // The AS 243 series and the direct integral are independent
        // computations of the same upper tail; they must agree wherever the
        // series is trustworthy (delta >= 0, and delta < 0 with mild cancellation).
        for &(t, dof, delta) in &[
            (0.5, 5.0, 1.5),
            (2.0, 10.0, 1.5),
            (6.0, 10.0, 1.5),
            (3.0, 30.0, 4.0),
            (1.0, 3.5, 0.0),
            (4.0, 12.0, -0.5),
            (2.0, 10.0, -1.0),
        ] {
            let (_, ln_series) = noncentral_t_tails(t, dof, delta).unwrap();
            let ln_quad = noncentral_t_upper_by_quadrature(t, dof, delta).unwrap();
            assert!(
                (ln_series - ln_quad).abs() < 1e-8,
                "t={t} dof={dof} delta={delta}: {ln_series} vs {ln_quad}"
            );
        }
    }

    #[test]
    fn noncentral_t_far_tail_on_the_wrong_side_stays_accurate() {
        // t = 1 with delta = -4 and 30 dof: P(T > 1) is tiny and the AS 243
        // difference cancels almost completely; the quadrature fallback must
        // keep it ordered, finite and consistent with its complement.
        let s = spec(StatKind::TtestNonc, &[30.0, -4.0]);
        let mut prev = f64::INFINITY;
        for t in [0.5, 1.0, 2.0, 4.0, 8.0, 20.0] {
            let up = p_value(&s, t, Tail::Upper).unwrap();
            let lo = p_value(&s, t, Tail::Lower).unwrap();
            assert!(up.ln_p() < prev, "not decreasing at {t}");
            prev = up.ln_p();
            assert!((up.p() + lo.p() - 1.0).abs() < 1e-9, "t={t}");
        }
        // And the mirror image: a lower tail on the wrong side.
        let mirror = spec(StatKind::TtestNonc, &[30.0, 4.0]);
        let a = p_value(&mirror, -1.0, Tail::Lower).unwrap().ln_p();
        let b = p_value(&s, 1.0, Tail::Upper).unwrap().ln_p();
        assert!((a - b).abs() < 1e-12);
    }

    #[test]
    fn noncentral_inverses_round_trip() {
        for s in [
            spec(StatKind::ChisqNonc, &[6.0, 4.0]),
            spec(StatKind::FtestNonc, &[4.0, 18.0, 5.0]),
            spec(StatKind::TtestNonc, &[15.0, 2.0]),
            spec(StatKind::TtestNonc, &[15.0, -2.0]),
        ] {
            for tail in [Tail::Lower, Tail::Upper] {
                for p in [1e-8, 1e-3, 0.05, 0.5, 0.95] {
                    let x = critical_value(&s, p, tail)
                        .unwrap_or_else(|e| panic!("{:?} {tail:?} p={p}: {e}", s.kind));
                    let back = p_value(&s, x, tail).unwrap().p();
                    assert!(
                        close(back, p, 1e-7),
                        "{:?} {tail:?} p={p}: x={x} -> {back}",
                        s.kind
                    );
                }
            }
        }
    }

    #[test]
    fn gamma_second_parameter_is_a_rate() {
        // Gamma(shape 3, rate 2) at x = 2 equals Gamma(3, 1) at x = 4.
        let a = p_value(&spec(StatKind::Gamma, &[3.0, 2.0]), 2.0, Tail::Upper)
            .unwrap()
            .p();
        let b = p_value(&spec(StatKind::Gamma, &[3.0, 1.0]), 4.0, Tail::Upper)
            .unwrap()
            .p();
        assert!(close(a, b, 1e-14));
    }
}
