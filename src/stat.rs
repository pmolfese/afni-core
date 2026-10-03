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
// `StatKind` (which statistical distribution a value follows) and `StatSpec`
// (a kind plus its parameters, such as `Ttest(23)`).
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * Originally written in `afni-io/src/stat.rs`. Under the type-ownership
//   decision in docs/ARCHITECTURE.md it moved here, and `afni-io` re-exports
//   it from the old path (`afni_io::stat::StatSpec`), so existing code keeps
//   compiling. Moving it avoids a second, drifting copy.
// * `dataset.rs` / `column.rs` attach a `StatSpec` to each data column.
// * Phase 2 (statistics) will add the math: p-values and critical values for
//   each kind. This file deliberately has NONE of that; it only says WHAT a
//   statistic is. In particular it does not decide one- vs two-sided tails.
//
// A note on the "statsym" text form (`Ttest(23)`): it is AFNI's canonical
// *name* for a statistic, used identically in .HEAD files, NIML datasets and
// NIfTI extensions. It is a naming convention rather than the syntax of any one
// file format, so `parse` / `to_statsym` live here. Splitting a whole
// `;`-separated *list* of them across columns is file layout and stays in
// `afni-io` (`afni_io::stat::parse_statsym_list`).
// ---------------------------------------------------------------------------

//! Statistic kinds and specifications (`StatKind`, `StatSpec`).
//!
//! AFNI's statistic codes are the NIfTI intent codes 2-24
//! (`NI_STAT_FIRSTCODE` to `NI_STAT_LASTCODE` in `niml.h`), so one [`StatSpec`]
//! describes a statistic wherever it came from.
//!
//! References: `afni/src/niml/niml_stat.c` (`NI_stat_decode`,
//! `NI_stat_encode`), `3ddata.h` (`FUNC_*_TYPE`).

use crate::error::{Error, Result};

/// A statistic's distribution: AFNI stat code = NIfTI intent code.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum StatKind {
    /// 2: correlation coefficient. Params: samples, fit parameters, ort parameters.
    Correl = 2,
    /// 3: Student t. Params: degrees of freedom.
    Ttest = 3,
    /// 4: F. Params: numerator and denominator degrees of freedom.
    Ftest = 4,
    /// 5: standard normal z.
    Zscore = 5,
    /// 6: chi-squared. Params: degrees of freedom.
    Chisq = 6,
    /// 7: beta. Params: a, b.
    Beta = 7,
    /// 8: binomial. Params: trials, probability per trial.
    Binom = 8,
    /// 9: gamma. Params: shape, scale.
    Gamma = 9,
    /// 10: Poisson. Params: mean.
    Poisson = 10,
    /// 11: normal. Params: mean, standard deviation.
    Normal = 11,
    /// 12: noncentral F. Params: numerator dof, denominator dof, noncentrality.
    FtestNonc = 12,
    /// 13: noncentral chi-squared. Params: dof, noncentrality.
    ChisqNonc = 13,
    /// 14: logistic. Params: location, scale.
    Logistic = 14,
    /// 15: Laplace. Params: location, scale.
    Laplace = 15,
    /// 16: uniform. Params: start, end.
    Uniform = 16,
    /// 17: noncentral t. Params: dof, noncentrality.
    TtestNonc = 17,
    /// 18: Weibull. Params: location, scale, power.
    Weibull = 18,
    /// 19: chi. Params: dof.
    Chi = 19,
    /// 20: inverse Gaussian. Params: mu, lambda.
    Invgauss = 20,
    /// 21: extreme value type I. Params: location, scale.
    Extval = 21,
    /// 22: the value is a p-value.
    Pval = 22,
    /// 23: the value encodes a probability as `-ln(p)`; `p = exp(-|value|)`.
    /// (NIfTI permits a signed log value but its reference library emits the
    /// positive `-ln(p)`, so the sign is ignored.)
    LogPval = 23,
    /// 24: the value encodes a probability as `-log10(p)`; `p = 10^(-|value|)`.
    Log10Pval = 24,
}

const ALL: [StatKind; 23] = [
    StatKind::Correl,
    StatKind::Ttest,
    StatKind::Ftest,
    StatKind::Zscore,
    StatKind::Chisq,
    StatKind::Beta,
    StatKind::Binom,
    StatKind::Gamma,
    StatKind::Poisson,
    StatKind::Normal,
    StatKind::FtestNonc,
    StatKind::ChisqNonc,
    StatKind::Logistic,
    StatKind::Laplace,
    StatKind::Uniform,
    StatKind::TtestNonc,
    StatKind::Weibull,
    StatKind::Chi,
    StatKind::Invgauss,
    StatKind::Extval,
    StatKind::Pval,
    StatKind::LogPval,
    StatKind::Log10Pval,
];

impl StatKind {
    /// Map a stat code / NIfTI intent code (2–24) onto a variant.
    pub fn from_code(code: i64) -> Option<Self> {
        ALL.iter().copied().find(|k| k.code() == code)
    }

    /// The stat code, which is also the NIfTI intent code.
    pub fn code(self) -> i64 {
        self as i64
    }

    /// The name used in `BRICK_STATSYM` / `COLMS_STATSYM` strings, e.g.
    /// `Ttest` (`distname` in `niml_stat.c`).
    pub fn name(self) -> &'static str {
        match self {
            Self::Correl => "Correl",
            Self::Ttest => "Ttest",
            Self::Ftest => "Ftest",
            Self::Zscore => "Zscore",
            Self::Chisq => "Chisq",
            Self::Beta => "Beta",
            Self::Binom => "Binom",
            Self::Gamma => "Gamma",
            Self::Poisson => "Poisson",
            Self::Normal => "Normal",
            Self::FtestNonc => "Ftest_nonc",
            Self::ChisqNonc => "Chisq_nonc",
            Self::Logistic => "Logistic",
            Self::Laplace => "Laplace",
            Self::Uniform => "Uniform",
            Self::TtestNonc => "Ttest_nonc",
            Self::Weibull => "Weibull",
            Self::Chi => "Chi",
            Self::Invgauss => "Invgauss",
            Self::Extval => "Extval",
            Self::Pval => "Pval",
            Self::LogPval => "LogPval",
            Self::Log10Pval => "Log10Pval",
        }
    }

    /// Number of distribution parameters (`numparam` in `niml_stat.c`).
    pub fn num_params(self) -> usize {
        match self {
            Self::Zscore | Self::Pval | Self::LogPval | Self::Log10Pval => 0,
            Self::Ttest | Self::Chisq | Self::Poisson | Self::Chi => 1,
            Self::Correl | Self::FtestNonc | Self::Weibull => 3,
            _ => 2,
        }
    }

    /// Whether AFNI's `.HEAD` loader keeps this code. It stores only the
    /// classic `FUNC_*_TYPE` codes 2–10 (`FUNC_IS_STAT`, `3ddata.h`) and
    /// drops the rest; this crate keeps all of them.
    pub fn is_afni_brick_stat(self) -> bool {
        self.code() <= 10
    }
}

/// A statistic and its parameters, e.g. `Ttest(23)`.
#[derive(Debug, Clone, PartialEq)]
pub struct StatSpec {
    /// The distribution.
    pub kind: StatKind,
    /// Exactly [`StatKind::num_params`] parameters.
    pub params: Vec<f64>,
}

impl StatSpec {
    /// Build a spec from however many parameters are available: extras are
    /// dropped and missing ones are set to `fill`. (AFNI fills with 1.0 when
    /// decoding a STATSYM string and with 0.0 for `BRICK_STATAUX`.)
    pub fn new(kind: StatKind, params: &[f64], fill: f64) -> Self {
        let n = kind.num_params();
        let mut params: Vec<f64> = params.iter().copied().take(n).collect();
        params.resize(n, fill);
        Self { kind, params }
    }

    /// Parse one STATSYM entry such as `Ttest(23)` or `Ftest(2,40)`, as
    /// `NI_stat_decode` does: the name is matched case-insensitively, and a
    /// missing or unreadable parameter becomes 1.0. `none`, an empty string or
    /// an unknown name gives `None`.
    pub fn parse(text: &str) -> Option<Self> {
        let text = text.trim();
        let bytes = text.as_bytes();
        let kind = ALL.iter().copied().find(|k| {
            let n = k.name().len();
            bytes.len() > n
                && bytes[..n].eq_ignore_ascii_case(k.name().as_bytes())
                && bytes[n] == b'('
        })?;
        let args = &text[kind.name().len() + 1..];
        let values: Vec<f64> = args
            .split([',', ')'])
            .take(kind.num_params())
            .map(|v| v.trim().parse().unwrap_or(1.0))
            .collect();
        Some(Self::new(kind, &values, 1.0))
    }

    /// Format as AFNI writes it (`NI_stat_encode`), e.g. `Ftest(2,40)` or
    /// `Zscore()`.
    pub fn to_statsym(&self) -> String {
        let params: Vec<String> = self
            .params
            .iter()
            .map(|&v| format_param(v as f32))
            .collect();
        format!("{}({})", self.kind.name(), params.join(","))
    }

    /// The statistic described by a NIfTI intent code and its three
    /// `intent_p1..3` values, interpreted according to who wrote them.
    ///
    /// Returns `Ok(None)` when `code` is not a statistic, and an error when the
    /// parameters cannot be read under the stated (or detected) origin.
    ///
    /// Only correlation (code 2) is parameterized differently by the two
    /// conventions, which is why [`IntentOrigin`] exists:
    ///
    /// * AFNI: `Correl(samples, nfit, nort)`, three parameters.
    /// * NIfTI standard: one parameter, the degrees of freedom `d`, with the
    ///   other two unused (zero). That is the same distribution as AFNI's
    ///   `Correl(d + 1, 1, 0)`, and this function returns it in that form, so
    ///   downstream code only ever sees one parameterization.
    ///
    /// With [`IntentOrigin::Unknown`], correlation is classified by structure
    /// rather than guessed: AFNI always has `nfit >= 1`, while the standard
    /// leaves the second and third parameters zero. Anything else is an error.
    pub fn from_intent(code: i64, params: [f64; 3], origin: IntentOrigin) -> Result<Option<Self>> {
        let Some(kind) = StatKind::from_code(code) else {
            return Ok(None);
        };
        if kind != StatKind::Correl {
            // Every other statistic means the same thing in both conventions.
            return Ok(Some(Self::new(kind, &params, 0.0)));
        }
        let origin = match origin {
            IntentOrigin::Unknown => {
                if params[1] >= 1.0 {
                    IntentOrigin::Afni
                } else if params[1] == 0.0 && params[2] == 0.0 {
                    IntentOrigin::Standard
                } else {
                    return Err(Error::InvalidParameter {
                        name: "correlation intent parameters".into(),
                        reason: format!(
                            "{params:?} fit neither AFNI's (samples, nfit >= 1, nort) nor the \
                             NIfTI standard's (dof, 0, 0); say which convention applies"
                        ),
                    });
                }
            }
            known => known,
        };
        Ok(Some(match origin {
            IntentOrigin::Standard => {
                if params[1] != 0.0 || params[2] != 0.0 {
                    return Err(Error::InvalidParameter {
                        name: "correlation intent parameters".into(),
                        reason: "the NIfTI standard defines one parameter (dof); \
                                 intent_p2 and intent_p3 must be zero"
                            .into(),
                    });
                }
                Self::new(kind, &[params[0] + 1.0, 1.0, 0.0], 0.0)
            }
            _ => Self::new(kind, &params, 0.0),
        }))
    }

    /// Like [`from_intent`](Self::from_intent) with [`IntentOrigin::Unknown`],
    /// but reports an ambiguous or malformed correlation as `None`, the same as
    /// "not a statistic". Kept for callers that cannot surface an error; prefer
    /// `from_intent`.
    pub fn from_nifti_intent(code: i64, params: [f64; 3]) -> Option<Self> {
        Self::from_intent(code, params, IntentOrigin::Unknown)
            .ok()
            .flatten()
    }
}

/// Who wrote a NIfTI/GIfTI `intent_p1..3` triple. Matters only for correlation;
/// see [`StatSpec::from_intent`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum IntentOrigin {
    /// Written by AFNI, which copies its own stat parameters verbatim.
    Afni,
    /// Written by a standards-following tool.
    Standard,
    /// Not known; classify by structure and fail if ambiguous. The default.
    #[default]
    Unknown,
}

/// `NI_fval_to_char` from `niml_stat.c`: integers print as integers, other
/// values with a precision chosen by magnitude and trailing zeros removed.
/// (AFNI pads its `%-9.0f` case with trailing blanks; they are trimmed here.)
fn format_param(v: f32) -> String {
    if v == 0.0 {
        return "0".into();
    }
    if v.abs() < 99_999_999.0 && v == v.trunc() {
        return format!("{}", v as i64);
    }
    let strip = |s: String| {
        let trimmed = s.trim_end_matches('0');
        if trimmed.len() > 1 {
            trimmed.to_string()
        } else {
            s
        }
    };
    let lv = (10.0001 + f64::from(v.abs()).log10()) as i64;
    let v64 = f64::from(v);
    match lv {
        6..=10 => strip(format!("{v64:.6}")),
        11 => strip(format!("{v64:.5}")),
        12 => strip(format!("{v64:.4}")),
        13 => strip(format!("{v64:.3}")),
        14 => strip(format!("{v64:.2}")),
        15 => strip(format!("{v64:.1}")),
        16 => format!("{v64:.0}"),
        _ => c_exponent(v64, if v > 0.0 { 6 } else { 5 }),
    }
}

/// C's `%.Ne`: mantissa with `digits` decimals and a signed exponent of at
/// least two digits, e.g. `1.500000e-05`.
fn c_exponent(v: f64, digits: usize) -> String {
    let s = format!("{v:.digits$e}");
    let (mantissa, exp) = s.split_once('e').unwrap_or((&s, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    let sign = if exp < 0 { '-' } else { '+' };
    format!("{mantissa}e{sign}{:02}", exp.abs())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_names_and_parameter_counts_match_niml_stat() {
        // numparam[] from niml_stat.c, codes 2..=24.
        let numparam = [
            3, 1, 2, 0, 1, 2, 2, 2, 1, 2, 3, 2, 2, 2, 2, 2, 3, 1, 2, 2, 0, 0, 0,
        ];
        for (i, kind) in ALL.iter().enumerate() {
            assert_eq!(kind.code(), i as i64 + 2);
            assert_eq!(kind.num_params(), numparam[i], "{kind:?}");
            assert_eq!(StatKind::from_code(kind.code()), Some(*kind));
        }
        assert_eq!(StatKind::from_code(1), None);
        assert_eq!(StatKind::from_code(25), None);
        assert!(StatKind::Poisson.is_afni_brick_stat() && !StatKind::Normal.is_afni_brick_stat());
    }

    #[test]
    fn parses_like_ni_stat_decode() {
        let t = StatSpec::parse("Ttest(23)").unwrap();
        assert_eq!(
            (t.kind, t.params.as_slice()),
            (StatKind::Ttest, &[23.0][..])
        );
        let f = StatSpec::parse(" ftest(2, 40) ").unwrap();
        assert_eq!(
            (f.kind, f.params.as_slice()),
            (StatKind::Ftest, &[2.0, 40.0][..])
        );
        // Ttest_nonc is not mistaken for Ttest; missing parameters become 1.0.
        let nc = StatSpec::parse("Ttest_nonc(12)").unwrap();
        assert_eq!(
            (nc.kind, nc.params.as_slice()),
            (StatKind::TtestNonc, &[12.0, 1.0][..])
        );
        assert_eq!(
            StatSpec::parse("Zscore()").unwrap().params,
            Vec::<f64>::new()
        );
        assert_eq!(StatSpec::parse("none"), None);
        assert_eq!(StatSpec::parse("Ttest"), None);
        assert_eq!(StatSpec::parse("Bogus(1)"), None);
    }

    #[test]
    fn encodes_like_ni_stat_encode() {
        let enc = |kind, params: &[f64]| StatSpec::new(kind, params, 0.0).to_statsym();
        assert_eq!(enc(StatKind::Ttest, &[23.0]), "Ttest(23)");
        assert_eq!(enc(StatKind::Ftest, &[2.0, 40.0]), "Ftest(2,40)");
        assert_eq!(enc(StatKind::Zscore, &[]), "Zscore()");
        assert_eq!(enc(StatKind::Correl, &[30.0, 2.0, 1.0]), "Correl(30,2,1)");
        assert_eq!(enc(StatKind::Binom, &[10.0, 0.25]), "Binom(10,0.25)");
        assert_eq!(enc(StatKind::Ttest, &[12.5]), "Ttest(12.5)");
        assert_eq!(
            enc(StatKind::Normal, &[0.00001, -2.5e-7]),
            "Normal(1.000000e-05,-2.50000e-07)"
        );
        for text in ["Ttest(23)", "Ftest(2,40)", "Binom(10,0.25)", "Zscore()"] {
            assert_eq!(StatSpec::parse(text).unwrap().to_statsym(), text);
        }
    }

    #[test]
    fn correlation_intent_parameters_are_not_confused() {
        // AFNI-written: (samples, nfit, nort) kept as is.
        let afni = StatSpec::from_intent(2, [30.0, 2.0, 1.0], IntentOrigin::Afni)
            .unwrap()
            .unwrap();
        assert_eq!(afni.params, [30.0, 2.0, 1.0]);
        // Standard: one dof parameter becomes the equivalent AFNI form.
        let std = StatSpec::from_intent(2, [18.0, 0.0, 0.0], IntentOrigin::Standard)
            .unwrap()
            .unwrap();
        assert_eq!(std.params, [19.0, 1.0, 0.0]);
        // Standard with stray extra parameters is an error, not a silent reinterpretation.
        assert!(StatSpec::from_intent(2, [18.0, 2.0, 0.0], IntentOrigin::Standard).is_err());
        // Unknown origin: classified by structure...
        assert_eq!(
            StatSpec::from_intent(2, [30.0, 2.0, 1.0], IntentOrigin::Unknown)
                .unwrap()
                .unwrap()
                .params,
            [30.0, 2.0, 1.0]
        );
        assert_eq!(
            StatSpec::from_intent(2, [18.0, 0.0, 0.0], IntentOrigin::Unknown)
                .unwrap()
                .unwrap()
                .params,
            [19.0, 1.0, 0.0]
        );
        // ...and an impossible mixture is rejected.
        assert!(StatSpec::from_intent(2, [18.0, 0.5, 0.0], IntentOrigin::Unknown).is_err());
        assert!(StatSpec::from_nifti_intent(2, [18.0, 0.5, 0.0]).is_none());
        // Other statistics do not depend on origin; non-statistics are Ok(None).
        let t = StatSpec::from_intent(3, [10.0, 0.0, 0.0], IntentOrigin::Standard)
            .unwrap()
            .unwrap();
        assert_eq!(
            (t.kind, t.params.as_slice()),
            (StatKind::Ttest, &[10.0][..])
        );
        assert_eq!(
            StatSpec::from_intent(1008, [0.0; 3], IntentOrigin::Afni).unwrap(),
            None
        );
    }
}
