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
// The shared test harness for comparing `afni-core` against AFNI. Integration
// tests pull it in with `mod common;`. It provides, in one standard place:
//
//   * loading committed AFNI "reference" fixtures (tests/data/conformance/),
//   * a tolerance-aware comparison (`assert_close`),
//   * an optional LIVE replay of the same cases against real AFNI programs,
//   * small property-test helpers (monotonicity, round trips).
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// Every later phase's conformance test (p-values, FDR curves, color maps,
// thresholds, clusters) should load its expected values through `Fixture` and
// compare with `assert_close`, so there is exactly one way to compare against
// AFNI. This file is test-only; it is not part of the published library.
//
// THE TWO MODES
//
//   cargo test                     Uses only committed fixtures. AFNI is NOT
//                                  needed on PATH. This is what CI runs.
//   AFNI_CORE_LIVE=1 cargo test    Additionally re-runs each fixture case with
//                                  the real AFNI program and checks the
//                                  committed value still matches. Use it when
//                                  regenerating fixtures or upgrading AFNI.
// ---------------------------------------------------------------------------

// Each integration-test file is compiled as its own crate and uses only some of
// these helpers, so unused-function warnings would be noise.
#![allow(dead_code)]

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::Command;

// ---------------------------------------------------------------------------
// Locating and loading fixtures
// ---------------------------------------------------------------------------

/// Path to a committed fixture, relative to `tests/data/`.
///
/// Panics if the file is missing: a committed fixture that has vanished is a
/// repository error, not a reason to silently skip a test.
pub fn data(relative: &str) -> PathBuf {
    // CARGO_MANIFEST_DIR is the afni-core crate directory, set by Cargo at
    // compile time, so this works no matter where `cargo test` is launched.
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(relative);
    assert!(
        path.exists(),
        "missing committed fixture {}",
        path.display()
    );
    path
}

/// One line of a conformance fixture: the AFNI program arguments and the value
/// AFNI printed for them.
#[derive(Debug, Clone, PartialEq)]
pub struct Case {
    /// Arguments exactly as they were passed to the AFNI program,
    /// e.g. `["-t2p", "fitt", "2.0", "10"]`.
    pub args: Vec<String>,
    /// The number AFNI produced.
    pub expected: f64,
    /// 1-based line number in the fixture, for error messages.
    pub line: usize,
}

/// A parsed conformance fixture file.
///
/// File format (see `tests/data/regenerate_conformance.sh`):
///
/// ```text
/// # key: value          <- provenance header, one per line
/// -t2p fitt 2.0 10 => 0.073388
/// ```
#[derive(Debug, Clone)]
pub struct Fixture {
    /// `# key: value` header lines. Must contain `afni_version`.
    pub provenance: BTreeMap<String, String>,
    /// The cases, in file order.
    pub cases: Vec<Case>,
}

impl Fixture {
    /// Load and parse `tests/data/<relative>`.
    ///
    /// Panics with a descriptive message on malformed input; tests should fail
    /// loudly on a corrupt fixture.
    pub fn load(relative: &str) -> Fixture {
        let text = std::fs::read_to_string(data(relative)).expect("read fixture");
        Fixture::parse(&text)
    }

    /// Parse fixture text (separated from `load` so it can be unit-tested).
    pub fn parse(text: &str) -> Fixture {
        let mut provenance = BTreeMap::new();
        let mut cases = Vec::new();

        for (i, raw) in text.lines().enumerate() {
            let line_no = i + 1;
            let line = raw.trim();
            if line.is_empty() {
                continue;
            }
            // Comment lines of the form "# key: value" are provenance;
            // other comments are ignored.
            if let Some(comment) = line.strip_prefix('#') {
                if let Some((key, value)) = comment.split_once(':') {
                    provenance.insert(key.trim().to_owned(), value.trim().to_owned());
                }
                continue;
            }
            // Otherwise: "<args> => <value>".
            let (args, value) = line
                .split_once("=>")
                .unwrap_or_else(|| panic!("line {line_no}: expected '<args> => <value>'"));
            let expected: f64 = value
                .trim()
                .parse()
                .unwrap_or_else(|_| panic!("line {line_no}: bad number {:?}", value.trim()));
            cases.push(Case {
                args: args.split_whitespace().map(str::to_owned).collect(),
                expected,
                line: line_no,
            });
        }
        Fixture { provenance, cases }
    }

    /// The AFNI version string that generated this fixture.
    pub fn afni_version(&self) -> &str {
        self.provenance
            .get("afni_version")
            .expect("fixture is missing the '# afni_version:' provenance line")
    }
}

// ---------------------------------------------------------------------------
// Comparing numbers
// ---------------------------------------------------------------------------

/// How close two numbers must be to count as equal.
///
/// A pair passes if EITHER the absolute difference is within `abs` OR the
/// relative difference is within `rel`. The absolute term makes comparisons
/// near zero sensible; the relative term handles large magnitudes.
#[derive(Debug, Clone, Copy)]
pub struct Tolerance {
    /// Maximum `|actual - expected|`.
    pub abs: f64,
    /// Maximum `|actual - expected| / |expected|`.
    pub rel: f64,
}

impl Tolerance {
    /// Tolerance suited to AFNI programs that print 6 significant digits.
    ///
    /// Rounding to 6 significant digits changes a value by up to 5e-6 of
    /// itself when its leading digit is 1 (e.g. 1.33437e-2 stands for anything
    /// in [1.334365e-2, 1.334375e-2]), hence the relative term.
    pub const SIX_DIGITS: Tolerance = Tolerance {
        abs: 1e-9,
        rel: 5.5e-6,
    };
    /// Tolerance for comparing two f64 computations of the same quantity.
    pub const TIGHT: Tolerance = Tolerance {
        abs: 1e-12,
        rel: 1e-10,
    };
}

/// Whether `actual` is within `tol` of `expected`. NaN is never close to
/// anything, including NaN.
pub fn is_close(actual: f64, expected: f64, tol: Tolerance) -> bool {
    let diff = (actual - expected).abs();
    diff <= tol.abs || diff <= tol.rel * expected.abs()
}

/// Panic with a readable message unless `actual` is close to `expected`.
/// `label` says which case failed (e.g. the fixture line or AFNI arguments).
#[track_caller]
pub fn assert_close(label: &str, actual: f64, expected: f64, tol: Tolerance) {
    assert!(
        is_close(actual, expected, tol),
        "{label}: got {actual:e}, expected {expected:e} (tolerance abs {:e}, rel {:e})",
        tol.abs,
        tol.rel
    );
}

// ---------------------------------------------------------------------------
// Live AFNI replay (opt-in)
// ---------------------------------------------------------------------------

/// True when the user asked for live AFNI comparisons with `AFNI_CORE_LIVE=1`.
pub fn live_enabled() -> bool {
    matches!(std::env::var("AFNI_CORE_LIVE").as_deref(), Ok("1"))
}

/// Run an AFNI program and return its stdout.
///
/// Returns `Err` with a message if the program is not on `PATH` or exits
/// unsuccessfully. Only call this behind [`live_enabled`].
pub fn run_afni(program: &str, args: &[String]) -> Result<String, String> {
    let output = Command::new(program)
        .args(args)
        .output()
        .map_err(|e| format!("could not run `{program}` (is AFNI on PATH?): {e}"))?;
    if !output.status.success() {
        return Err(format!("`{program}` exited with {}", output.status));
    }
    Ok(String::from_utf8_lossy(&output.stdout).into_owned())
}

/// Extract the number from the first line of `cdf` output, which looks like
/// `p = 0.073388` or `t = 2.2414`.
pub fn parse_cdf_output(stdout: &str) -> Option<f64> {
    stdout
        .lines()
        .next()?
        .rsplit('=')
        .next()?
        .trim()
        .parse()
        .ok()
}

// ---------------------------------------------------------------------------
// Property-test helpers
// ---------------------------------------------------------------------------

/// Assert `f` is non-decreasing (or non-increasing when `increasing` is false)
/// over `xs`, which must be sorted ascending. Catches sign errors and
/// discontinuities in tail probabilities and color ramps.
#[track_caller]
pub fn assert_monotonic(name: &str, xs: &[f64], increasing: bool, f: impl Fn(f64) -> f64) {
    let mut prev: Option<(f64, f64)> = None;
    for &x in xs {
        let y = f(x);
        if let Some((px, py)) = prev {
            let ok = if increasing { y >= py } else { y <= py };
            assert!(
                ok,
                "{name} is not monotonic: f({px}) = {py:e} then f({x}) = {y:e}"
            );
        }
        prev = Some((x, y));
    }
}

/// Assert `inverse(forward(x)) == x` (within `tol`) for every `x` in `xs`.
/// This is the standard forward/inverse round-trip check for p-value and
/// critical-value pairs.
#[track_caller]
pub fn assert_round_trip(
    name: &str,
    xs: &[f64],
    tol: Tolerance,
    forward: impl Fn(f64) -> f64,
    inverse: impl Fn(f64) -> f64,
) {
    for &x in xs {
        let back = inverse(forward(x));
        assert_close(&format!("{name} round trip at {x}"), back, x, tol);
    }
}
