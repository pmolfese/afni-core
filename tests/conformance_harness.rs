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
// Tests of the conformance *harness* itself (tests/common/mod.rs) and of the
// committed AFNI `cdf` fixture. There is no statistics code yet (Phase 2), so
// these tests prove three things Phase 2 will rely on:
//
//   1. fixtures parse, and always record which AFNI version produced them;
//   2. the fixture values are internally sensible (known identities such as
//      "two-sided p at z = 0 is 1");
//   3. when AFNI_CORE_LIVE=1, the committed values still match real AFNI.
//
// When Phase 2 lands, its tests will iterate over the same `Fixture` and call
// the Rust implementation instead of the identity checks below.
// ---------------------------------------------------------------------------

mod common;

use common::{assert_close, assert_monotonic, assert_round_trip, Fixture, Tolerance};

const CDF_FIXTURE: &str = "conformance/cdf.ref";

#[test]
fn fixture_records_afni_provenance() {
    let fx = Fixture::load(CDF_FIXTURE);
    // The version string must name a real AFNI release so a reader can tell
    // which behavior the numbers encode.
    assert!(
        fx.afni_version().contains("AFNI_"),
        "unexpected afni_version: {}",
        fx.afni_version()
    );
    assert!(fx.provenance.contains_key("generated"));
    assert!(!fx.cases.is_empty());
}

#[test]
fn parser_handles_comments_blanks_and_provenance() {
    let fx =
        Fixture::parse("# afni_version: AFNI_TEST\n\n# just a note\n-t2p fitt 2.0 10 => 0.5\n");
    assert_eq!(fx.afni_version(), "AFNI_TEST");
    assert_eq!(fx.cases.len(), 1);
    assert_eq!(fx.cases[0].args, ["-t2p", "fitt", "2.0", "10"]);
    assert_eq!(fx.cases[0].expected, 0.5);
    assert_eq!(fx.cases[0].line, 4);
}

#[test]
#[should_panic(expected = "expected '<args> => <value>'")]
fn parser_rejects_malformed_line() {
    Fixture::parse("-t2p fizt 1.0 0.3\n");
}

#[test]
fn tolerance_semantics() {
    let tol = Tolerance::SIX_DIGITS;
    assert!(common::is_close(0.073388, 0.0733880004, tol));
    // Near zero the absolute term applies...
    assert!(common::is_close(1e-12, 0.0, tol));
    // ...and a genuine disagreement fails.
    assert!(!common::is_close(0.0734, 0.073388, tol));
    // NaN never matches, not even itself.
    assert!(!common::is_close(f64::NAN, f64::NAN, tol));
}

#[test]
fn cdf_fixture_satisfies_known_identities() {
    let fx = Fixture::load(CDF_FIXTURE);
    for case in &fx.cases {
        let label = format!("cdf.ref line {} ({})", case.line, case.args.join(" "));
        let a: Vec<&str> = case.args.iter().map(String::as_str).collect();
        match a.as_slice() {
            // Two-sided p for a symmetric statistic is exactly 1 at zero.
            ["-t2p", "fizt" | "fitt", "0.0", ..] => {
                assert_close(&label, case.expected, 1.0, Tolerance::SIX_DIGITS)
            }
            // Every probability lies in [0, 1].
            ["-t2p", ..] => assert!(
                (0.0..=1.0).contains(&case.expected),
                "{label}: p out of range"
            ),
            // Critical values for p < 1 are positive.
            ["-p2t", ..] => assert!(case.expected > 0.0, "{label}: critical value <= 0"),
            other => panic!("{label}: unrecognised case {other:?}"),
        }
    }
}

#[test]
fn cdf_fixture_probabilities_decrease_with_the_statistic() {
    // Property: larger |t| means a smaller tail probability. Checked on the
    // committed z-test points, sorted ascending.
    let fx = Fixture::load(CDF_FIXTURE);
    let mut pts: Vec<(f64, f64)> = fx
        .cases
        .iter()
        .filter(|c| c.args.starts_with(&["-t2p".into(), "fizt".into()]))
        .map(|c| (c.args[2].parse().unwrap(), c.expected))
        .collect();
    pts.sort_by(|a, b| a.0.total_cmp(&b.0));
    let xs: Vec<f64> = pts.iter().map(|p| p.0).collect();
    assert_monotonic("two-sided z p-value", &xs, false, |x| {
        pts.iter().find(|p| p.0 == x).unwrap().1
    });
}

#[test]
fn round_trip_helper_detects_both_success_and_failure() {
    let xs = [0.5, 1.0, 2.0];
    // A perfect inverse passes...
    assert_round_trip(
        "double/halve",
        &xs,
        Tolerance::TIGHT,
        |x| 2.0 * x,
        |y| y / 2.0,
    );
    // ...and a broken one must be caught.
    let broken = std::panic::catch_unwind(|| {
        assert_round_trip("broken", &[1.0], Tolerance::TIGHT, |x| 2.0 * x, |y| y / 3.0)
    });
    assert!(broken.is_err());
}

/// Opt-in: replay every committed case through the real `cdf` program.
/// Skipped (and says so) unless `AFNI_CORE_LIVE=1`.
#[test]
fn live_afni_matches_committed_cdf_fixture() {
    if !common::live_enabled() {
        eprintln!("skipping: set AFNI_CORE_LIVE=1 to compare against live AFNI");
        return;
    }
    let fx = Fixture::load(CDF_FIXTURE);
    for case in &fx.cases {
        let out = common::run_afni("cdf", &case.args).expect("run cdf");
        let live = common::parse_cdf_output(&out).expect("parse cdf output");
        assert_close(
            &format!("live cdf {}", case.args.join(" ")),
            live,
            case.expected,
            Tolerance::SIX_DIGITS,
        );
    }
}
