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
// Conformance of `afni_core::instacorr` against the steps SUMA's InstaCorr takes,
// executed with AFNI's own libmri functions (tests/data/regenerate_instacorr_refs.sh
// documents exactly which calls, in which order, and why SUMA itself cannot be driven
// from the command line). tests/data/conformance/instacorr.ref has 6 cases of 12 series
// each: the removed-dimension count, the cleaned unit-length rows, and the correlations
// for a single-row seed, an ROI (mean of three rows) seed and an external seed.
//
// Cases WITHOUT polynomial regressors must match to 32-bit accuracy. Cases WITH them
// differ slightly on purpose: AFNI's pseudo-inverse also projects out the rounding-
// noise directions that the filter leaves behind from the constant and linear
// regressors (see instacorr.rs); the difference is measured and bounded here, and the
// removed-dimension count must still match exactly.
// ---------------------------------------------------------------------------

mod common;

use afni_core::instacorr::{prepare_rows, InstaCorrOptions};

struct Case {
    name: String,
    options: InstaCorrOptions,
    rows: Vec<Vec<f64>>,
    seed_rows: Vec<usize>,
    external: Vec<f64>,
    ndof: usize,
    prepared: Vec<Vec<f64>>,
    corr_row: Vec<f64>,
    corr_roi: Vec<f64>,
    corr_ext: Vec<f64>,
}

fn numbers(s: &str) -> Vec<f64> {
    s.split_whitespace().map(|v| v.parse().unwrap()).collect()
}

fn parse() -> Vec<Case> {
    let text = std::fs::read_to_string(common::data("conformance/instacorr.ref")).unwrap();
    let mut out = Vec::new();
    let mut cur: Option<Case> = None;
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let (key, rest) = line.split_once(' ').unwrap_or((line, ""));
        match key {
            "case" => {
                cur = Some(Case {
                    name: rest.to_owned(),
                    options: InstaCorrOptions::default(),
                    rows: vec![],
                    seed_rows: vec![],
                    external: vec![],
                    ndof: 0,
                    prepared: vec![],
                    corr_row: vec![],
                    corr_roi: vec![],
                    corr_ext: vec![],
                })
            }
            "params" => {
                let p: Vec<&str> = rest.split_whitespace().collect();
                let c = cur.as_mut().unwrap();
                let (dt, fbot, ftop): (f64, f64, f64) = (
                    p[2].parse().unwrap(),
                    p[3].parse().unwrap(),
                    p[4].parse().unwrap(),
                );
                c.options = InstaCorrOptions {
                    tr_seconds: Some(dt),
                    band: (p[7] == "1").then_some((fbot, ftop)),
                    polort: p[5].parse().unwrap(),
                    normalize: true,
                    orts: vec![],
                };
            }
            "row" => cur.as_mut().unwrap().rows.push(numbers(rest)),
            "ort" => cur.as_mut().unwrap().options.orts.push(numbers(rest)),
            "seed" => {
                cur.as_mut().unwrap().seed_rows = rest
                    .split_whitespace()
                    .map(|v| v.parse().unwrap())
                    .collect()
            }
            "ext" => cur.as_mut().unwrap().external = numbers(rest),
            "ndof" => cur.as_mut().unwrap().ndof = rest.trim().parse().unwrap(),
            "prep" => cur.as_mut().unwrap().prepared.push(numbers(rest)),
            "corr_row" => cur.as_mut().unwrap().corr_row = numbers(rest),
            "corr_roi" => cur.as_mut().unwrap().corr_roi = numbers(rest),
            "corr_ext" => cur.as_mut().unwrap().corr_ext = numbers(rest),
            "end" => out.push(cur.take().unwrap()),
            other => panic!("unknown key {other}"),
        }
    }
    out
}

fn max_diff(ours: &[f32], theirs: &[f64]) -> f64 {
    ours.iter()
        .zip(theirs)
        .map(|(&a, &b)| (f64::from(a) - b).abs())
        .fold(0.0, f64::max)
}

#[test]
fn matches_suma_steps_run_with_afni_functions() {
    let cases = parse();
    assert_eq!(cases.len(), 6);
    for c in &cases {
        let degenerate = c.options.polort >= 0; // Legendre regressors under a linear detrend
        let p = prepare_rows(&c.rows, &c.options).unwrap();
        // The removed-dimension count is AFNI's, exactly, in every case.
        assert_eq!(p.removed_dof(), c.ndof, "{}: removed dimensions", c.name);

        let seed = c.seed_rows[0];
        let single = p.correlate_row(seed).unwrap();
        let roi = p.correlate_rows(&c.seed_rows).unwrap();
        let ext = p.correlate_external(&c.external).unwrap();
        let diffs = [
            ("single seed", max_diff(&single.values, &c.corr_row)),
            ("ROI seed", max_diff(&roi.values, &c.corr_roi)),
            ("external seed", max_diff(&ext.values, &c.corr_ext)),
        ];
        // Without polynomial regressors: agreement to 32-bit accuracy. With them: the
        // effect of AFNI also removing two rounding-noise directions.
        let limit = if degenerate { 0.06 } else { 2e-5 };
        for (what, d) in diffs {
            eprintln!("{} {what}: max correlation difference {d:.2e}", c.name);
            assert!(d < limit, "{} {what}: {d} (limit {limit})", c.name);
        }
        if !degenerate {
            // The cleaned rows themselves agree too.
            for (r, want) in c.prepared.iter().enumerate() {
                let got = p.row(r).unwrap();
                let worst = got
                    .iter()
                    .zip(want)
                    .map(|(&a, &b)| (f64::from(a) - b).abs())
                    .fold(0.0, f64::max);
                assert!(worst < 2e-5, "{} row {r}: {worst}", c.name);
            }
        }
        // The statistic carries Correl(samples, 1, removed dimensions).
        let stat = single.stat.as_ref().unwrap();
        assert_eq!(
            stat.params,
            vec![c.rows[0].len() as f64, 1.0, c.ndof as f64],
            "{}",
            c.name
        );
    }
}

#[test]
fn the_cases_include_both_kinds_and_real_structure() {
    let cases = parse();
    assert!(cases
        .iter()
        .any(|c| c.options.polort < 0 && c.options.band.is_some()));
    assert!(cases.iter().any(|c| c.options.polort >= 0));
    assert!(cases.iter().any(|c| !c.options.orts.is_empty()));
    assert!(cases.iter().any(|c| c.options.band.is_none()));
    // The reference correlations are not all near zero (so a constant could not pass).
    let c = cases
        .iter()
        .find(|c| c.name == "bandpass_no_polort")
        .unwrap();
    assert!(c.corr_row.iter().filter(|v| v.abs() > 0.5).count() >= 3);
    // An odd-length series is included (FFT padding).
    assert!(cases.iter().any(|c| c.rows[0].len() % 2 == 1));
}

#[test]
fn live_harness_agrees_with_the_committed_cases() {
    if !common::live_enabled() {
        eprintln!("skipping: set AFNI_CORE_LIVE=1 to rebuild the harness and re-run AFNI");
        return;
    }
    let committed = std::fs::read_to_string(common::data("conformance/instacorr.ref")).unwrap();
    let target = std::env::temp_dir().join(format!("afni_core_live_{}.ref", std::process::id()));
    let out = std::process::Command::new("bash")
        .arg(common::data("regenerate_instacorr_refs.sh"))
        .env("REGEN_OUT", &target)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let regenerated = std::fs::read_to_string(&target).unwrap();
    let _ = std::fs::remove_file(&target);
    let body = |s: &str| {
        s.lines()
            .filter(|l| !l.starts_with("# afni_version"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(body(&regenerated), body(&committed));
}
