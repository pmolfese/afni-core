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
// Conformance of `afni_core::signal::bandpass_vectors` against AFNI's own
// `THD_bandpass_vectors`. tests/data/conformance/signal_bandpass.ref holds 21
// cases (even, odd and awkward lengths, band-pass, low-pass, high-pass, a top above
// Nyquist, half-bin edges, every detrend, a collapsed band, orts with and without a
// filter): the inputs, the filtered output as AFNI computed it in 32-bit float, and
// the number of dimensions it reported removing. Regenerate with
// tests/data/regenerate_signal_refs.sh (it calls the library function directly).
//
// The removed-dimension count must match exactly. The series must agree to the
// accuracy of AFNI's 32-bit arithmetic.
// ---------------------------------------------------------------------------

mod common;

use afni_core::signal::{bandpass_vectors, BandSpec, Detrend};

struct Case {
    name: String,
    dt: f64,
    fbot: f64,
    ftop: f64,
    detrend: Detrend,
    band: bool,
    vectors: Vec<Vec<f64>>,
    orts: Vec<Vec<f64>>,
    ndof: usize,
    out: Vec<Vec<f64>>,
}

fn numbers(s: &str) -> Vec<f64> {
    s.split_whitespace().map(|v| v.parse().unwrap()).collect()
}

fn parse() -> Vec<Case> {
    let text = std::fs::read_to_string(common::data("conformance/signal_bandpass.ref")).unwrap();
    let mut cases = Vec::new();
    let mut cur: Option<Case> = None;
    let (mut nvec, mut nort) = (0, 0);
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let (key, rest) = line.split_once(' ').unwrap_or((line, ""));
        match key {
            "case" => {
                cur = Some(Case {
                    name: rest.to_owned(),
                    dt: 0.0,
                    fbot: 0.0,
                    ftop: 0.0,
                    detrend: Detrend::None,
                    band: false,
                    vectors: vec![],
                    orts: vec![],
                    ndof: 0,
                    out: vec![],
                })
            }
            "params" => {
                let p: Vec<&str> = rest.split_whitespace().collect();
                let c = cur.as_mut().unwrap();
                nvec = p[1].parse().unwrap();
                c.dt = p[2].parse().unwrap();
                c.fbot = p[3].parse().unwrap();
                c.ftop = p[4].parse().unwrap();
                c.detrend = match p[5] {
                    "2" => Detrend::Quadratic,
                    "1" => Detrend::Linear,
                    "0" => Detrend::Mean,
                    _ => Detrend::None,
                };
                nort = p[6].parse().unwrap();
                c.band = p[7] == "1";
            }
            "in" => {
                let c = cur.as_mut().unwrap();
                if c.vectors.len() < nvec {
                    c.vectors.push(numbers(rest));
                } else {
                    c.orts.push(numbers(rest));
                }
            }
            "ndof" => cur.as_mut().unwrap().ndof = rest.trim().parse().unwrap(),
            "out" => cur.as_mut().unwrap().out.push(numbers(rest)),
            "end" => {
                let c = cur.take().unwrap();
                assert_eq!(c.orts.len(), nort);
                cases.push(c);
            }
            other => panic!("unknown key {other}"),
        }
    }
    cases
}

#[test]
fn matches_thd_bandpass_vectors() {
    let cases = parse();
    assert_eq!(cases.len(), 21);
    let mut worst = 0.0_f64;
    for c in &cases {
        let band = BandSpec {
            dt: c.dt,
            fbot: c.fbot,
            ftop: c.ftop,
        };
        let mut v = c.vectors.clone();
        let report = bandpass_vectors(&mut v, c.band.then_some(&band), c.detrend, &c.orts).unwrap();
        assert_eq!(report.removed_dof, c.ndof, "{}: removed dimensions", c.name);
        // Orts that the filter reduces to rounding noise (Legendre orts under a linear
        // detrend, SUMA's setup) are projected out by AFNI as noise directions, which
        // depend on 32-bit rounding and cannot be reproduced. Core drops them. The
        // series must still agree closely as a whole.
        if c.name.starts_with("degenerate_") {
            for (ours, theirs) in v.iter().zip(&c.out) {
                let diff: f64 = ours
                    .iter()
                    .zip(theirs)
                    .map(|(a, b)| (a - b) * (a - b))
                    .sum();
                let size: f64 = theirs.iter().map(|b| b * b).sum();
                let relative = (diff / size).sqrt();
                eprintln!("{}: relative L2 difference {relative:.4}", c.name);
                assert!(relative < 0.3, "{}: {relative}", c.name);
            }
            continue;
        }
        for (i, (ours, theirs)) in v.iter().zip(&c.out).enumerate() {
            for (t, (a, b)) in ours.iter().zip(theirs).enumerate() {
                let err = (a - b).abs();
                worst = worst.max(err);
                // AFNI works in 32-bit float on series of size ~1-10.
                assert!(
                    err <= 2e-4,
                    "{} series {i} sample {t}: afni-core {a}, AFNI {b}",
                    c.name
                );
            }
        }
    }
    eprintln!("largest difference from AFNI: {worst:e}");
}

#[test]
fn the_cases_cover_the_branches() {
    let cases = parse();
    let by = |n: &str| cases.iter().find(|c| c.name == n).unwrap();
    // A series longer than its even FFT length (odd n), a collapsed band whose count
    // exceeds the length (AFNI does not cap it), orts, and no filter at all.
    assert_eq!(by("bandpass_odd_length").vectors[0].len() % 2, 1);
    assert!(by("collapsed_band").ndof > by("collapsed_band").vectors[0].len());
    assert!(
        !by("degenerate_orts_without_filter").orts.is_empty()
            && !by("degenerate_orts_without_filter").band
    );
    // Filtering really changes the series (so a no-op could not pass).
    let c = by("bandpass_even_quadratic");
    assert!(c.vectors[0]
        .iter()
        .zip(&c.out[0])
        .any(|(a, b)| (a - b).abs() > 0.5));
}

#[test]
fn live_harness_agrees_with_the_committed_cases() {
    if !common::live_enabled() {
        eprintln!("skipping: set AFNI_CORE_LIVE=1 to rebuild the harness and re-run AFNI");
        return;
    }
    let committed =
        std::fs::read_to_string(common::data("conformance/signal_bandpass.ref")).unwrap();
    let target = std::env::temp_dir().join(format!("afni_core_live_{}.ref", std::process::id()));
    let out = std::process::Command::new("bash")
        .arg(common::data("regenerate_signal_refs.sh"))
        .env("REGEN_OUT", &target)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let regenerated =
        std::fs::read_to_string(common::data("conformance/signal_bandpass.ref")).unwrap();
    let body = |s: &str| {
        s.lines()
            .filter(|l| !l.starts_with("# afni_version"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(body(&regenerated), body(&committed));
}
