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
// Conformance of `afni_core::suma_colormaps` against AFNI's own
// `MakeColorMap -std <name>`. tests/data/conformance/suma_colormaps.ref holds the
// printed maps (regenerate with tests/data/regenerate_suma_colormaps.sh). The
// program prints two decimals, so each channel must agree within half a unit of
// the last digit. With AFNI_CORE_LIVE=1 the maps are also re-run live.
// ---------------------------------------------------------------------------

mod common;

use afni_core::suma_colormaps::SumaColorMap;

/// Two decimals are printed, so the true value is within 0.005 of the printed one
/// (plus a hair for f32 rounding).
const PRINTED_TOLERANCE: f32 = 0.0051;

/// `name -> rows`, parsed from the fixture.
fn parse_fixture() -> Vec<(String, Vec<[f32; 3]>)> {
    let text = std::fs::read_to_string(common::data("conformance/suma_colormaps.ref")).unwrap();
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|line| {
            let (name, rows) = line.split_once(" | ").unwrap();
            let rows = rows
                .split_whitespace()
                .map(|row| {
                    let v: Vec<f32> = row.split(',').map(|x| x.parse().unwrap()).collect();
                    [v[0], v[1], v[2]]
                })
                .collect();
            (name.to_owned(), rows)
        })
        .collect()
}

fn assert_close(name: &str, ours: &[[f32; 3]], theirs: &[[f32; 3]]) {
    assert_eq!(ours.len(), theirs.len(), "{name}: number of colors");
    for (i, (a, b)) in ours.iter().zip(theirs).enumerate() {
        for c in 0..3 {
            assert!(
                (a[c] - b[c]).abs() <= PRINTED_TOLERANCE,
                "{name} entry {i} channel {c}: afni-core {}, SUMA {}",
                a[c],
                b[c]
            );
        }
    }
}

#[test]
fn every_suma_map_matches_make_color_map() {
    let fixture = parse_fixture();
    // The fixture covers every map we define, and nothing else.
    assert_eq!(fixture.len(), SumaColorMap::ALL.len());
    for (name, rows) in &fixture {
        let map = SumaColorMap::from_name(name).unwrap_or_else(|| panic!("no map named {name}"));
        assert_close(name, &map.colors(), rows);
    }
}

#[test]
fn live_make_color_map_agrees_with_the_committed_maps() {
    if !common::live_enabled() {
        eprintln!("skipping: set AFNI_CORE_LIVE=1 to re-run MakeColorMap");
        return;
    }
    for (name, rows) in parse_fixture() {
        let out = common::run_afni("MakeColorMap", &["-std".to_string(), name.clone()]).unwrap();
        let live: Vec<[f32; 3]> = out
            .lines()
            .filter_map(|l| {
                let v: Vec<f32> = l
                    .split_whitespace()
                    .filter_map(|x| x.parse().ok())
                    .collect();
                (v.len() == 3).then(|| [v[0], v[1], v[2]])
            })
            .collect();
        assert_eq!(live, rows, "{name}: live output differs from the fixture");
    }
}
