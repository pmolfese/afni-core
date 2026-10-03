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
// Conformance of `afni_core::volume_cluster` against AFNI's `3dClusterize`.
// tests/data/conformance/volume_clusters.ref holds 18 synthetic cases: the volume,
// the 3dClusterize options, the cluster map the program wrote (`-pref_map`) and the
// rows of its report (regenerate with tests/data/regenerate_volume_clusters.sh).
//
// Each case is re-run through `cluster_volume` and compared:
// * the cluster MAP must be identical, voxel for voxel (this pins connectivity,
//   thresholds, the zero-data rule, masks, size limits and the order of equal-size
//   clusters);
// * every report column must agree to the precision the program prints (one
//   decimal for coordinates, about six significant digits for mean/SEM/peak).
//
// One case is oblique. 3dClusterize still reports cardinal grid coordinates, and so
// does this crate when given the cardinal matrix; the case proves connectivity and
// the report ignore the oblique header.
//
// With AFNI_CORE_LIVE=1 every case is also re-run through the live program.
// ---------------------------------------------------------------------------

mod common;

use afni_core::domain::VolumeDomain;
use afni_core::volume_cluster::{
    cluster_volume, Tails, VolumeClusterInput, VolumeClusterParams, VolumeClusters,
    VoxelConnectivity, VoxelThreshold,
};

const FIXTURE: &str = "conformance/volume_clusters.ref";

/// One parsed case.
struct Case {
    name: String,
    dims: [usize; 3],
    affine: [[f64; 4]; 4],
    options: Vec<String>,
    thr: Vec<f64>,
    dat: Option<Vec<f64>>,
    mask: Option<Vec<bool>>,
    map: Vec<u32>,
    /// Each row: the 16 printed numbers of the report table.
    rows: Vec<Vec<f64>>,
}

fn numbers(text: &str) -> Vec<f64> {
    text.split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect()
}

fn parse_cases() -> Vec<Case> {
    let text = std::fs::read_to_string(common::data(FIXTURE)).unwrap();
    let mut cases = Vec::new();
    let mut current: Option<Case> = None;
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let (key, rest) = line.split_once(' ').unwrap_or((line, ""));
        match key {
            "case" => {
                current = Some(Case {
                    name: rest.to_owned(),
                    dims: [0; 3],
                    affine: [[0.0; 4]; 4],
                    options: Vec::new(),
                    thr: Vec::new(),
                    dat: None,
                    mask: None,
                    map: Vec::new(),
                    rows: Vec::new(),
                });
            }
            "end" => cases.push(current.take().unwrap()),
            _ => {
                let c = current.as_mut().expect("line outside a case");
                match key {
                    "dims" => {
                        let v = numbers(rest);
                        c.dims = [v[0] as usize, v[1] as usize, v[2] as usize];
                    }
                    "affine" => {
                        let v = numbers(rest);
                        for r in 0..3 {
                            for k in 0..4 {
                                c.affine[r][k] = v[r * 4 + k];
                            }
                        }
                        c.affine[3] = [0.0, 0.0, 0.0, 1.0];
                    }
                    "options" => c.options = rest.split_whitespace().map(str::to_owned).collect(),
                    "thr" => c.thr = numbers(rest),
                    "dat" => c.dat = (rest != "-").then(|| numbers(rest)),
                    "mask" => {
                        c.mask =
                            (rest != "-").then(|| numbers(rest).iter().map(|&m| m != 0.0).collect())
                    }
                    "map" => c.map = numbers(rest).iter().map(|&v| v as u32).collect(),
                    "row" => c.rows.push(numbers(rest)),
                    other => panic!("unknown fixture key {other}"),
                }
            }
        }
    }
    cases
}

/// AFNI reads thresholds as 32-bit floats and compares them with 32-bit data; do
/// the same so a threshold such as 1.9 cannot differ in the last bits.
fn f32_exact(text: &str) -> f64 {
    f64::from(text.parse::<f32>().unwrap())
}

/// The cluster parameters equivalent to a set of `3dClusterize` options.
fn params_for(case: &Case) -> VolumeClusterParams {
    let mut nn = 1;
    let mut threshold = None;
    let mut tails = Tails::Separate;
    let mut min_voxels = None;
    let mut it = case.options.iter();
    while let Some(opt) = it.next() {
        let mut next = || it.next().expect("option value").as_str();
        match opt.as_str() {
            "-NN" => nn = next().parse().unwrap(),
            "-1sided" => {
                let side = next();
                let t = f32_exact(next());
                threshold = Some(match side {
                    "RIGHT_TAIL" => VoxelThreshold::RightTail(t),
                    "LEFT_TAIL" => VoxelThreshold::LeftTail(t),
                    other => panic!("side {other}"),
                });
            }
            "-2sided" | "-bisided" => {
                tails = if opt == "-2sided" {
                    Tails::Merged
                } else {
                    Tails::Separate
                };
                let (l, r) = (f32_exact(next()), f32_exact(next()));
                threshold = Some(VoxelThreshold::TwoSided {
                    left_upper: l,
                    right_lower: r,
                });
            }
            "-within_range" => {
                let (lo, hi) = (f32_exact(next()), f32_exact(next()));
                threshold = Some(VoxelThreshold::WithinRange { lo, hi });
            }
            "-clust_nvox" => min_voxels = Some(next().parse().unwrap()),
            other => panic!("{}: unhandled option {other}", case.name),
        }
    }
    let mut params = VolumeClusterParams::new(
        VoxelConnectivity::from_nn(nn).unwrap(),
        threshold.expect("a threshold option"),
    );
    params.tails = tails;
    params.min_voxels = min_voxels;
    params
}

fn run(case: &Case) -> VolumeClusters {
    let domain = VolumeDomain::new(None, case.dims, Some(case.affine)).unwrap();
    cluster_volume(
        &VolumeClusterInput {
            domain: &domain,
            threshold_values: &case.thr,
            data_values: case.dat.as_deref(),
            mask: case.mask.as_deref(),
        },
        &params_for(case),
    )
    .unwrap_or_else(|e| panic!("{}: {e}", case.name))
}

/// Report numbers print with one decimal (coordinates) or about six significant
/// digits (mean, SEM, peak value).
const COORD_TOLERANCE: f64 = 0.0501;
const VALUE_TOLERANCE: f64 = 5.0e-4;

fn close(what: &str, case: &str, row: usize, ours: f64, theirs: f64, tolerance: f64) {
    assert!(
        (ours - theirs).abs() <= tolerance,
        "{case} cluster {}: {what}: afni-core {ours}, 3dClusterize {theirs}",
        row + 1
    );
}

fn compare(case: &Case, got: &VolumeClusters) {
    // The cluster map: exact.
    assert_eq!(got.labels, case.map, "{}: cluster map", case.name);
    // The report: one row per cluster, in rank order.
    assert_eq!(
        got.clusters.len(),
        case.rows.len(),
        "{}: clusters",
        case.name
    );
    for (r, (c, row)) in got.clusters.iter().zip(&case.rows).enumerate() {
        let name = case.name.as_str();
        // Column 0 is labelled "Volume" but 3dClusterize prints the VOXEL COUNT.
        assert_eq!(
            c.voxel_count as f64,
            row[0],
            "{name} cluster {}: voxels",
            r + 1
        );
        let com = c.center_of_mass_world.unwrap();
        let (lo, hi) = c.bounds_world.unwrap();
        let peak = c.peak_world.unwrap();
        for axis in 0..3 {
            close(
                "center of mass",
                name,
                r,
                com[axis],
                row[1 + axis],
                COORD_TOLERANCE,
            );
            close(
                "minimum",
                name,
                r,
                lo[axis],
                row[4 + 2 * axis],
                COORD_TOLERANCE,
            );
            close(
                "maximum",
                name,
                r,
                hi[axis],
                row[5 + 2 * axis],
                COORD_TOLERANCE,
            );
            close(
                "peak position",
                name,
                r,
                peak[axis],
                row[13 + axis],
                COORD_TOLERANCE,
            );
        }
        close("mean", name, r, c.mean, row[10], VALUE_TOLERANCE);
        close("SEM", name, r, c.std_error, row[11], VALUE_TOLERANCE);
        close("peak value", name, r, c.peak.1, row[12], VALUE_TOLERANCE);
    }
}

#[test]
fn matches_3dclusterize_map_and_report() {
    let cases = parse_cases();
    assert_eq!(cases.len(), 18);
    let mut clusters = 0;
    for case in &cases {
        let got = run(case);
        compare(case, &got);
        clusters += got.clusters.len();
    }
    assert!(clusters > 500, "only {clusters} clusters were compared");
}

#[test]
fn the_cases_exercise_every_option() {
    let cases = parse_cases();
    for needle in [
        "right_nn1",
        "right_nn3",
        "left_nn1",
        "twosided_nn2",
        "bisided_nn2",
        "bisided_nvox",
        "within_range",
        "idat",
        "mask",
        "oblique",
        "ties_bisided",
        "ties_twosided",
        "no_clusters",
    ] {
        assert!(
            cases.iter().any(|c| c.name.starts_with(needle)),
            "no case named {needle}*"
        );
    }
    // Larger neighborhoods really do change the result in the reference.
    let by = |name: &str| cases.iter().find(|c| c.name == name).unwrap();
    assert_ne!(by("right_nn1").map, by("right_nn3").map);
    // Separate tails really differ from merged tails.
    assert_ne!(by("twosided_nn2").map, by("bisided_nn2").map);
    // The no-cluster case really has none.
    assert!(by("no_clusters").rows.is_empty());
    // The data and mask cases really use their extra inputs.
    assert!(by("idat").dat.is_some() && by("mask").mask.is_some());
}

#[test]
fn live_3dclusterize_agrees_with_the_committed_cases() {
    if !common::live_enabled() {
        eprintln!("skipping: set AFNI_CORE_LIVE=1 to re-run 3dClusterize");
        return;
    }
    // Re-running needs the volumes rebuilt with 3dUndump; the regenerate script does
    // exactly that, so run it into a temporary file (REGEN_OUT) and compare with the committed one.
    let script = common::data("regenerate_volume_clusters.sh");
    let committed = std::fs::read_to_string(common::data(FIXTURE)).unwrap();
    let target = std::env::temp_dir().join(format!("afni_core_live_{}.ref", std::process::id()));
    let out = std::process::Command::new("bash")
        .arg(&script)
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
    // Ignore the version line, which differs between AFNI builds.
    let body = |s: &str| {
        s.lines()
            .filter(|l| !l.starts_with("# afni_version"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(body(&regenerated), body(&committed));
}
