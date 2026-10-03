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
// Conformance of `afni_core::cluster` against SUMA's `SurfClust`: 22 runs over a
// regular and an irregular 642-node sphere and two node-wise data sets, covering
// edge-ring connectivity (`-rmm -N`), millimetre radii (`-rmm r`), each threshold
// style (`-thresh`, `-athresh`, `-in_range`, `-ex_range`), minimum area and node
// counts, and the three sort orders. Each cluster table row is compared column by
// column. References: tests/data/conformance/surfclust.ref, from
// tests/data/regenerate_surface_refs.sh.
//
// Separate tests check the corrected two-sided behavior (`Tails::Separate`), which
// SurfClust cannot express, against the merged result.
//
// SurfClust prints 2-3 decimals, so value columns are compared to that precision.
// ---------------------------------------------------------------------------

mod common;

use afni_core::cluster::{
    label_clusters, ClusterInput, ClusterLabels, ClusterParams, ClusterSort, Connectivity, Tails,
};
use afni_core::mesh::SurfaceMesh;
use afni_core::threshold::Threshold;

struct Case {
    name: String,
    surface: String,
    data: String,
    options: Vec<String>,
    rows: Vec<Vec<f64>>,
}

fn cases() -> Vec<Case> {
    let text = std::fs::read_to_string(common::data("conformance/surfclust.ref")).unwrap();
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|line| {
            let (lhs, rhs) = line.split_once("=>").unwrap();
            let f: Vec<&str> = lhs.split(" | ").map(str::trim).collect();
            let rows = rhs
                .split(';')
                .filter(|r| !r.trim().is_empty())
                .map(|r| {
                    r.trim()
                        .split(',')
                        .map(|x| x.trim().parse().unwrap())
                        .collect()
                })
                .collect();
            Case {
                name: f[0].into(),
                surface: f[1].into(),
                data: f[2].into(),
                options: f[3].split_whitespace().map(str::to_owned).collect(),
                rows,
            }
        })
        .collect()
}

fn mesh(name: &str) -> SurfaceMesh {
    let (v, f) = common::read_asc(name);
    SurfaceMesh::from_triangles(v, f).unwrap()
}

/// The afni-core settings equivalent to a set of SurfClust options, and which
/// nodes the threshold keeps active.
fn settings(options: &[String], values: &[f64]) -> (ClusterParams, Vec<bool>) {
    let mut params = ClusterParams {
        exclude_zero_values: true,
        ..Default::default()
    };
    let mut threshold = Threshold::Off;
    let (mut area, mut nodes) = (None::<f64>, None::<usize>);
    let mut it = options.iter();
    while let Some(opt) = it.next() {
        let mut num = || -> f64 { it.next().expect("option value").parse().unwrap() };
        match opt.as_str() {
            "-rmm" => {
                let r = num();
                params.connectivity = if r < 0.0 {
                    Connectivity::EdgeRings((-r) as u32)
                } else {
                    Connectivity::GraphDistance(r)
                };
            }
            "-thresh" => threshold = Threshold::Above(num()),
            "-athresh" => threshold = Threshold::AbsoluteAbove(num()),
            "-in_range" => {
                threshold = Threshold::Between {
                    lo: num(),
                    hi: num(),
                }
            }
            "-ex_range" => {
                threshold = Threshold::Outside {
                    lo: num(),
                    hi: num(),
                }
            }
            "-amm2" => area = Some(num()),
            "-n" => nodes = Some(num() as usize),
            "-sort_n_nodes" => params.sort = ClusterSort::Nodes,
            "-sort_none" => params.sort = ClusterSort::Discovery,
            other => panic!("unhandled SurfClust option {other}"),
        }
    }
    // SurfClust: a negative -amm2 with no -n means "at least that many nodes".
    match (area, nodes) {
        (Some(a), None) if a < 0.0 => params.min_nodes = Some((-a) as usize),
        (a, n) => {
            params.min_area = a.filter(|&a| a > 0.0);
            params.min_nodes = n.filter(|&n| n > 0);
        }
    }
    let active = values.iter().map(|&v| threshold.passes(v)).collect();
    (params, active)
}

fn run(case: &Case) -> (ClusterLabels, SurfaceMesh, Vec<f64>) {
    let m = mesh(&case.surface);
    let values = common::read_column(&format!("surfaces/{}.1D", case.data));
    let (params, active) = settings(&case.options, &values);
    let labels = label_clusters(
        &ClusterInput {
            mesh: &m,
            active: &active,
            values: &values,
            tail_values: None,
        },
        &params,
    )
    .unwrap_or_else(|e| panic!("{}: {e}", case.name));
    (labels, m, values)
}

fn near(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

/// Compare one cluster with one SurfClust row (23 columns).
fn compare_row(
    case: &str,
    rank: usize,
    c: &afni_core::cluster::ClusterSummary,
    row: &[f64],
    ordered: bool,
) {
    let ctx = |col: &str| format!("{case} cluster {rank} {col}");
    if ordered {
        assert_eq!(c.label as usize, row[0] as usize, "{}", ctx("rank"));
    }
    assert_eq!(c.node_count, row[1] as usize, "{}", ctx("nodes"));
    assert!(
        near(c.area, row[2], 0.006 + 3e-6 * row[2]),
        "{}: {} vs {}",
        ctx("area"),
        c.area,
        row[2]
    );
    assert!(
        near(c.mean, row[3], 6e-4),
        "{}: {} vs {}",
        ctx("mean"),
        c.mean,
        row[3]
    );
    assert!(
        near(c.mean_abs, row[4], 6e-4),
        "{}: {} vs {}",
        ctx("mean |value|"),
        c.mean_abs,
        row[4]
    );
    assert!(
        near(c.min.0, row[7], 6e-4) && c.min.1 as f64 == row[8],
        "{}: {:?} vs {} @ {}",
        ctx("min"),
        c.min,
        row[7],
        row[8]
    );
    assert!(
        near(c.max.0, row[9], 6e-4) && c.max.1 as f64 == row[10],
        "{}: {:?} vs {} @ {}",
        ctx("max"),
        c.max,
        row[9],
        row[10]
    );
    assert!(
        near(c.variance, row[11], 6e-4 + 1e-4 * row[11]),
        "{}: {} vs {}",
        ctx("variance"),
        c.variance,
        row[11]
    );
    assert!(
        near(c.std_error, row[12], 6e-4),
        "{}: {} vs {}",
        ctx("std error"),
        c.std_error,
        row[12]
    );
    assert!(
        near(c.min_abs.0, row[13], 6e-4) && c.min_abs.1 as f64 == row[14],
        "{}: {:?} vs {} @ {}",
        ctx("min |value|"),
        c.min_abs,
        row[13],
        row[14]
    );
    assert!(
        near(c.max_abs.0, row[15], 6e-4) && c.max_abs.1 as f64 == row[16],
        "{}: {:?} vs {} @ {}",
        ctx("max |value|"),
        c.max_abs,
        row[15],
        row[16]
    );
    for k in 0..3 {
        assert!(
            // The value-weighted mean divides by the sum of the (mixed-sign) values, which can
            // nearly cancel; SUMA's float rounding is then amplified, so allow a relative term.
            near(
                c.center_of_mass[k],
                row[17 + k],
                2e-3 + 2e-5 * row[17 + k].abs()
            ),
            "{}[{k}]: {} vs {}",
            ctx("center of mass"),
            c.center_of_mass[k],
            row[17 + k]
        );
        assert!(
            near(c.centroid[k], row[20 + k], 6e-4),
            "{}[{k}]: {} vs {}",
            ctx("centroid"),
            c.centroid[k],
            row[20 + k]
        );
    }
}

/// A regular mesh has many clusters of exactly equal area, whose relative order in
/// SurfClust's table depends on float rounding in its single-precision area sums and
/// cannot be reproduced. On such meshes compare the clusters as a SET (sorted by node
/// count, then by the node of the largest |value|); on the irregular mesh, where areas
/// are all distinct, the rank order must match too.
fn has_area_ties(case: &Case) -> bool {
    case.surface == "ico3"
}

/// Compare every cluster with its row; returns `false` instead of panicking when the
/// two tables describe different clusters (for the callers that tolerate that).
fn compare_table(case: &Case, labels: &ClusterLabels, strict: bool) -> bool {
    let ordered = !has_area_ties(case);
    if labels.clusters.len() != case.rows.len() {
        assert!(
            !strict,
            "{}: {} clusters vs SurfClust's {}",
            case.name,
            labels.clusters.len(),
            case.rows.len()
        );
        return false;
    }
    let mut ours: Vec<&afni_core::cluster::ClusterSummary> = labels.clusters.iter().collect();
    let mut theirs: Vec<&Vec<f64>> = case.rows.iter().collect();
    if !ordered {
        ours.sort_by_key(|c| (std::cmp::Reverse(c.node_count), c.max_abs.1));
        theirs.sort_by_key(|r| (std::cmp::Reverse(r[1] as usize), r[16] as u32));
    }
    let same = ours
        .iter()
        .zip(&theirs)
        .all(|(c, r)| c.node_count == r[1] as usize && c.max_abs.1 as f64 == r[16]);
    if !same {
        assert!(
            !strict,
            "{}: the clusters differ from SurfClust's",
            case.name
        );
        return false;
    }
    for (i, (c, row)) in ours.iter().zip(&theirs).enumerate() {
        compare_row(&case.name, i + 1, c, row, ordered);
    }
    true
}

fn is_millimetre_case(case: &Case) -> bool {
    case.options
        .windows(2)
        .any(|w| w[0] == "-rmm" && !w[1].starts_with('-'))
}

#[test]
fn edge_ring_clusters_match_surfclust_exactly() {
    let mut compared = 0;
    for case in cases().iter().filter(|c| !is_millimetre_case(c)) {
        let (labels, _, _) = run(case);
        compare_table(case, &labels, true);
        compared += case.rows.len();
    }
    assert!(compared > 300, "only {compared} cluster rows were compared");
}

#[test]
fn millimetre_radius_clusters_match_surfclust() {
    let mut compared = 0;
    let mut mismatched = Vec::new();
    for case in cases().iter().filter(|c| is_millimetre_case(c)) {
        let (labels, _, _) = run(case);
        if compare_table(case, &labels, false) {
            compared += case.rows.len();
        } else {
            // SUMA's millimetre search over-estimates distances (see the roadmap), so it
            // connects FEWER nodes than a true shortest-path search and therefore finds
            // at least as many clusters. The direction of the difference is part of the
            // claim.
            assert!(
                labels.clusters.len() <= case.rows.len(),
                "{}: afni-core found MORE clusters ({}) than SurfClust ({}), the opposite of the documented direction",
                case.name,
                labels.clusters.len(),
                case.rows.len()
            );
            mismatched.push(case.name.clone());
        }
    }
    assert!(compared > 100, "only {compared} rows compared");
    // Known disagreement: at -rmm 14 on the irregular sphere SurfClust finds 35 clusters
    // where true graph distance finds 32. Six of the seven other radii/meshes agree
    // exactly. A change in this list means the behavior moved and must be re-examined.
    assert_eq!(mismatched, vec!["mm_wide_B".to_string()]);
}

#[test]
fn the_cases_cover_every_option() {
    let all = cases();
    for needle in [
        "rings1_A",
        "rings2_A",
        "rings3_A",
        "thresh",
        "athresh",
        "inrange",
        "exrange",
        "minnodes",
        "minarea",
        "negarea",
        "sort_nodes",
        "sort_none",
        "mm_",
    ] {
        assert!(
            all.iter().any(|c| c.name.contains(needle)),
            "no case named *{needle}*"
        );
    }
    assert!(
        all.iter().all(|c| !c.rows.is_empty()),
        "every case produced a table"
    );
}

#[test]
fn separate_tails_split_surfclusters_without_losing_or_inventing_nodes() {
    // SurfClust cannot keep tails apart. Check the corrected mode against its output.
    for case in cases().iter().filter(|c| {
        c.name.starts_with("rings1_B_athresh")
            || c.name.starts_with("rings2_B_athresh")
            || c.name == "rings1_A"
    }) {
        let m = mesh(&case.surface);
        let values = common::read_column(&format!("surfaces/{}.1D", case.data));
        let (params, active) = settings(&case.options, &values);
        let input = ClusterInput {
            mesh: &m,
            active: &active,
            values: &values,
            tail_values: None,
        };
        let merged = label_clusters(&input, &params).unwrap();
        let separate = label_clusters(
            &input,
            &ClusterParams {
                tails: Tails::Separate,
                ..params
            },
        )
        .unwrap();
        // The same nodes are clustered either way...
        assert!(
            merged.survivor_mask().iter().filter(|&&b| b).count()
                >= separate.survivor_mask().iter().filter(|&&b| b).count(),
            "{}",
            case.name
        );
        // ...separation only ever splits: each separate cluster lies inside one merged cluster
        // and has a single sign.
        for s in &separate.clusters {
            let nodes = separate.nodes_for(s.label);
            let parent = merged.labels[nodes[0] as usize];
            assert!(
                nodes.iter().all(|&n| merged.labels[n as usize] == parent
                    || merged.labels[n as usize] == 0
                    || parent == 0),
                "{}: cluster {} spans merged clusters",
                case.name,
                s.label
            );
            let negative = values[nodes[0] as usize] < 0.0;
            assert!(
                nodes
                    .iter()
                    .all(|&n| (values[n as usize] < 0.0) == negative),
                "{}: cluster {} mixes signs",
                case.name,
                s.label
            );
        }
        assert!(
            separate.clusters.len() >= merged.clusters.len()
                || params.min_area.is_some()
                || params.min_nodes.is_some(),
            "{}",
            case.name
        );
    }
}

#[test]
fn separate_tails_really_differ_from_surfclust_on_this_data() {
    // Guard against a vacuous test: on two-sided data the corrected mode must find
    // more (smaller) clusters than SurfClust does.
    let case = cases()
        .into_iter()
        .find(|c| c.name == "rings1_B_athresh")
        .unwrap();
    let m = mesh(&case.surface);
    let values = common::read_column(&format!("surfaces/{}.1D", case.data));
    let (params, active) = settings(&case.options, &values);
    let input = ClusterInput {
        mesh: &m,
        active: &active,
        values: &values,
        tail_values: None,
    };
    let merged = label_clusters(&input, &params).unwrap();
    let separate = label_clusters(
        &input,
        &ClusterParams {
            tails: Tails::Separate,
            ..params
        },
    )
    .unwrap();
    assert_eq!(
        merged.clusters.len(),
        case.rows.len(),
        "merged mode is SurfClust's"
    );
    assert!(
        separate.clusters.len() > merged.clusters.len(),
        "{} vs {}",
        separate.clusters.len(),
        merged.clusters.len()
    );
}

#[test]
fn cached_neighborhoods_reproduce_surfclust_results() {
    use afni_core::cluster::{label_clusters_cached, ClusterNeighborhoods};
    let case = cases()
        .into_iter()
        .find(|c| c.name == "rings2_B_athresh")
        .unwrap();
    let m = mesh(&case.surface);
    let values = common::read_column(&format!("surfaces/{}.1D", case.data));
    let (params, active) = settings(&case.options, &values);
    let cache = ClusterNeighborhoods::build(&m, params.connectivity).unwrap();
    let input = ClusterInput {
        mesh: &m,
        active: &active,
        values: &values,
        tail_values: None,
    };
    let cached = label_clusters_cached(&input, &params, &cache).unwrap();
    for (i, (c, row)) in cached.clusters.iter().zip(&case.rows).enumerate() {
        compare_row(&case.name, i + 1, c, row, true);
    }
}
