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
// Conformance of `afni_core::roi_ops` against SUMA: shortest-path distances
// against `SurfDist` (120 node pairs on three meshes) and ROI growth by distance
// against `ROIgrow -lim` (9 cases). The reference lines are in
// tests/data/conformance/roi_ops.ref (regenerate with
// tests/data/regenerate_roi_refs.sh).
// ---------------------------------------------------------------------------

mod common;

use afni_core::mesh::SurfaceMesh;
use afni_core::roi::NodeSet;
use afni_core::roi_ops::{dilate_by_distance, path_length, shortest_path};

fn mesh(name: &str) -> SurfaceMesh {
    // `common` has no surface reader (reading files is afni-io's job), so parse
    // SUMA's tiny ASCII format here: "n_nodes n_faces", then n_nodes lines of
    // "x y z flag", then n_faces lines of "a b c flag".
    let text = std::fs::read_to_string(common::data(&format!("surfaces/{name}.asc"))).unwrap();
    let mut lines = text.lines().filter(|l| !l.starts_with('#'));
    let counts: Vec<usize> = lines
        .next()
        .unwrap()
        .split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect();
    let vertices: Vec<[f32; 3]> = (0..counts[0])
        .map(|_| {
            let v: Vec<f32> = lines
                .next()
                .unwrap()
                .split_whitespace()
                .map(|x| x.parse().unwrap())
                .collect();
            [v[0], v[1], v[2]]
        })
        .collect();
    let faces: Vec<[u32; 3]> = (0..counts[1])
        .map(|_| {
            let v: Vec<u32> = lines
                .next()
                .unwrap()
                .split_whitespace()
                .map(|x| x.parse().unwrap())
                .collect();
            [v[0], v[1], v[2]]
        })
        .collect();
    SurfaceMesh::from_triangles(vertices, faces).unwrap()
}

fn reference_lines(kind: &str) -> Vec<Vec<String>> {
    std::fs::read_to_string(common::data("conformance/roi_ops.ref"))
        .unwrap()
        .lines()
        .filter(|l| l.starts_with(kind))
        .map(|l| l.split_whitespace().map(str::to_owned).collect())
        .collect()
}

#[test]
fn shortest_path_lengths_match_surfdist() {
    let lines = reference_lines("distance");
    assert_eq!(lines.len(), 120);
    let meshes: Vec<(String, SurfaceMesh)> = ["ico3", "irr3", "irr2"]
        .iter()
        .map(|n| (n.to_string(), mesh(n)))
        .collect();
    for l in &lines {
        let m = &meshes.iter().find(|(n, _)| *n == l[1]).unwrap().1;
        let (from, to): (u32, u32) = (l[2].parse().unwrap(), l[3].parse().unwrap());
        let theirs: f64 = l[4].parse().unwrap();
        let path = shortest_path(m, from, to).unwrap();
        assert_eq!((path[0], *path.last().unwrap()), (from, to));
        let ours = path_length(m, &path).unwrap();
        // SurfDist prints two decimals.
        assert!(
            (ours - theirs).abs() <= 0.0051,
            "{} {from}->{to}: afni-core {ours}, SurfDist {theirs}",
            l[1]
        );
    }
}

/// How many nodes the true shortest-path reach has beyond `ROIgrow`'s, per case
/// (mesh, limit, seeds). On the regular icosahedron the two agree exactly; on the
/// irregular meshes `ROIgrow` stops short of some nodes that a shortest path along
/// the edges reaches within the limit. It is the same approximation `SurfClust -rmm`
/// makes (see the roadmap log, Phase 6), so the difference is pinned, not fixed.
const EXTRA_NODES: [(&str, f64, &str, usize); 9] = [
    ("ico3", 8.0, "10", 0),
    ("ico3", 15.0, "10,300", 0),
    ("ico3", 25.0, "5,6,200", 0),
    ("irr3", 8.0, "10", 0),
    ("irr3", 15.0, "10,300", 0),
    ("irr3", 25.0, "5,6,200", 6),
    ("irr3", 40.0, "77", 11),
    ("irr2", 20.0, "3", 0),
    ("irr2", 30.0, "3,100", 1),
];

#[test]
fn growth_by_distance_matches_roigrow_up_to_its_known_shortfall() {
    let lines = reference_lines("grow");
    assert_eq!(lines.len(), EXTRA_NODES.len());
    for (l, &(mesh_name, lim_pinned, seeds_pinned, extra)) in lines.iter().zip(&EXTRA_NODES) {
        let m = mesh(&l[1]);
        let lim: f64 = l[2].parse().unwrap();
        assert_eq!(
            (l[1].as_str(), lim, l[3].as_str()),
            (mesh_name, lim_pinned, seeds_pinned)
        );
        let seeds = NodeSet::new(l[3].split(',').map(|s| s.parse::<u32>().unwrap()));
        let theirs = NodeSet::new(l[5..].iter().map(|s| s.parse::<u32>().unwrap()));
        let ours = dilate_by_distance(&m, &seeds, lim).unwrap();
        // Everything ROIgrow reaches, a true shortest path reaches too.
        assert!(
            theirs.difference(&ours).is_empty(),
            "{} lim {lim}: ROIgrow reached {:?} but afni-core did not",
            l[1],
            theirs.difference(&ours).as_slice()
        );
        // And each extra node really is within the limit (by the verified path length).
        let extras = ours.difference(&theirs);
        assert_eq!(extras.len(), extra, "{} lim {lim}: extra nodes", l[1]);
        for node in extras.iter() {
            let nearest = seeds
                .iter()
                .map(|s| path_length(&m, &shortest_path(&m, s, node).unwrap()).unwrap())
                .fold(f64::INFINITY, f64::min);
            assert!(
                nearest <= lim + 1e-6,
                "{} node {node}: {nearest} > {lim}",
                l[1]
            );
        }
    }
}

#[test]
fn live_surfdist_and_roigrow_agree_with_the_committed_references() {
    if !common::live_enabled() {
        eprintln!("skipping: set AFNI_CORE_LIVE=1 to re-run SurfDist and ROIgrow");
        return;
    }
    let committed = std::fs::read_to_string(common::data("conformance/roi_ops.ref")).unwrap();
    let target = std::env::temp_dir().join(format!("afni_core_live_{}.ref", std::process::id()));
    let out = std::process::Command::new("bash")
        .arg(common::data("regenerate_roi_refs.sh"))
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
    // The AFNI version line may differ between builds.
    let body = |s: &str| {
        s.lines()
            .filter(|l| !l.starts_with("# afni_version"))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(body(&regenerated), body(&committed));
}
