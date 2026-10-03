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
// A timing harness for surface clustering on realistic mesh sizes (the roadmap's
// "benchmark wide ring/radius searches" item). It builds icosahedron spheres of
// increasing resolution, marks about a third of the nodes active with
// deterministic pseudo-random values, and times:
//
//   * one clustering with neighborhoods computed on the fly, and
//   * building a `ClusterNeighborhoods` cache once, then re-clustering twenty
//     times (the situation while a threshold slider moves).
//
//     cargo run --release --example cluster_bench
//
// Always run it with --release; debug builds are an order of magnitude slower.
// ---------------------------------------------------------------------------

use std::collections::HashMap;
use std::time::Instant;

use afni_core::cluster::{
    label_clusters, label_clusters_cached, ClusterInput, ClusterNeighborhoods, ClusterParams,
    Connectivity,
};
use afni_core::mesh::SurfaceMesh;

/// A unit-ish icosahedron subdivided `levels` times and scaled to `radius`.
fn icosphere(levels: u32, radius: f32) -> SurfaceMesh {
    let t = (1.0 + 5.0_f32.sqrt()) / 2.0;
    let mut v: Vec<[f32; 3]> = vec![
        [-1.0, t, 0.0],
        [1.0, t, 0.0],
        [-1.0, -t, 0.0],
        [1.0, -t, 0.0],
        [0.0, -1.0, t],
        [0.0, 1.0, t],
        [0.0, -1.0, -t],
        [0.0, 1.0, -t],
        [t, 0.0, -1.0],
        [t, 0.0, 1.0],
        [-t, 0.0, -1.0],
        [-t, 0.0, 1.0],
    ];
    let mut f: Vec<[u32; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for _ in 0..levels {
        let mut midpoint: HashMap<(u32, u32), u32> = HashMap::new();
        let mut next = Vec::with_capacity(f.len() * 4);
        let mut mid = |a: u32, b: u32, v: &mut Vec<[f32; 3]>| -> u32 {
            *midpoint.entry((a.min(b), a.max(b))).or_insert_with(|| {
                let (p, q) = (v[a as usize], v[b as usize]);
                v.push([
                    (p[0] + q[0]) / 2.0,
                    (p[1] + q[1]) / 2.0,
                    (p[2] + q[2]) / 2.0,
                ]);
                v.len() as u32 - 1
            })
        };
        for tri in &f {
            let (a, b, c) = (
                mid(tri[0], tri[1], &mut v),
                mid(tri[1], tri[2], &mut v),
                mid(tri[2], tri[0], &mut v),
            );
            next.extend([[tri[0], a, c], [tri[1], b, a], [tri[2], c, b], [a, b, c]]);
        }
        f = next;
    }
    for p in &mut v {
        let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        *p = [p[0] / n * radius, p[1] / n * radius, p[2] / n * radius];
    }
    SurfaceMesh::from_triangles(v, f).expect("icosphere is a valid mesh")
}

fn main() {
    println!(
        "{:>8} {:>22} {:>12} {:>12} {:>14} {:>9}",
        "nodes", "connectivity", "one run", "cache build", "20 cached runs", "clusters"
    );
    for levels in [4, 5, 6] {
        let mesh = icosphere(levels, 100.0);
        let n = mesh.vertices().len();
        // Deterministic pseudo-random values; a third of the nodes are active.
        let mut s: u64 = 42;
        let mut rng = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (s >> 33) as f64 / (1u64 << 31) as f64
        };
        let values: Vec<f64> = (0..n).map(|_| (rng() - 0.5) * 8.0).collect();
        let active: Vec<bool> = values.iter().map(|v| v.abs() > 2.6).collect();
        let input = ClusterInput {
            mesh: &mesh,
            active: &active,
            values: &values,
            tail_values: None,
        };
        // Mean edge length sets a sensible millimetre radius.
        let mean_edge = mesh.edge_lengths().iter().sum::<f64>() / mesh.edge_lengths().len() as f64;
        for (label, conn) in [
            ("1 ring", Connectivity::EdgeRings(1)),
            ("3 rings", Connectivity::EdgeRings(3)),
            ("6 rings", Connectivity::EdgeRings(6)),
            (
                "radius 2 edges",
                Connectivity::GraphDistance(2.0 * mean_edge),
            ),
            (
                "radius 5 edges",
                Connectivity::GraphDistance(5.0 * mean_edge),
            ),
        ] {
            let params = ClusterParams {
                connectivity: conn,
                ..Default::default()
            };
            let t0 = Instant::now();
            let first = label_clusters(&input, &params).unwrap();
            let one = t0.elapsed();
            let t1 = Instant::now();
            let cache = ClusterNeighborhoods::build(&mesh, conn).unwrap();
            let build = t1.elapsed();
            let t2 = Instant::now();
            for _ in 0..20 {
                let r = label_clusters_cached(&input, &params, &cache).unwrap();
                assert_eq!(r.clusters.len(), first.clusters.len());
            }
            let repeated = t2.elapsed();
            println!(
                "{:>8} {:>22} {:>10.1?} {:>12.1?} {:>14.1?} {:>9}",
                n,
                label,
                one,
                build,
                repeated,
                first.clusters.len()
            );
        }
    }
}
