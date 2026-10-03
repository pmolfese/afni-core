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
// Tractography results: tracts (ordered 3D points along a fiber), grouped into
// bundles (usually all the tracts that join one pair of regions), with the numbers a
// viewer or a script needs from them: polyline length, tangent directions, a
// bounding box, and selection by bundle, id, length, point count or region.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `afni-io` reads and writes FATCAT's `.niml.tract` files (the
//   `TAYLOR_TRACT_DATUM` records of a `<network>`); its `adapt` module converts to
//   and from `TractSet`. Nothing here knows about NIML.
// * Length selection uses `threshold::Threshold`, the same type that selects surface
//   nodes and graph edges.
// * `graph.rs` holds the network these tracts connect. A graph file may link to a
//   tract file; that link is file business and stays in `afni-io`.
//
// COORDINATES
//
// Points are stored as the file has them: AFNI's DICOM frame (+x toward the patient's
// left, +y posterior, which AFNI calls "RAI"), in millimetres. `to_ras` converts to
// the frame surfaces and NIfTI use (see `domain::flip_dicom_ras`). Length and
// direction do not depend on the frame; a bounding box does.
//
// `length` matches AFNI's `Tract_Length` (TrackIO.c): the sum of the straight
// distances between consecutive points.
// ---------------------------------------------------------------------------

//! Tracts, bundles and the selections a viewer needs.

use crate::domain::flip_dicom_ras;
use crate::error::{Error, Result};
use crate::threshold::Threshold;

/// One fiber: an id and its points in order.
#[derive(Debug, Clone, PartialEq)]
pub struct Tract {
    /// The tract's id (FATCAT numbers them within a network).
    pub id: i32,
    /// The points, in order along the fiber. At least one, all finite.
    pub points: Vec<[f32; 3]>,
}

impl Tract {
    /// A tract; errors for no points or a non-finite coordinate.
    pub fn new(id: i32, points: Vec<[f32; 3]>) -> Result<Self> {
        if points.is_empty() {
            return Err(Error::Empty(format!("tract {id} points")));
        }
        for p in &points {
            for &c in p {
                crate::numeric::ensure_finite("tract coordinate", f64::from(c))?;
            }
        }
        Ok(Self { id, points })
    }

    /// Length of the polyline through the points (AFNI's `Tract_Length`); 0 for a
    /// single point. Summed in `f64`.
    pub fn length(&self) -> f64 {
        self.points.windows(2).map(|w| distance(w[0], w[1])).sum()
    }

    /// The distance along the tract to each point (the first is 0).
    pub fn arc_lengths(&self) -> Vec<f64> {
        let mut total = 0.0;
        let mut out = Vec::with_capacity(self.points.len());
        out.push(0.0);
        for w in self.points.windows(2) {
            total += distance(w[0], w[1]);
            out.push(total);
        }
        out
    }

    /// Unit direction of travel at each point: the central difference of its
    /// neighbors (one-sided at the ends). A point where the tract does not move (a
    /// repeated point, or a single-point tract) has direction `[0, 0, 0]`.
    pub fn tangents(&self) -> Vec<[f32; 3]> {
        let n = self.points.len();
        (0..n)
            .map(|i| {
                let (a, b) = (
                    self.points[i.saturating_sub(1)],
                    self.points[(i + 1).min(n - 1)],
                );
                let d = [
                    f64::from(b[0]) - f64::from(a[0]),
                    f64::from(b[1]) - f64::from(a[1]),
                    f64::from(b[2]) - f64::from(a[2]),
                ];
                let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                if len == 0.0 {
                    [0.0; 3]
                } else {
                    [
                        (d[0] / len) as f32,
                        (d[1] / len) as f32,
                        (d[2] / len) as f32,
                    ]
                }
            })
            .collect()
    }
}

fn distance(a: [f32; 3], b: [f32; 3]) -> f64 {
    let d: [f64; 3] = std::array::from_fn(|k| f64::from(a[k]) - f64::from(b[k]));
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

/// The tracts that belong together (usually those joining one pair of regions).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TractBundle {
    /// The bundle's tag (FATCAT's `Bundle_Tag`: the network edge number).
    pub tag: Option<i32>,
    /// A second tag (`Bundle_Alt_Tag`).
    pub alt_tag: Option<i32>,
    /// Labels of the two end regions (`Bundle_Ends`).
    pub ends: Option<String>,
    /// The tracts.
    pub tracts: Vec<Tract>,
}

/// A box around some points, plus the sphere (about its center) that holds them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpatialBounds {
    /// Smallest coordinate on each axis.
    pub min: [f32; 3],
    /// Largest coordinate on each axis.
    pub max: [f32; 3],
    /// The middle of the box.
    pub center: [f32; 3],
    /// The largest half-extent over the three axes (at least `f32::EPSILON`, so it is
    /// safe to divide by).
    pub radius: f32,
}

impl SpatialBounds {
    /// The bounds of some points, or `None` if there are none.
    pub fn from_points<'a>(points: impl IntoIterator<Item = &'a [f32; 3]>) -> Option<Self> {
        let mut it = points.into_iter();
        let first = *it.next()?;
        let (mut min, mut max) = (first, first);
        for p in it {
            for k in 0..3 {
                min[k] = min[k].min(p[k]);
                max[k] = max[k].max(p[k]);
            }
        }
        let center: [f32; 3] = std::array::from_fn(|k| (min[k] + max[k]) * 0.5);
        let radius = (0..3)
            .map(|k| (max[k] - center[k]).abs())
            .fold(0.0, f32::max)
            .max(f32::EPSILON);
        Some(Self {
            min,
            max,
            center,
            radius,
        })
    }
}

/// A set of bundles: what one `.niml.tract` network file holds.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TractSet {
    /// The bundles, in file order.
    pub bundles: Vec<TractBundle>,
}

impl TractSet {
    /// Number of tracts in all bundles.
    pub fn tract_count(&self) -> usize {
        self.bundles.iter().map(|b| b.tracts.len()).sum()
    }

    /// Number of points in all tracts.
    pub fn point_count(&self) -> usize {
        self.iter_tracts().map(|t| t.points.len()).sum()
    }

    /// Every tract, bundle by bundle.
    pub fn iter_tracts(&self) -> impl Iterator<Item = &Tract> {
        self.bundles.iter().flat_map(|b| b.tracts.iter())
    }

    /// The bounds of every point, or `None` if there are none.
    pub fn bounds(&self) -> Option<SpatialBounds> {
        SpatialBounds::from_points(self.iter_tracts().flat_map(|t| t.points.iter()))
    }

    /// The same tracts with every point converted between AFNI's DICOM frame and RAS
    /// (the conversion is its own inverse).
    pub fn flipped(&self) -> TractSet {
        let mut out = self.clone();
        for t in out.bundles.iter_mut().flat_map(|b| b.tracts.iter_mut()) {
            t.points.iter_mut().for_each(|p| *p = flip_dicom_ras(*p));
        }
        out
    }

    /// Keep the tracts for which `keep` is true (bundles keep their order and labels;
    /// a bundle left empty stays, so tags are not lost; call
    /// [`without_empty_bundles`](Self::without_empty_bundles) to drop them).
    pub fn filter(&self, mut keep: impl FnMut(&TractBundle, &Tract) -> bool) -> TractSet {
        TractSet {
            bundles: self
                .bundles
                .iter()
                .map(|b| TractBundle {
                    tracts: b.tracts.iter().filter(|t| keep(b, t)).cloned().collect(),
                    ..b.clone()
                })
                .collect(),
        }
    }

    /// Drop bundles that have no tracts.
    pub fn without_empty_bundles(mut self) -> TractSet {
        self.bundles.retain(|b| !b.tracts.is_empty());
        self
    }

    /// Only the bundles with the given tag.
    pub fn with_bundle_tag(&self, tag: i32) -> TractSet {
        TractSet {
            bundles: self
                .bundles
                .iter()
                .filter(|b| b.tag == Some(tag))
                .cloned()
                .collect(),
        }
    }

    /// Only the tracts with one of the given ids.
    pub fn with_ids(&self, ids: &[i32]) -> TractSet {
        self.filter(|_, t| ids.contains(&t.id))
    }

    /// Only the tracts whose length passes `threshold` (for example
    /// `Threshold::Between { lo: 20.0, hi: 120.0 }`). Lengths are in the units of the
    /// coordinates, millimetres for FATCAT files.
    pub fn with_length(&self, threshold: &Threshold) -> Result<TractSet> {
        threshold.validate()?;
        Ok(self.filter(|_, t| threshold.passes(t.length())))
    }

    /// Only the tracts with at least `min_points` points.
    pub fn with_min_points(&self, min_points: usize) -> TractSet {
        self.filter(|_, t| t.points.len() >= min_points)
    }

    /// Only the tracts that pass within `radius` of `center` (any point inside the
    /// sphere, ends included; points exactly `radius` away count). Errors for a
    /// negative or non-finite radius.
    pub fn through_sphere(&self, center: [f32; 3], radius: f64) -> Result<TractSet> {
        crate::numeric::ensure_finite("radius", radius)?;
        if radius < 0.0 {
            return Err(Error::InvalidParameter {
                name: "radius".into(),
                reason: format!("{radius} is negative"),
            });
        }
        Ok(self.filter(|_, t| t.points.iter().any(|&p| distance(p, center) <= radius)))
    }

    /// Only the tracts with a point inside the box `[min, max]` (ends included).
    pub fn through_box(&self, min: [f32; 3], max: [f32; 3]) -> TractSet {
        self.filter(|_, t| {
            t.points
                .iter()
                .any(|p| (0..3).all(|k| p[k] >= min[k] && p[k] <= max[k]))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tract(id: i32, points: &[[f32; 3]]) -> Tract {
        Tract::new(id, points.to_vec()).unwrap()
    }

    fn set() -> TractSet {
        TractSet {
            bundles: vec![
                TractBundle {
                    tag: Some(7),
                    alt_tag: Some(3),
                    ends: Some("A-B".into()),
                    tracts: vec![
                        tract(0, &[[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]]), // length 5
                        tract(1, &[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, 3.0]]), // 3
                    ],
                },
                TractBundle {
                    tag: Some(8),
                    alt_tag: None,
                    ends: None,
                    tracts: vec![tract(2, &[[10.0, 10.0, 10.0]])],
                },
            ],
        }
    }

    #[test]
    fn length_arc_length_and_tangents() {
        let s = set();
        let t = &s.bundles[0].tracts[1];
        assert_eq!(t.length(), 3.0);
        assert_eq!(t.arc_lengths(), vec![0.0, 1.0, 3.0]);
        assert_eq!(t.tangents(), vec![[0.0, 0.0, 1.0]; 3]);
        // The polyline length of AFNI's own example (see tract_conformance).
        let one = tract(9, &[[0.0, 0.0, 0.0], [1.0, 2.0, 0.25], [2.0, 4.0, 1.0]]);
        assert!((one.length() - 4.608495).abs() < 1e-5);
        // A single point has no length and no direction; a repeated point has none there.
        let single = &s.bundles[1].tracts[0];
        assert_eq!((single.length(), single.tangents()), (0.0, vec![[0.0; 3]]));
        let rep = tract(5, &[[0.0; 3], [0.0; 3]]);
        assert_eq!(rep.tangents(), vec![[0.0; 3]; 2]);
        // Central differences at an interior bend are the unit average direction.
        let bend = tract(6, &[[0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]);
        let tg = bend.tangents();
        let h = std::f32::consts::FRAC_1_SQRT_2;
        assert!((tg[1][0] - h).abs() < 1e-6 && (tg[1][1] - h).abs() < 1e-6);
    }

    #[test]
    fn counts_bounds_and_frame_flip() {
        let s = set();
        assert_eq!((s.tract_count(), s.point_count()), (3, 6));
        let b = s.bounds().unwrap();
        assert_eq!((b.min, b.max), ([0.0, 0.0, 0.0], [10.0, 10.0, 10.0]));
        assert_eq!((b.center, b.radius), ([5.0, 5.0, 5.0], 5.0));
        assert!(TractSet::default().bounds().is_none());
        let f = s.flipped();
        assert_eq!(f.bundles[0].tracts[0].points[1], [-3.0, -4.0, 0.0]);
        assert_eq!(f.flipped(), s, "the flip is its own inverse");
        // Lengths do not depend on the frame; bounds do.
        assert_eq!(f.bundles[0].tracts[0].length(), 5.0);
        assert_eq!(f.bounds().unwrap().min, [-10.0, -10.0, 0.0]);
    }

    #[test]
    fn selections() {
        let s = set();
        assert_eq!(s.with_bundle_tag(8).tract_count(), 1);
        assert_eq!(s.with_bundle_tag(99).tract_count(), 0);
        assert_eq!(s.with_ids(&[1, 2]).tract_count(), 2);
        // Lengths 5, 3, 0: keep 3 <= length <= 5.
        let mid = s
            .with_length(&Threshold::Between { lo: 3.0, hi: 5.0 })
            .unwrap();
        assert_eq!(
            mid.iter_tracts().map(|t| t.id).collect::<Vec<_>>(),
            vec![0, 1]
        );
        assert!(s.with_length(&Threshold::AbsoluteAbove(-1.0)).is_err());
        assert_eq!(s.with_min_points(2).tract_count(), 2);
        // A bundle emptied by a filter keeps its tag until asked to drop.
        let only_long = s.with_length(&Threshold::Above(4.0)).unwrap();
        assert_eq!(only_long.bundles.len(), 2);
        assert_eq!(only_long.bundles[1].tag, Some(8));
        assert_eq!(only_long.without_empty_bundles().bundles.len(), 1);
        // Spheres and boxes: tract 1 reaches z = 3; the radius is inclusive.
        let near = s.through_sphere([0.0, 0.0, 3.0], 0.0).unwrap();
        assert_eq!(
            near.iter_tracts().map(|t| t.id).collect::<Vec<_>>(),
            vec![1]
        );
        assert!(s.through_sphere([0.0; 3], -1.0).is_err());
        assert!(s.through_sphere([0.0; 3], f64::NAN).is_err());
        let boxed = s.through_box([9.0, 9.0, 9.0], [11.0, 11.0, 11.0]);
        assert_eq!(
            boxed.iter_tracts().map(|t| t.id).collect::<Vec<_>>(),
            vec![2]
        );
    }

    #[test]
    fn construction_checks() {
        assert!(Tract::new(0, vec![]).is_err());
        assert!(Tract::new(0, vec![[0.0, f32::NAN, 0.0]]).is_err());
        assert!(Tract::new(0, vec![[0.0, f32::INFINITY, 0.0]]).is_err());
    }
}
