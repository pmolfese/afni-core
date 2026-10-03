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
// The color model: `Rgba` (a color), `ColorStop`, `ContinuousColorMap` (colors
// at positions 0..1 with explicit interpolation), and `ColorMap`, which is
// either a continuous map or a label-color map. Everything is plain data with
// pure functions; no GUI or GPU types.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `afni_colors.rs` builds AFNI's built-in scales as `ContinuousColorMap`s.
// * `labels.rs` supplies `LabelColorMap` (integer keys to colors) for the
//   `ColorMap::Labels` variant.
// * `column.rs` supplies `ColumnRange`, used to turn a data value into a position.
// * Phase 5 (overlays) will combine a column, a threshold, and one of these maps
//   into a buffer of colors; `lookup_table` produces the table a GPU shader
//   would sample, keeping shader code itself in the viewer.
//
// CONVENTIONS (decisions, documented because they affect every pixel)
//
// * Channels are `f32` in 0..=1, straight (non-premultiplied) alpha.
// * Interpolation defaults to ENCODED RGB (the numbers as stored), which is
//   what AFNI and SUMA do and so what parity requires. Interpolating in linear
//   light is available but is a separate, explicitly named option.
// * Positions are 0..=1, 0 = lowest data value. A position outside 0..=1 is
//   clamped; a NaN position is an error (or the caller's "missing" color).
// * DUPLICATE STOPS make a hard edge. At exactly the shared position the EARLIER
//   stop's color is returned; just above it, the later stop's.
// ---------------------------------------------------------------------------

//! Colors, color stops, and continuous color maps.

use crate::column::ColumnRange;
use crate::error::{Error, Result};
use crate::labels::LabelColorMap;
use crate::numeric::ensure_finite;

/// A color with `f32` channels in `[0, 1]` and straight (non-premultiplied)
/// alpha.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rgba {
    /// Red, 0 to 1.
    pub r: f32,
    /// Green, 0 to 1.
    pub g: f32,
    /// Blue, 0 to 1.
    pub b: f32,
    /// Opacity, 0 (transparent) to 1 (opaque).
    pub a: f32,
}

impl Rgba {
    /// Fully transparent black.
    pub const TRANSPARENT: Rgba = Rgba {
        r: 0.0,
        g: 0.0,
        b: 0.0,
        a: 0.0,
    };
    /// Opaque black.
    pub const BLACK: Rgba = Rgba {
        r: 0.0,
        g: 0.0,
        b: 0.0,
        a: 1.0,
    };
    /// Opaque white.
    pub const WHITE: Rgba = Rgba {
        r: 1.0,
        g: 1.0,
        b: 1.0,
        a: 1.0,
    };

    /// A validated color: every channel finite and within `[0, 1]`.
    pub fn new(r: f32, g: f32, b: f32, a: f32) -> Result<Self> {
        for (name, v) in [("red", r), ("green", g), ("blue", b), ("alpha", a)] {
            ensure_finite(&format!("{name} channel"), f64::from(v))?;
            if !(0.0..=1.0).contains(&v) {
                return Err(Error::InvalidParameter {
                    name: format!("{name} channel"),
                    reason: format!("{v} is outside [0, 1]"),
                });
            }
        }
        Ok(Self { r, g, b, a })
    }

    /// From 8-bit channels, exactly `value / 255`.
    pub fn from_u8(r: u8, g: u8, b: u8, a: u8) -> Self {
        let f = |v: u8| f32::from(v) / 255.0;
        Self {
            r: f(r),
            g: f(g),
            b: f(b),
            a: f(a),
        }
    }

    /// To 8-bit channels, rounding to the nearest. `from_u8` followed by `to_u8`
    /// is always the identity.
    pub fn to_u8(self) -> [u8; 4] {
        let q = |v: f32| (v.clamp(0.0, 1.0) * 255.0).round() as u8;
        [q(self.r), q(self.g), q(self.b), q(self.a)]
    }

    /// The channels as `[r, g, b, a]`, the order a GPU texture expects.
    pub fn to_array(self) -> [f32; 4] {
        [self.r, self.g, self.b, self.a]
    }

    /// The same color with opacity `a` (clamped to `[0, 1]`; NaN becomes 0).
    pub fn with_alpha(self, a: f32) -> Self {
        Self {
            a: if a.is_nan() { 0.0 } else { a.clamp(0.0, 1.0) },
            ..self
        }
    }

    /// Brightness scaling: multiply the color channels by `factor` (clamped to
    /// `[0, 1]`), leaving alpha alone. `factor = 1` is unchanged; `0` is black.
    /// A NaN factor is treated as 0, i.e. black: callers should decide what a
    /// missing brightness means before getting here.
    pub fn scaled(self, factor: f32) -> Self {
        let f = if factor.is_nan() {
            0.0
        } else {
            factor.max(0.0)
        };
        let c = |v: f32| (v * f).clamp(0.0, 1.0);
        Self {
            r: c(self.r),
            g: c(self.g),
            b: c(self.b),
            a: self.a,
        }
    }

    /// Linear interpolation of all four channels as stored (encoded RGB), with
    /// `t` clamped to `[0, 1]`.
    pub fn lerp(self, other: Rgba, t: f32) -> Rgba {
        let t = t.clamp(0.0, 1.0);
        let mix = |a: f32, b: f32| a + (b - a) * t;
        Rgba {
            r: mix(self.r, other.r),
            g: mix(self.g, other.g),
            b: mix(self.b, other.b),
            a: mix(self.a, other.a),
        }
    }

    /// Interpolation in LINEAR LIGHT: the sRGB-encoded color channels are
    /// converted to linear intensity, mixed, and converted back. Alpha mixes
    /// linearly. Differs visibly from [`lerp`](Self::lerp) (midpoints are
    /// brighter); AFNI and SUMA do not do this.
    pub fn lerp_linear_light(self, other: Rgba, t: f32) -> Rgba {
        let t = t.clamp(0.0, 1.0);
        let mix = |a: f32, b: f32| {
            let (la, lb) = (srgb_to_linear(a), srgb_to_linear(b));
            linear_to_srgb(la + (lb - la) * t)
        };
        Rgba {
            r: mix(self.r, other.r),
            g: mix(self.g, other.g),
            b: mix(self.b, other.b),
            a: self.a + (other.a - self.a) * t,
        }
    }
}

/// The sRGB transfer function's decoding half: encoded `[0, 1]` to linear light.
fn srgb_to_linear(c: f32) -> f32 {
    if c <= 0.04045 {
        c / 12.92
    } else {
        ((c + 0.055) / 1.055).powf(2.4)
    }
}

/// The sRGB encoding half: linear light to encoded `[0, 1]`.
fn linear_to_srgb(l: f32) -> f32 {
    if l <= 0.003_130_8 {
        l * 12.92
    } else {
        1.055 * l.powf(1.0 / 2.4) - 0.055
    }
}

/// A color at a position in `[0, 1]` along a color map.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ColorStop {
    /// Position, 0 (lowest value) to 1 (highest).
    pub position: f64,
    /// The color at this position.
    pub color: Rgba,
}

impl ColorStop {
    /// A stop; `position` must be finite and within `[0, 1]`.
    pub fn new(position: f64, color: Rgba) -> Result<Self> {
        ensure_finite("color stop position", position)?;
        if !(0.0..=1.0).contains(&position) {
            return Err(Error::InvalidParameter {
                name: "color stop position".into(),
                reason: format!("{position} is outside [0, 1]"),
            });
        }
        Ok(Self { position, color })
    }
}

/// How colors between two stops are chosen.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Interpolation {
    /// Mix the two neighbouring stops linearly (the default).
    #[default]
    Continuous,
    /// "Direct"/banded: the color of the nearest stop AT OR BELOW the position
    /// (the first stop's color below the first position). Positions exactly on a
    /// stop get that stop. This is how an N-pane color bar behaves.
    Stepped,
    /// The color of the stop closest to the position; ties go to the lower one.
    Nearest,
}

/// The color space in which [`Interpolation::Continuous`] mixes colors.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum InterpolationSpace {
    /// The channel values as stored, i.e. gamma-encoded RGB. What AFNI and SUMA
    /// do, and the default for parity.
    #[default]
    EncodedRgb,
    /// Linear light (a perceptually motivated option, not AFNI's behavior).
    LinearLight,
}

/// A named continuous color map: stops sorted by position, plus how to
/// interpolate between them.
#[derive(Debug, Clone, PartialEq)]
pub struct ContinuousColorMap {
    name: String,
    stops: Vec<ColorStop>,
    interpolation: Interpolation,
    space: InterpolationSpace,
}

impl ContinuousColorMap {
    /// Build a map. The name must not be blank, there must be at least one
    /// stop, and positions must be in ascending order. Equal positions are
    /// allowed and make a hard edge (see the module docs).
    pub fn new(name: impl Into<String>, stops: Vec<ColorStop>) -> Result<Self> {
        let name = name.into();
        if name.trim().is_empty() {
            return Err(Error::Empty("color map name".into()));
        }
        if stops.is_empty() {
            return Err(Error::Empty("color map stops".into()));
        }
        for pair in stops.windows(2) {
            if pair[0].position > pair[1].position {
                return Err(Error::InvalidParameter {
                    name: "color stops".into(),
                    reason: "stops must be sorted by position".into(),
                });
            }
        }
        Ok(Self {
            name,
            stops,
            interpolation: Interpolation::Continuous,
            space: InterpolationSpace::EncodedRgb,
        })
    }

    /// A two-stop black-to-white ramp. Not an AFNI-defined map.
    pub fn grayscale() -> Self {
        Self::new(
            "Grayscale",
            vec![
                ColorStop {
                    position: 0.0,
                    color: Rgba::BLACK,
                },
                ColorStop {
                    position: 1.0,
                    color: Rgba::WHITE,
                },
            ],
        )
        .expect("static stops are valid")
    }

    /// Choose how colors between stops are picked.
    pub fn with_interpolation(mut self, interpolation: Interpolation) -> Self {
        self.interpolation = interpolation;
        self
    }

    /// Choose the color space for continuous interpolation.
    pub fn with_space(mut self, space: InterpolationSpace) -> Self {
        self.space = space;
        self
    }

    /// The map's name.
    pub fn name(&self) -> &str {
        &self.name
    }
    /// The stops, ascending by position.
    pub fn stops(&self) -> &[ColorStop] {
        &self.stops
    }
    /// The interpolation mode.
    pub fn interpolation(&self) -> Interpolation {
        self.interpolation
    }
    /// The interpolation color space.
    pub fn space(&self) -> InterpolationSpace {
        self.space
    }

    /// The color at `position`. Positions outside `[0, 1]` are clamped to the
    /// ends. A NaN position gives the first stop's color: use
    /// [`sample_checked`](Self::sample_checked) or [`sample_or`](Self::sample_or)
    /// where missing data must be handled explicitly.
    pub fn sample(&self, position: f64) -> Rgba {
        let first = self.stops[0];
        let last = self.stops[self.stops.len() - 1];
        if position.is_nan() || position <= first.position {
            return first.color;
        }
        if position >= last.position {
            return last.color;
        }
        match self.interpolation {
            Interpolation::Stepped => {
                // Last stop whose position is <= `position`.
                let at_or_below = self.stops.partition_point(|s| s.position <= position);
                self.stops[at_or_below - 1].color
            }
            Interpolation::Nearest => {
                let above = self.stops.partition_point(|s| s.position < position);
                let (lo, hi) = (self.stops[above - 1], self.stops[above]);
                // Ties go to the lower stop.
                if position - lo.position <= hi.position - position {
                    lo.color
                } else {
                    hi.color
                }
            }
            Interpolation::Continuous => {
                // First window whose right end is at or above `position`; at an
                // exact duplicate position this returns the EARLIER stop.
                let right = self.stops.partition_point(|s| s.position < position);
                let (lo, hi) = (self.stops[right - 1], self.stops[right]);
                let span = hi.position - lo.position;
                let t = if span <= 0.0 {
                    1.0
                } else {
                    ((position - lo.position) / span) as f32
                };
                match self.space {
                    InterpolationSpace::EncodedRgb => lo.color.lerp(hi.color, t),
                    InterpolationSpace::LinearLight => lo.color.lerp_linear_light(hi.color, t),
                }
            }
        }
    }

    /// Like [`sample`](Self::sample), but a non-finite position is an error.
    pub fn sample_checked(&self, position: f64) -> Result<Rgba> {
        ensure_finite("color map position", position)?;
        Ok(self.sample(position))
    }

    /// The color at `position`, or `missing` if it is not finite.
    pub fn sample_or(&self, position: f64, missing: Rgba) -> Rgba {
        if position.is_finite() {
            self.sample(position)
        } else {
            missing
        }
    }

    /// The color for a data `value` given the display `range`: the value is
    /// normalized into a position (clamped; a zero-width range maps to the
    /// middle) and sampled. A non-finite value gives `missing`.
    pub fn color_for_value(&self, value: f64, range: &ColumnRange, missing: Rgba) -> Rgba {
        if !value.is_finite() {
            return missing;
        }
        self.sample(range.normalized(value))
    }

    /// `entries` colors sampled evenly from position 0 to 1 (inclusive): the
    /// table a viewer uploads for GPU lookup. Needs at least 2 entries.
    pub fn lookup_table(&self, entries: usize) -> Result<Vec<Rgba>> {
        if entries < 2 {
            return Err(Error::InvalidParameter {
                name: "entries".into(),
                reason: format!("a lookup table needs at least 2 entries, got {entries}"),
            });
        }
        let last = (entries - 1) as f64;
        Ok((0..entries).map(|i| self.sample(i as f64 / last)).collect())
    }
}

/// How data are turned into colors: along a continuous map, or by integer label.
#[derive(Debug, Clone, PartialEq)]
pub enum ColorMap {
    /// Colors by position along a gradient.
    Continuous(ContinuousColorMap),
    /// Colors by integer label key.
    Labels(LabelColorMap),
}

impl ColorMap {
    /// The continuous map, if this is one.
    pub fn as_continuous(&self) -> Option<&ContinuousColorMap> {
        match self {
            ColorMap::Continuous(c) => Some(c),
            ColorMap::Labels(_) => None,
        }
    }

    /// The label-color map, if this is one.
    pub fn as_labels(&self) -> Option<&LabelColorMap> {
        match self {
            ColorMap::Labels(l) => Some(l),
            ColorMap::Continuous(_) => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gray(v: f32) -> Rgba {
        Rgba::new(v, v, v, 1.0).unwrap()
    }

    fn map(stops: &[(f64, Rgba)]) -> ContinuousColorMap {
        ContinuousColorMap::new(
            "t",
            stops
                .iter()
                .map(|&(p, c)| ColorStop::new(p, c).unwrap())
                .collect(),
        )
        .unwrap()
    }

    #[test]
    fn rgba_validation_and_byte_round_trip() {
        assert!(Rgba::new(0.0, 0.5, 1.0, 1.0).is_ok());
        assert!(Rgba::new(1.1, 0.0, 0.0, 1.0).is_err());
        assert!(Rgba::new(f32::NAN, 0.0, 0.0, 1.0).is_err());
        assert!(Rgba::new(0.0, 0.0, 0.0, -0.1).is_err());
        for v in 0..=255_u8 {
            assert_eq!(Rgba::from_u8(v, v, v, v).to_u8(), [v, v, v, v]);
        }
        assert_eq!(
            Rgba::from_u8(255, 0, 0, 255).to_array(),
            [1.0, 0.0, 0.0, 1.0]
        );
    }

    #[test]
    fn alpha_and_brightness_helpers() {
        let c = Rgba::new(0.8, 0.4, 0.2, 0.5).unwrap();
        assert_eq!(c.with_alpha(2.0).a, 1.0);
        assert_eq!(c.with_alpha(f32::NAN).a, 0.0);
        let dim = c.scaled(0.5);
        assert_eq!((dim.r, dim.g, dim.b, dim.a), (0.4, 0.2, 0.1, 0.5));
        assert_eq!(c.scaled(10.0).r, 1.0, "clamped, not wrapped");
        assert_eq!(c.scaled(1.0), c);
        assert_eq!(c.scaled(0.0).to_u8()[..3], [0, 0, 0]);
        assert_eq!(c.scaled(f32::NAN).to_u8()[..3], [0, 0, 0]);
    }

    #[test]
    fn construction_rules() {
        assert!(ContinuousColorMap::new(
            " ",
            vec![ColorStop {
                position: 0.0,
                color: Rgba::BLACK
            }]
        )
        .is_err());
        assert!(ContinuousColorMap::new("x", vec![]).is_err());
        let unsorted = vec![
            ColorStop::new(1.0, Rgba::BLACK).unwrap(),
            ColorStop::new(0.0, Rgba::BLACK).unwrap(),
        ];
        assert!(ContinuousColorMap::new("x", unsorted).is_err());
        assert!(ColorStop::new(1.5, Rgba::BLACK).is_err());
        assert!(ColorStop::new(f64::NAN, Rgba::BLACK).is_err());
    }

    #[test]
    fn continuous_interpolates_and_clamps() {
        let m = ContinuousColorMap::grayscale();
        assert_eq!(m.sample(0.25), gray(0.25));
        assert_eq!(m.sample(-3.0), gray(0.0));
        assert_eq!(m.sample(7.0), gray(1.0));
        assert_eq!(m.sample(0.0), Rgba::BLACK);
        assert_eq!(m.sample(1.0), Rgba::WHITE);
        // A single stop is a constant color.
        let one = map(&[(0.3, gray(0.4))]);
        assert_eq!(one.sample(0.0), gray(0.4));
        assert_eq!(one.sample(1.0), gray(0.4));
    }

    #[test]
    fn duplicate_stops_make_a_hard_edge_with_the_earlier_color_on_the_edge() {
        let red = Rgba::from_u8(255, 0, 0, 255);
        let blue = Rgba::from_u8(0, 0, 255, 255);
        let m = map(&[(0.0, red), (0.5, red), (0.5, blue), (1.0, blue)]);
        assert_eq!(m.sample(0.25), red);
        assert_eq!(m.sample(0.5), red, "exactly on the edge: the earlier stop");
        assert_eq!(m.sample(0.5000001), blue);
        assert_eq!(m.sample(0.75), blue);
        // Duplicates at the very ends.
        let ends = map(&[(0.0, red), (0.0, blue), (1.0, blue)]);
        assert_eq!(ends.sample(0.0), red);
        assert_eq!(ends.sample(0.1), blue);
    }

    #[test]
    fn stepped_and_nearest_modes() {
        let (a, b, c) = (gray(0.0), gray(0.5), gray(1.0));
        let stops = [(0.0, a), (0.4, b), (1.0, c)];
        let stepped = map(&stops).with_interpolation(Interpolation::Stepped);
        assert_eq!(stepped.sample(0.0), a);
        assert_eq!(stepped.sample(0.39), a);
        assert_eq!(
            stepped.sample(0.4),
            b,
            "a position exactly on a stop takes that stop"
        );
        assert_eq!(stepped.sample(0.99), b);
        assert_eq!(stepped.sample(1.0), c);
        let nearest = map(&stops).with_interpolation(Interpolation::Nearest);
        assert_eq!(nearest.sample(0.19), a);
        assert_eq!(nearest.sample(0.2), a, "a tie goes to the lower stop");
        assert_eq!(nearest.sample(0.21), b);
        assert_eq!(nearest.sample(0.69), b);
        assert_eq!(nearest.sample(0.71), c);
    }

    #[test]
    fn interpolation_space_is_explicit_and_encoded_is_the_default() {
        let m = ContinuousColorMap::grayscale();
        assert_eq!(m.space(), InterpolationSpace::EncodedRgb);
        let linear = m.clone().with_space(InterpolationSpace::LinearLight);
        let encoded_mid = m.sample(0.5).r;
        let linear_mid = linear.sample(0.5).r;
        assert_eq!(encoded_mid, 0.5);
        // Mixing in linear light gives the familiar brighter midpoint (~0.735).
        assert!((linear_mid - 0.735).abs() < 0.005, "{linear_mid}");
        // Ends are unaffected, and the transfer functions invert each other.
        assert_eq!(linear.sample(0.0), Rgba::BLACK);
        assert_eq!(linear.sample(1.0).to_u8(), [255, 255, 255, 255]);
        for v in [0.0_f32, 0.02, 0.1, 0.5, 0.9, 1.0] {
            assert!((linear_to_srgb(srgb_to_linear(v)) - v).abs() < 1e-6, "{v}");
        }
    }

    #[test]
    fn missing_and_nonfinite_values() {
        let m = ContinuousColorMap::grayscale();
        let missing = Rgba::TRANSPARENT;
        let range = ColumnRange::new(0.0, 10.0).unwrap();
        assert_eq!(m.color_for_value(5.0, &range, missing), gray(0.5));
        assert_eq!(m.color_for_value(f64::NAN, &range, missing), missing);
        assert_eq!(m.color_for_value(f64::INFINITY, &range, missing), missing);
        assert_eq!(
            m.color_for_value(99.0, &range, missing),
            Rgba::WHITE,
            "clamped"
        );
        // Zero-width range maps to the middle of the map.
        let flat = ColumnRange::new(3.0, 3.0).unwrap();
        assert_eq!(m.color_for_value(3.0, &flat, missing), gray(0.5));
        assert!(m.sample_checked(f64::NAN).is_err());
        assert_eq!(m.sample_or(f64::NAN, missing), missing);
        assert_eq!(m.sample_or(0.5, missing), gray(0.5));
    }

    #[test]
    fn lookup_tables_are_evenly_sampled_and_validated() {
        let t = ContinuousColorMap::grayscale().lookup_table(5).unwrap();
        assert_eq!(t.len(), 5);
        assert_eq!(t[0], Rgba::BLACK);
        assert_eq!(t[2], gray(0.5));
        assert_eq!(t[4], Rgba::WHITE);
        assert!(ContinuousColorMap::grayscale().lookup_table(1).is_err());
    }

    #[test]
    fn color_map_enum_accessors() {
        let c = ColorMap::Continuous(ContinuousColorMap::grayscale());
        assert!(c.as_continuous().is_some() && c.as_labels().is_none());
    }
}
