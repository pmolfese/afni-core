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
// AFNI's built-in "big" color scales, reproduced exactly: the two hue-to-RGB
// functions `DC_spectrum_AJJ` and `DC_spectrum_ZSS`, and the nine 256-entry
// scales built from them (`NJ_bigmaps_init` in AFNI's display.c): the red-to-blue
// spectrum with and without a black gap, yellow-to-cyan with and without a gap,
// yellow-to-red, two color circles, and two "reds and blues" scales.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `color.rs` supplies `Rgba` and `ContinuousColorMap`; `AfniColorScale::to_color_map`
//   converts a table made here into one.
// * Phase 5 (overlays) takes a color scale, maps data to positions, and looks up
//   colors; this file only says what the scales ARE.
// * Golden tables generated from AFNI's own C code are committed in
//   tests/data/conformance/afni_colorscales.ref, and tests/colorscale_conformance.rs
//   compares every entry of every scale.
//
// ORIENTATION: AFNI stores a scale with index 0 at the TOP of the color bar (the
// highest data value) and the last index at the bottom. Tables returned by
// `AfniColorScale::table` keep AFNI's index order so they can be compared byte
// for byte; `to_color_map` flips them so that position 0 is the LOWEST value.
//
// EXACTNESS: AFNI computes channels as `(int)(255 * pow(x, gamma) + 0.5)` where
// its private `pow` is `exp(gamma * log(x))` (and 0 for x <= 0). That is not
// bit-identical to the math library's `pow`, and the difference can change a
// byte at a rounding boundary, so the same formula is used here.
// ---------------------------------------------------------------------------

//! AFNI's built-in color scales (`display.c`).

use crate::color::{ColorStop, ContinuousColorMap, Rgba};
use crate::error::{Error, Result};
use crate::overlay::ColorTable;

/// AFNI's private power function (`mypow` in display.c): 0 for `x <= 0`, the
/// identity for `y == 1`, otherwise `exp(y ln x)`.
fn afni_pow(x: f64, y: f64) -> f64 {
    if x <= 0.0 {
        0.0
    } else if y == 1.0 {
        x
    } else {
        (y * x.ln()).exp()
    }
}

/// One color channel from a unit-range value: `(int)(255 * x^gamma + 0.5)`.
/// The cast truncates toward zero, as in C.
fn channel(x: f64, gamma: f64) -> i32 {
    (255.0 * afni_pow(x, gamma) + 0.5) as i32
}

/// Clamp a computed channel into a byte. AFNI stores these in `byte` fields; the
/// formulas never leave `0..=255`, so this never changes a value in practice.
fn byte(value: i32) -> u8 {
    value.clamp(0, 255) as u8
}

/// Bring a hue angle into `[0, 360]` the way both AFNI spectrum functions do.
fn wrap_degrees(mut angle: f64) -> f64 {
    while angle < 0.0 {
        angle += 360.0;
    }
    while angle > 360.0 {
        angle -= 360.0;
    }
    angle
}

/// AFNI's `DC_spectrum_AJJ`: a hue angle in degrees (0 red, 60 yellow, 120
/// green, 180 cyan, 240 blue, 300 purple) to RGB bytes, after Andrzej
/// Jesmanowicz with RW Cox's constants (`ak = ab = 5`). A `gamma <= 0` means 1.
///
/// ```
/// use afni_core::afni_colors::spectrum_ajj;
/// assert_eq!(spectrum_ajj(60.0, 0.8), [255, 255, 0]); // yellow
/// ```
pub fn spectrum_ajj(angle: f64, gamma: f64) -> [u8; 3] {
    let gamma = if gamma <= 0.0 { 1.0 } else { gamma };
    // RWC's choices: the floor of each channel's ramp is 5 of 255.
    let (ak, ab) = (5.0_f64, 5.0_f64);
    let s = 255.0 - ak;
    let c = s / 60.0;
    let sb = 255.0 - ab;
    let cb = sb / 60.0;
    let an = wrap_degrees(angle);
    let (r, g, b);
    if an < 120.0 {
        r = channel((ak + s.min((120.0 - an) * c)) / 255.0, gamma);
        g = channel((ak + s.min(an * c)) / 255.0, gamma);
        b = 0;
    } else if an < 240.0 {
        r = 0;
        g = channel((ak + s.min((240.0 - an) * c)) / 255.0, gamma);
        b = channel((ab + sb.min((an - 120.0) * cb)) / 255.0, gamma);
    } else {
        r = channel((ak + s.min((an - 240.0) * c)) / 255.0, gamma);
        g = 0;
        // AFNI uses `s` (not `sb`) in this branch; both equal 250 here.
        b = channel((ab + s.min((360.0 - an) * cb)) / 255.0, gamma);
    }
    [byte(r), byte(g), byte(b)]
}

/// AFNI's `DC_spectrum_ZSS` (after Ziad Saad): hue angle in degrees to RGB
/// bytes through four quadrants. A `gamma <= 0` means 1.
pub fn spectrum_zss(angle: f64, gamma: f64) -> [u8; 3] {
    let gamma = if gamma <= 0.0 { 1.0 } else { gamma };
    let an = wrap_degrees(angle) / 90.0;
    let (r, g, b);
    if an <= 1.0 {
        r = channel(1.0 - an, gamma);
        g = channel(0.5 * an, gamma);
        b = channel(an, gamma);
    } else if an <= 2.0 {
        r = 0;
        g = channel(0.5 * an, gamma);
        b = channel(2.0 - an, gamma);
    } else if an <= 3.0 {
        r = channel(an - 2.0, gamma);
        g = 255;
        b = 0;
    } else {
        r = 255;
        g = channel(4.0 - an, gamma);
        b = 0;
    }
    [byte(r), byte(g), byte(b)]
}

// Hue fiducials from display.h.
const AJJ_YEL: f64 = 60.0;
const AJJ_CYN: f64 = 180.0;
const AJJ_BLU: f64 = 240.0;

/// AFNI's nine built-in 256-entry ("big") color scales, in AFNI's own order
/// (`BIGMAP_NAMES` in display.h).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AfniColorScale {
    /// `Spectrum:red_to_blue`: red at the top through yellow, green, cyan to blue.
    SpectrumRedToBlue,
    /// `Spectrum:red_to_blue+gap`: red-yellow, a black gap, then cyan-blue.
    SpectrumRedToBlueGap,
    /// `Spectrum:yellow_to_cyan`: yellow-red, a purple bridge, blue-cyan.
    SpectrumYellowToCyan,
    /// `Spectrum:yellow_to_cyan+gap`: the same with a black gap in the middle.
    SpectrumYellowToCyanGap,
    /// `Spectrum:yellow_to_red` (gamma 0.7).
    SpectrumYellowToRed,
    /// `Color_circle_AJJ`: a full hue circle.
    ColorCircleAjj,
    /// `Color_circle_ZSS`: a full hue circle from the ZSS function.
    ColorCircleZss,
    /// `Reds_and_Blues`: yellow-red over blue-cyan, no green.
    RedsAndBlues,
    /// `Reds_and_Blues_w_Green`: `RedsAndBlues` with two green entries in the middle.
    RedsAndBluesWithGreen,
}

impl AfniColorScale {
    /// Every scale, in AFNI's order.
    pub const ALL: [AfniColorScale; 9] = [
        AfniColorScale::SpectrumRedToBlue,
        AfniColorScale::SpectrumRedToBlueGap,
        AfniColorScale::SpectrumYellowToCyan,
        AfniColorScale::SpectrumYellowToCyanGap,
        AfniColorScale::SpectrumYellowToRed,
        AfniColorScale::ColorCircleAjj,
        AfniColorScale::ColorCircleZss,
        AfniColorScale::RedsAndBlues,
        AfniColorScale::RedsAndBluesWithGreen,
    ];

    /// AFNI's name for the scale (`BIGMAP_NAMES`).
    pub fn name(self) -> &'static str {
        match self {
            Self::SpectrumRedToBlue => "Spectrum:red_to_blue",
            Self::SpectrumRedToBlueGap => "Spectrum:red_to_blue+gap",
            Self::SpectrumYellowToCyan => "Spectrum:yellow_to_cyan",
            Self::SpectrumYellowToCyanGap => "Spectrum:yellow_to_cyan+gap",
            Self::SpectrumYellowToRed => "Spectrum:yellow_to_red",
            Self::ColorCircleAjj => "Color_circle_AJJ",
            Self::ColorCircleZss => "Color_circle_ZSS",
            Self::RedsAndBlues => "Reds_and_Blues",
            Self::RedsAndBluesWithGreen => "Reds_and_Blues_w_Green",
        }
    }

    /// Find a scale by AFNI's name (case-insensitive).
    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL
            .into_iter()
            .find(|s| s.name().eq_ignore_ascii_case(name.trim()))
    }

    /// The scale's colors as RGB bytes in AFNI's index order (index 0 is the
    /// TOP of the color bar, the highest value), for a scale of `entries`
    /// colors. AFNI's default is 256 (`NPANE_BIG`); it can be changed with an
    /// environment variable, and every formula scales with it.
    ///
    /// `entries` must be at least 32 (below that the black "gap" is empty) and at
    /// most 2048 (AFNI's `NPANE_BIGGEST`).
    pub fn table(self, entries: usize) -> Result<Vec<[u8; 3]>> {
        if !(32..=2048).contains(&entries) {
            return Err(Error::InvalidParameter {
                name: "entries".into(),
                reason: format!("AFNI color scales have 32 to 2048 entries, got {entries}"),
            });
        }
        let n = entries as i64;
        // Integer geometry from display.h; `/` is integer division, as in C.
        let gap = n / 32; // NBIG_GAP
        let mbot = n / 2 - gap; // NBIG_MBOT
        let mtop = n / 2 + gap; // NBIG_MTOP
        let half = n / 2; // NPB_2
        let nf = n as f64;
        let table = (0..n)
            .map(|ii| {
                let i = ii as f64;
                match self {
                    Self::SpectrumRedToBlue => {
                        spectrum_ajj(i * ((AJJ_BLU + 8.0) / (nf - 1.0)) - 4.0, 0.8)
                    }
                    Self::SpectrumYellowToRed => {
                        spectrum_ajj(60.0 - i * (AJJ_YEL / (nf - 1.0)), 0.7)
                    }
                    Self::ColorCircleAjj => spectrum_ajj(i * (360.0 / (nf - 1.0)), 0.8),
                    Self::ColorCircleZss => spectrum_zss(360.0 - i * (360.0 / (nf - 1.0)), 1.0),
                    Self::SpectrumRedToBlueGap => {
                        if ii < mbot {
                            spectrum_ajj(i * (AJJ_YEL / (mbot as f64 - 1.0)), 0.8)
                        } else if ii > mtop {
                            spectrum_ajj(
                                AJJ_CYN
                                    + (ii - mtop - 1) as f64 * (60.0 / (nf - mtop as f64 - 2.0)),
                                0.8,
                            )
                        } else {
                            [0, 0, 0]
                        }
                    }
                    Self::SpectrumYellowToCyan | Self::SpectrumYellowToCyanGap => {
                        if ii < mbot {
                            spectrum_ajj(AJJ_YEL - i * (AJJ_YEL / (mbot as f64 - 1.0)), 0.8)
                        } else if ii > mtop {
                            spectrum_ajj(
                                AJJ_BLU
                                    - (ii - mtop - 1) as f64 * (60.0 / (nf - mtop as f64 - 2.0)),
                                0.8,
                            )
                        } else if self == Self::SpectrumYellowToCyanGap {
                            [0, 0, 0]
                        } else {
                            spectrum_ajj(
                                360.0
                                    - (ii - mbot + 1) as f64
                                        * (120.0 / ((mtop - mbot) as f64 + 2.0)),
                                0.8,
                            )
                        }
                    }
                    Self::RedsAndBlues | Self::RedsAndBluesWithGreen => {
                        if ii < half {
                            spectrum_ajj(AJJ_YEL - i * (AJJ_YEL / (half as f64 - 1.0)), 0.8)
                        } else {
                            // AFNI subtracts MTOP + 1 (not the half-way point) here, so
                            // the first upper entries use a NEGATIVE offset, which makes the
                            // hue start above 240 degrees. Kept: it is AFNI's scale.
                            spectrum_ajj(
                                AJJ_BLU
                                    - (ii - mtop - 1) as f64 * (60.0 / (nf - half as f64 - 2.0)),
                                0.8,
                            )
                        }
                    }
                }
            })
            .collect::<Vec<_>>();
        let mut table = table;
        if self == Self::RedsAndBluesWithGreen {
            // Two green entries at the seam (indices half-1 and half).
            let green = spectrum_ajj(half as f64 * ((AJJ_BLU + 8.0) / (nf - 1.0)) - 4.0, 0.8);
            table[(half - 1) as usize] = green;
            table[half as usize] = green;
        }
        Ok(table)
    }

    /// The scale as a [`ColorTable`] for pane-style lookup: entry 0 is the LOWEST
    /// value (AFNI's last index) and the last entry the highest (AFNI's index 0).
    pub fn to_color_table(self, entries: usize) -> Result<ColorTable> {
        let mut table = self.table(entries)?;
        table.reverse(); // AFNI's index 0 is the top of the bar
        ColorTable::new(
            table
                .iter()
                .map(|&[r, g, b]| Rgba::from_u8(r, g, b, 255))
                .collect(),
        )
    }

    /// The scale as a [`ContinuousColorMap`] with position 0 at the LOWEST value
    /// (AFNI's index `entries - 1`) and position 1 at the highest (index 0).
    /// One stop per table entry, colors exactly the table's bytes.
    pub fn to_color_map(self, entries: usize) -> Result<ContinuousColorMap> {
        let table = self.table(entries)?;
        let last = (table.len() - 1) as f64;
        let stops = table
            .iter()
            .enumerate()
            .rev() // lowest value first
            .map(|(index, &[r, g, b])| {
                ColorStop::new((last - index as f64) / last, Rgba::from_u8(r, g, b, 255))
            })
            .collect::<Result<Vec<_>>>()?;
        ContinuousColorMap::new(self.name(), stops)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spectrum_fiducials() {
        // Pure hues: the ramps top out at 255 for gamma 1 and a full-strength hue.
        assert_eq!(spectrum_ajj(120.0, 1.0), [0, 255, 5]); // green; blue keeps the 5/255 floor
        assert_eq!(spectrum_ajj(60.0, 0.8), [255, 255, 0]);
        assert_eq!(spectrum_ajj(0.0, 1.0), [255, 5, 0]);
        assert_eq!(
            spectrum_ajj(-30.0, 0.8),
            spectrum_ajj(330.0, 0.8),
            "negative hues wrap"
        );
        assert_eq!(
            spectrum_ajj(390.0, 0.8),
            spectrum_ajj(30.0, 0.8),
            "large hues wrap"
        );
        assert_eq!(
            spectrum_ajj(100.0, 0.0),
            spectrum_ajj(100.0, 1.0),
            "gamma <= 0 means 1"
        );
        assert_eq!(spectrum_zss(0.0, 1.0), [255, 0, 0]);
        assert_eq!(spectrum_zss(180.0, 1.0), [0, 255, 0]);
        assert_eq!(spectrum_zss(270.0, 1.0), [255, 255, 0]);
    }

    #[test]
    fn names_round_trip_and_order_matches_afni() {
        let names: Vec<_> = AfniColorScale::ALL.iter().map(|s| s.name()).collect();
        assert_eq!(names[0], "Spectrum:red_to_blue");
        assert_eq!(names[8], "Reds_and_Blues_w_Green");
        for s in AfniColorScale::ALL {
            assert_eq!(AfniColorScale::from_name(s.name()), Some(s));
            assert_eq!(AfniColorScale::from_name(&s.name().to_uppercase()), Some(s));
        }
        assert_eq!(AfniColorScale::from_name("nope"), None);
    }

    #[test]
    fn tables_have_the_requested_size_and_validate_it() {
        for s in AfniColorScale::ALL {
            assert_eq!(s.table(256).unwrap().len(), 256);
            assert_eq!(s.table(128).unwrap().len(), 128);
        }
        assert!(AfniColorScale::ColorCircleAjj.table(31).is_err());
        assert!(AfniColorScale::ColorCircleAjj.table(2049).is_err());
    }

    #[test]
    fn gapped_scales_have_a_black_band_and_others_do_not() {
        let gap = AfniColorScale::SpectrumRedToBlueGap.table(256).unwrap();
        // NBIG_GAP = 8, so indices 120..=136 are black.
        assert!(gap[120..=136].iter().all(|c| *c == [0, 0, 0]));
        assert_ne!(gap[119], [0, 0, 0]);
        assert_ne!(gap[137], [0, 0, 0]);
        let plain = AfniColorScale::SpectrumYellowToCyan.table(256).unwrap();
        assert!(plain[120..=136].iter().all(|c| *c != [0, 0, 0]));
        let green = AfniColorScale::RedsAndBluesWithGreen.table(256).unwrap();
        let base = AfniColorScale::RedsAndBlues.table(256).unwrap();
        assert_eq!(green[127], green[128]);
        assert_ne!(green[127], base[127]);
        let differing = (0..256).filter(|&i| green[i] != base[i]).count();
        assert_eq!(differing, 2, "only the two seam entries differ");
    }

    #[test]
    fn color_map_is_oriented_low_to_high_with_exact_bytes() {
        let table = AfniColorScale::SpectrumRedToBlue.table(256).unwrap();
        let map = AfniColorScale::SpectrumRedToBlue.to_color_map(256).unwrap();
        assert_eq!(map.stops().len(), 256);
        // AFNI index 0 (red) is the TOP: position 1.0. The last index is position 0.
        let top = map.stops().last().unwrap();
        assert_eq!(top.position, 1.0);
        assert_eq!(
            top.color.to_u8(),
            [table[0][0], table[0][1], table[0][2], 255]
        );
        let bottom = map.stops().first().unwrap();
        assert_eq!(bottom.position, 0.0);
        let last = table[255];
        assert_eq!(bottom.color.to_u8(), [last[0], last[1], last[2], 255]);
        // Quantization is lossless: every stop converts back to the exact table byte.
        for (k, stop) in map.stops().iter().enumerate() {
            let t = table[255 - k];
            assert_eq!(stop.color.to_u8()[..3], t, "stop {k}");
        }
    }
}
