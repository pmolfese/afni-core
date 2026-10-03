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
// SUMA's own named color maps (`rgybr20`, `bgyr19`, `gray20`, `byr64`, ...), the
// ones listed by `SUMA_MakeStandardMap` in SUMA/SUMA_Color.c. They are different
// from AFNI's nine 256-entry "big" scales in `afni_colors.rs`: SUMA builds each
// one by linear interpolation between a few "fiducial" colors, giving a small
// table (2 to 64 entries) that SUMA then maps values onto in panes.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// Each map produces a `ColorTable` (see `overlay.rs`), which is exactly what the
// pane rules (`PaneRule::SumaInterpolated` and friends) consume. Entry 0 is the
// LOWEST value, matching SUMA (the opposite of AFNI's big scales, whose index 0 is
// the top of the bar).
//
// EXACTNESS
//
// SUMA computes with C `float`, and so do we: each entry is
// `fiducial + j * ((next - fiducial) / gap)` in `f32`, the same expression SUMA
// evaluates. `tests/suma_colormaps_conformance.rs` compares every map with AFNI's
// `MakeColorMap -std <name>`, which prints two decimals.
//
// The maps are SUMA's, not sumaru's. sumaru's own `fire`, `afni_p2_spanned` and
// `amber_monochrome` are separate (see the roadmap, R11).
// ---------------------------------------------------------------------------

//! SUMA's standard named color maps (`rgybr20`, `bgyr19`, `byr64`, ...).

use crate::color::Rgba;
use crate::error::Result;
use crate::overlay::ColorTable;

/// One of SUMA's standard color maps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SumaColorMap {
    /// `rgybr20`: red, green, blue, yellow, back toward red (20 colors; the last
    /// red is left out so the map wraps).
    Rgybr20,
    /// `bgyr19`: blue, green, yellow, red (19 colors).
    Bgyr19,
    /// `gray02`: two grays, 0.4 then 0.7.
    Gray02,
    /// `gray_i02`: `gray02` reversed, 0.7 then 0.4.
    GrayI02,
    /// `gray20`: 20 grays from 0.3 to 0.8.
    Gray20,
    /// `ngray20`: SUMA defines it identically to `gray20` (the "n" is historical).
    NGray20,
    /// `bw20`: 20 grays from black to white.
    Bw20,
    /// `byr64` (also `matlab_default_byr64`): MATLAB's old default "jet" map.
    Byr64,
    /// `bgyr64`: like `byr64` with a green section (a "rainbow").
    Bgyr64,
}

impl SumaColorMap {
    /// Every map, in the order SUMA's source lists them.
    pub const ALL: [SumaColorMap; 9] = [
        SumaColorMap::Rgybr20,
        SumaColorMap::Bgyr19,
        SumaColorMap::Gray02,
        SumaColorMap::GrayI02,
        SumaColorMap::Gray20,
        SumaColorMap::NGray20,
        SumaColorMap::Bw20,
        SumaColorMap::Byr64,
        SumaColorMap::Bgyr64,
    ];

    /// SUMA's name for the map.
    pub fn name(self) -> &'static str {
        match self {
            Self::Rgybr20 => "rgybr20",
            Self::Bgyr19 => "bgyr19",
            Self::Gray02 => "gray02",
            Self::GrayI02 => "gray_i02",
            Self::Gray20 => "gray20",
            Self::NGray20 => "ngray20",
            Self::Bw20 => "bw20",
            Self::Byr64 => "byr64",
            Self::Bgyr64 => "bgyr64",
        }
    }

    /// Find a map by SUMA's name (case-insensitive). `matlab_default_byr64` is an
    /// alias of `byr64`, as in SUMA.
    pub fn from_name(name: &str) -> Option<Self> {
        let name = name.trim();
        if name.eq_ignore_ascii_case("matlab_default_byr64") {
            return Some(Self::Byr64);
        }
        Self::ALL
            .into_iter()
            .find(|m| m.name().eq_ignore_ascii_case(name))
    }

    /// The map's colors as `[r, g, b]` in `0..=1`, lowest value first.
    pub fn colors(self) -> Vec<[f32; 3]> {
        const GRAY: fn(f32) -> [f32; 3] = |g| [g, g, g];
        match self {
            // `SUMA_MakeColorMap(fiducials, n, rgba=0, n_colors, skip_last, name)`:
            // equal gaps between fiducials.
            Self::Rgybr20 => equal_gaps(
                &[
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [1.0, 1.0, 0.0],
                    [1.0, 0.0, 0.0],
                ],
                20,
                true,
            ),
            Self::Bgyr19 => equal_gaps(
                &[
                    [0.0, 0.0, 1.0],
                    [0.0, 1.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [1.0, 0.0, 0.0],
                ],
                19,
                false,
            ),
            Self::Gray02 => equal_gaps(&[GRAY(0.4), GRAY(0.7)], 2, false),
            Self::GrayI02 => equal_gaps(&[GRAY(0.7), GRAY(0.4)], 2, false),
            Self::Gray20 | Self::NGray20 => equal_gaps(&[GRAY(0.3), GRAY(0.8)], 20, false),
            Self::Bw20 => equal_gaps(&[GRAY(0.0), GRAY(1.0)], 20, false),
            // `SUMA_MakeColorMap_v2`: each fiducial sits at a stated index.
            Self::Byr64 => at_indices(&[
                (0, [0.0, 0.0, 0.5625]),
                (7, [0.0, 0.0, 1.0]),
                (15, [0.0, 0.5, 1.0]),
                (23, [0.0, 1.0, 1.0]),
                (31, [0.5, 1.0, 0.5625]),
                (32, [0.5625, 1.0, 0.5]),
                (40, [1.0, 1.0, 0.0]),
                (48, [1.0, 0.5, 0.0]),
                (56, [1.0, 0.0, 0.0]),
                (63, [0.5625, 0.0, 0.0]),
            ]),
            Self::Bgyr64 => at_indices(&[
                (0, [0.0, 0.0, 0.5625]),
                (7, [0.0, 0.0, 1.0]),
                (15, [0.0, 0.5, 1.0]),
                (18, [0.0, 1.0, 1.0]),
                (24, [0.0, 0.5, 0.0]),
                (32, [0.0, 1.0, 0.0]),
                (43, [1.0, 1.0, 0.0]),
                (48, [1.0, 0.5, 0.0]),
                (56, [1.0, 0.0, 0.0]),
                (63, [0.5625, 0.0, 0.0]),
            ]),
        }
    }

    /// The map as a [`ColorTable`] for the pane rules (entry 0 is the lowest value).
    pub fn to_color_table(self) -> Result<ColorTable> {
        ColorTable::new(
            self.colors()
                .into_iter()
                .map(|[r, g, b]| Rgba { r, g, b, a: 1.0 })
                .collect(),
        )
    }
}

/// `SUMA_MakeColorMap`: `n_colors` colors with the same number of steps between
/// successive fiducials. With `skip_last` the final fiducial is not itself a
/// color (it is what the map would reach next, which lets a cyclic map wrap).
///
/// The sizes are fixed in this file, so the divisibility SUMA checks at run time
/// holds by construction (and is asserted in debug builds).
fn equal_gaps(fiducials: &[[f32; 3]], n_colors: usize, skip_last: bool) -> Vec<[f32; 3]> {
    let gaps = fiducials.len() - 1;
    let intermediates = n_colors - if skip_last { gaps } else { fiducials.len() };
    debug_assert_eq!(intermediates % gaps, 0, "uneven gaps");
    let steps = intermediates / gaps + 1; // colors per gap, counting the start
    let mut colors = Vec::with_capacity(n_colors);
    for pair in fiducials.windows(2) {
        // f32 throughout, like SUMA: delta = (b - a) / steps, color = a + j * delta.
        let delta: [f32; 3] = std::array::from_fn(|c| (pair[1][c] - pair[0][c]) / steps as f32);
        for j in 0..steps {
            colors.push(std::array::from_fn(|c| pair[0][c] + j as f32 * delta[c]));
        }
    }
    if !skip_last {
        colors.push(fiducials[gaps]);
    }
    colors
}

/// `SUMA_MakeColorMap_v2`: fiducials placed at given indices, interpolated
/// linearly between them. The first index must be 0; the last is the final entry.
fn at_indices(fiducials: &[(usize, [f32; 3])]) -> Vec<[f32; 3]> {
    debug_assert_eq!(fiducials[0].0, 0);
    let mut colors = Vec::with_capacity(fiducials[fiducials.len() - 1].0 + 1);
    for pair in fiducials.windows(2) {
        let ((i0, a), (i1, b)) = (pair[0], pair[1]);
        let delta: [f32; 3] = std::array::from_fn(|c| (b[c] - a[c]) / (i1 - i0) as f32);
        for j in 0..(i1 - i0) {
            colors.push(std::array::from_fn(|c| a[c] + j as f32 * delta[c]));
        }
    }
    colors.push(fiducials[fiducials.len() - 1].1);
    colors
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizes_match_suma() {
        let sizes = [20, 19, 2, 2, 20, 20, 20, 64, 64];
        for (map, size) in SumaColorMap::ALL.into_iter().zip(sizes) {
            assert_eq!(map.colors().len(), size, "{}", map.name());
            assert_eq!(map.to_color_table().unwrap().len(), size);
        }
    }

    #[test]
    fn endpoints_and_cycle() {
        // bgyr19 runs blue (lowest) to red (highest) and includes both ends.
        let bgyr = SumaColorMap::Bgyr19.colors();
        assert_eq!((bgyr[0], bgyr[18]), ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]));
        // rgybr20 leaves out the final red, so it starts at red but never ends on it.
        let rgybr = SumaColorMap::Rgybr20.colors();
        assert_eq!(rgybr[0], [1.0, 0.0, 0.0]);
        assert_ne!(rgybr[19], [1.0, 0.0, 0.0]);
    }

    #[test]
    fn names_round_trip_and_alias() {
        for map in SumaColorMap::ALL {
            assert_eq!(SumaColorMap::from_name(map.name()), Some(map));
        }
        assert_eq!(
            SumaColorMap::from_name("Matlab_Default_BYR64"),
            Some(SumaColorMap::Byr64)
        );
        assert_eq!(SumaColorMap::from_name("viridis"), None);
    }
}
