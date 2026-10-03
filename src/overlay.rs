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
// Overlay evaluation: turn data (an intensity column, optionally a separate
// threshold column and a brightness column) plus an immutable display intent
// (`OverlaySpec`) into colors, a pass/fail mask, and diagnostics
// (`OverlayEvaluation`). This is the pure, CPU reference for what a viewer shows;
// a GPU path must reproduce it (see `GpuOverlayParams`).
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * Consumes `Dataset`/`DataColumn` (dataset.rs, column.rs), color maps
//   (color.rs, afni_colors.rs, labels.rs) and thresholds (threshold.rs).
// * Produces plain colors that `composite.rs` layers over an underlay.
// * Ported from sumaru's `overlay.rs` with its viewer state removed: there is no
//   mutable cache here. Callers may cache an `OverlayEvaluation`; invalidation is
//   theirs, because the spec is immutable data.
//
// ORDER OF OPERATIONS (each step is a decision; this is the contract)
//
//   1. A missing (NaN/infinite) intensity gets `missing_color` and is not "passed".
//   2. The value is mapped to a color (continuous, pane lookup, or label key).
//   3. The color map brightness factor (SUMA `-br`) is applied to the map.
//   4. Brightness modulation by a second column multiplies the color channels.
//   5. Intensity-range clipping ("hide outside range") and "hide zero" clear alpha.
//   6. The threshold: passing samples are untouched; failing samples are hidden,
//      dimmed, or faded according to `FailedThreshold`.
//   7. Cluster survivors: a sample outside every surviving cluster is hidden,
//      unconditionally (it passed the sample-wise threshold, so the fade ramp
//      would otherwise leave it opaque).
//   8. The overlay opacity multiplies alpha.
//
// SINGLE PRECISION WHERE AFNI/SUMA ARE: pane lookups compute the color index in
// `f32`, because the index is a truncation and a value exactly on a pane boundary
// must land in the same pane as in SUMA/AFNI (computed in `float`).
// ---------------------------------------------------------------------------

//! Overlay evaluation: thresholds, color maps, brightness, and opacity applied to data.

use crate::color::{ContinuousColorMap, Rgba};
use crate::column::{ColumnRange, DataColumn};
use crate::dataset::Dataset;
use crate::error::{Error, Result};
use crate::labels::LabelColorMap;
use crate::numeric::ensure_finite;
use crate::threshold::{FadeModel, MissingThreshold, Threshold};

// ---------------------------------------------------------------------------
// Color tables and pane rules
// ---------------------------------------------------------------------------

/// An ordered list of colors for pane-style lookup. Entry 0 is the color for the
/// LOWEST value and the last entry for the highest.
#[derive(Debug, Clone, PartialEq)]
pub struct ColorTable {
    colors: Vec<Rgba>,
}

impl ColorTable {
    /// A table; needs at least two colors.
    pub fn new(colors: Vec<Rgba>) -> Result<Self> {
        if colors.len() < 2 {
            return Err(Error::InvalidParameter {
                name: "color table".into(),
                reason: format!("needs at least 2 colors, got {}", colors.len()),
            });
        }
        Ok(Self { colors })
    }

    /// A table from opaque 8-bit colors, lowest value first.
    pub fn from_rgb_bytes(rgb: &[[u8; 3]]) -> Result<Self> {
        Self::new(
            rgb.iter()
                .map(|&[r, g, b]| Rgba::from_u8(r, g, b, 255))
                .collect(),
        )
    }

    /// Number of colors.
    pub fn len(&self) -> usize {
        self.colors.len()
    }

    /// Always false (a table has at least two colors).
    pub fn is_empty(&self) -> bool {
        false
    }

    /// The colors, lowest value first.
    pub fn colors(&self) -> &[Rgba] {
        &self.colors
    }

    /// A copy with every color's RGB multiplied by `factor` and clamped to 1:
    /// SUMA's "brightness factor" applied to a color map (`-br`).
    fn brightened(&self, factor: f32) -> ColorTable {
        let scale = |c: Rgba| Rgba {
            r: (c.r * factor).min(1.0),
            g: (c.g * factor).min(1.0),
            b: (c.b * factor).min(1.0),
            a: c.a,
        };
        ColorTable {
            colors: self.colors.iter().copied().map(scale).collect(),
        }
    }
}

/// How a value is turned into an index into a [`ColorTable`] of `N` panes.
///
/// "Panes" is AFNI/SUMA's model of a color bar: `N` equal-height colored cells.
/// Entry `i` sits at the BOTTOM of pane `i`, so with interpolation the last color
/// is reached at `(N - 1) / N` of the range and held to the top. That differs from
/// a continuous map with `N` stops at `i / (N - 1)`, which is why panes are a
/// separate mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PaneRule {
    /// SUMA `SUMA_INTERP`: `s = (v - min)/(max - min) * N` clamped to `[0, N]`,
    /// pane `i = min(floor(s), N-1)`, color mixed between entry `i` and `i + 1`
    /// by `s - i` (the last pane has no neighbor and keeps its color). A range of
    /// zero width gives the middle entry `(N - 1) / 2`.
    SumaInterpolated,
    /// SUMA `SUMA_NO_INTERP`: the same pane index, no mixing.
    SumaBanded,
    /// AFNI's volume overlay: pane `j = floor(N/(max - min) * (max - v))` clamped to
    /// `[0, N-1]` counted from the TOP of the bar; no mixing. Differs from
    /// `SumaBanded` exactly on pane boundaries.
    AfniBanded,
    /// SUMA `SUMA_DIRECT`: the value, truncated to an integer, is the table index
    /// (negative values give 0; values past the end give the last entry). The
    /// range is ignored. For ROI-index style data.
    Direct,
}

impl PaneRule {
    /// The color for `value` in `table` over `range` (all arithmetic in `f32`,
    /// as in SUMA/AFNI).
    fn lookup(self, table: &ColorTable, value: f64, range: ColumnRange) -> Result<Rgba> {
        let n = table.len();
        let colors = table.colors();
        let (v, lo, hi) = (value as f32, range.min as f32, range.max as f32);
        Ok(match self {
            PaneRule::Direct => {
                let index = v as i64; // truncation toward zero, like a C cast
                colors[index.clamp(0, n as i64 - 1) as usize]
            }
            PaneRule::AfniBanded => {
                if lo >= hi {
                    return Err(Error::InvalidParameter {
                        name: "intensity range".into(),
                        reason: "AFNI pane lookup needs min < max".into(),
                    });
                }
                let factor = n as f32 / (hi - lo);
                let from_top = (factor * (hi - v)) as i64; // `(int)` truncation
                let j = from_top.clamp(0, n as i64 - 1) as usize;
                colors[n - 1 - j] // AFNI counts panes from the top
            }
            PaneRule::SumaInterpolated | PaneRule::SumaBanded => {
                let span = hi - lo;
                if span < 0.0 {
                    return Err(Error::InvalidParameter {
                        name: "intensity range".into(),
                        reason: format!("max {hi} is below min {lo}"),
                    });
                }
                if span == 0.0 {
                    return Ok(colors[(n - 1) / 2]);
                }
                let scaled = ((v - lo) / span * n as f32).clamp(0.0, n as f32);
                let i0 = (scaled as usize).min(n - 1);
                if self == PaneRule::SumaBanded || i0 + 1 >= n {
                    colors[i0]
                } else {
                    colors[i0].lerp(colors[i0 + 1], scaled - i0 as f32)
                }
            }
        })
    }
}

/// How intensity values become colors.
#[derive(Debug, Clone, PartialEq)]
pub enum OverlayColors {
    /// A continuous map sampled at the value's position in the display range.
    Continuous(ContinuousColorMap),
    /// A table of panes with one of the AFNI/SUMA lookup rules.
    Panes {
        /// The pane colors, lowest value first.
        table: ColorTable,
        /// How a value picks and mixes panes.
        rule: PaneRule,
    },
    /// Integer label keys looked up in a label table; the value must be an
    /// integer, otherwise it counts as missing.
    Labels(LabelColorMap),
}

impl OverlayColors {
    /// Whether this mapping needs an intensity range (labels and direct panes do not).
    fn needs_range(&self) -> bool {
        !matches!(
            self,
            OverlayColors::Labels(_)
                | OverlayColors::Panes {
                    rule: PaneRule::Direct,
                    ..
                }
        )
    }

    /// The color for `value`, or `None` if the mapping cannot color it.
    fn color_for(&self, value: f64, range: ColumnRange, brightness: f32) -> Result<Option<Rgba>> {
        Ok(match self {
            OverlayColors::Continuous(map) => Some(
                map.sample(range.normalized(value))
                    .scaled_clamped(brightness),
            ),
            OverlayColors::Panes { table, rule } => {
                // The map-level brightness factor is applied to the table first, as
                // SUMA does, so mixing sees the brightened entries.
                if brightness == 1.0 {
                    Some(rule.lookup(table, value, range)?)
                } else {
                    Some(rule.lookup(&table.brightened(brightness), value, range)?)
                }
            }
            OverlayColors::Labels(map) => {
                // Only exact integers are keys; anything else is not a label.
                if value.fract() == 0.0 && value.abs() < 9.0e15 {
                    Some(map.color_for_key(value as i64).scaled_clamped(brightness))
                } else {
                    None
                }
            }
        })
    }
}

impl Rgba {
    /// RGB scaled and clamped to 1 (so a factor above 1 brightens); alpha untouched.
    fn scaled_clamped(self, factor: f32) -> Rgba {
        if factor == 1.0 {
            self
        } else {
            Rgba {
                r: (self.r * factor).min(1.0),
                g: (self.g * factor).min(1.0),
                b: (self.b * factor).min(1.0),
                a: self.a,
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The display intent
// ---------------------------------------------------------------------------

/// Where a display range comes from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RangeSelection {
    /// The min and max of the finite values, computed now.
    Computed,
    /// The range recorded in the file (`COLMS_RANGE`); an error if there is none.
    /// Kept distinct from `Computed` because a stale header and the data can disagree.
    Recorded,
    /// A range the user chose.
    Manual(ColumnRange),
}

/// What happens to values outside the intensity range.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ClipMode {
    /// Take the end color (SUMA's and AFNI's behavior).
    #[default]
    Clamp,
    /// Hide the sample (sumaru's option; AFNI's `ZBELOW`/`ZABOVE` flags).
    Hide,
}

/// Modulate color brightness by a second column (SUMA's brightness column).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrightnessModulation {
    /// Range of the brightness column. `Manual` also CLIPS the column to the range,
    /// as SUMA does when both clip values are given.
    pub range: RangeSelection,
    /// Multipliers for the color channels at the low and high end of the range:
    /// `factor = scale[0] + position * (scale[1] - scale[0])`. A missing
    /// brightness value leaves the color unmodulated.
    pub scale: [f32; 2],
}

/// Extra, non-AFNI styling for faded samples (sumaru enhancements). All zero by
/// default, which is AFNI's behavior of changing opacity only.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct FadeStyle {
    /// Pull failing colors toward their own luminance by `shortfall * desaturate`.
    pub desaturate: f32,
    /// Darken failing colors by `shortfall * darken`.
    pub darken: f32,
    /// Push PASSING colors away from their luminance by `boost` (limited by gamut).
    pub boost: f32,
    /// A ceiling on the opacity of failing samples, keeping a visible step at the
    /// threshold even where the ramp is nearly flat. `None` for no ceiling.
    pub max_alpha: Option<f32>,
}

/// What to do with samples that fail the threshold.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FailedThreshold {
    /// Draw them anyway (the threshold only affects `passed`).
    Keep,
    /// Make them transparent.
    Hide,
    /// Multiply their color channels by this factor in `[0, 1]`.
    Dim(f32),
    /// Fade their opacity (and optionally style them) per `model`.
    Fade {
        /// The opacity ramp (AFNI, SUMA, or generalized).
        model: FadeModel,
        /// Optional non-AFNI styling.
        style: FadeStyle,
    },
}

/// Immutable display intent for one overlay.
#[derive(Debug, Clone, PartialEq)]
pub struct OverlaySpec {
    /// How values become colors.
    pub colors: OverlayColors,
    /// Where the display range comes from.
    pub intensity_range: RangeSelection,
    /// Replace the range by `[-e, e]` with `e = max(|min|, |max|)` (SUMA `-anr`).
    pub symmetric_range: bool,
    /// Treatment of values outside the range.
    pub clip: ClipMode,
    /// The display threshold (on the threshold column, or the intensity column if
    /// none is given).
    pub threshold: Threshold,
    /// What failing samples look like.
    pub failed: FailedThreshold,
    /// What a missing threshold value means.
    pub missing_threshold: MissingThreshold,
    /// Whether an intensity of exactly 0 is drawn (SUMA `shw_0`; false masks zeros).
    pub show_zero: bool,
    /// Overall opacity in `[0, 1]`.
    pub opacity: f32,
    /// Optional brightness modulation by another column.
    pub brightness: Option<BrightnessModulation>,
    /// Multiply the COLOR MAP by this factor in `(0, 2]` (SUMA `-br`); 1 = none.
    pub colormap_brightness: f32,
    /// The color of a sample whose intensity is missing or cannot be colored.
    pub missing_color: Rgba,
}

impl OverlaySpec {
    /// A spec with neutral defaults: computed range, clamping, no threshold,
    /// zeros shown, full opacity, no brightness changes, transparent missing color.
    pub fn new(colors: OverlayColors) -> Self {
        Self {
            colors,
            intensity_range: RangeSelection::Computed,
            symmetric_range: false,
            clip: ClipMode::Clamp,
            threshold: Threshold::Off,
            failed: FailedThreshold::Hide,
            missing_threshold: MissingThreshold::Hide,
            show_zero: true,
            opacity: 1.0,
            brightness: None,
            colormap_brightness: 1.0,
            missing_color: Rgba::TRANSPARENT,
        }
    }

    /// Check every setting before any data is touched.
    pub fn validate(&self) -> Result<()> {
        self.threshold.validate()?;
        if let FailedThreshold::Fade { model, style } = &self.failed {
            model.validate(&self.threshold)?;
            for (name, v) in [
                ("desaturate", style.desaturate),
                ("darken", style.darken),
                ("boost", style.boost),
            ] {
                ensure_finite(name, f64::from(v))?;
                if !(0.0..=1.0).contains(&v) {
                    return Err(Error::InvalidParameter {
                        name: name.into(),
                        reason: format!("{v} is outside [0, 1]"),
                    });
                }
            }
            if let Some(m) = style.max_alpha {
                ensure_finite("max_alpha", f64::from(m))?;
            }
        }
        if let FailedThreshold::Dim(f) = self.failed {
            ensure_finite("dim factor", f64::from(f))?;
            if !(0.0..=1.0).contains(&f) {
                return Err(Error::InvalidParameter {
                    name: "dim factor".into(),
                    reason: format!("{f} is outside [0, 1]"),
                });
            }
        }
        ensure_finite("opacity", f64::from(self.opacity))?;
        if !(0.0..=1.0).contains(&self.opacity) {
            return Err(Error::InvalidParameter {
                name: "opacity".into(),
                reason: format!("{} is outside [0, 1]", self.opacity),
            });
        }
        ensure_finite("colormap brightness", f64::from(self.colormap_brightness))?;
        if !(self.colormap_brightness > 0.0 && self.colormap_brightness <= 2.0) {
            return Err(Error::InvalidParameter {
                name: "colormap brightness".into(),
                reason: format!("{} is outside (0, 2]", self.colormap_brightness),
            });
        }
        if let Some(b) = &self.brightness {
            for s in b.scale {
                ensure_finite("brightness scale", f64::from(s))?;
            }
            if let RangeSelection::Manual(r) = b.range {
                ColumnRange::new(r.min, r.max)?;
            }
        }
        if let RangeSelection::Manual(r) = self.intensity_range {
            ColumnRange::new(r.min, r.max)?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Inputs and outputs
// ---------------------------------------------------------------------------

/// Data for one evaluation, as parallel slices (one entry per ROW).
#[derive(Debug, Clone, Copy, Default)]
pub struct OverlayInputs<'a> {
    /// The values to color.
    pub intensity: &'a [f64],
    /// The values to threshold. `None` thresholds the intensity itself.
    pub threshold: Option<&'a [f64]>,
    /// The values that modulate brightness, if any.
    pub brightness: Option<&'a [f64]>,
    /// `true` for rows that belong to a surviving cluster. `None` means no cluster
    /// restriction.
    pub cluster_survivors: Option<&'a [bool]>,
    /// The intensity range recorded in the file, for `RangeSelection::Recorded`.
    pub recorded_intensity_range: Option<ColumnRange>,
    /// The brightness range recorded in the file.
    pub recorded_brightness_range: Option<ColumnRange>,
}

/// What an evaluation found, for display and debugging.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct OverlayDiagnostics {
    /// Rows evaluated.
    pub rows: usize,
    /// Rows whose intensity was missing or could not be colored.
    pub missing_intensity: usize,
    /// Rows hidden because the intensity was exactly zero and zeros are not shown.
    pub hidden_zero: usize,
    /// Rows hidden because the intensity was outside the range (`ClipMode::Hide`).
    pub hidden_out_of_range: usize,
    /// Rows that failed the threshold.
    pub failed_threshold: usize,
    /// Rows that PASSED the threshold but were hidden because they are not in a
    /// surviving cluster.
    pub rejected_by_cluster: usize,
    /// The intensity range actually used, if one was needed.
    pub intensity_range: Option<ColumnRange>,
    /// The brightness range actually used, if brightness modulation was on.
    pub brightness_range: Option<ColumnRange>,
}

/// Colors and pass/fail flags per row, plus diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct RowEvaluation {
    /// One color per row (transparent where nothing is drawn).
    pub colors: Vec<Rgba>,
    /// `true` for rows that are drawn at full strength: finite intensity, passing
    /// the threshold, in a surviving cluster, not hidden by zero or range rules.
    pub passed: Vec<bool>,
    /// Counts and the ranges used.
    pub diagnostics: OverlayDiagnostics,
}

fn finite_range(values: &[f64]) -> Option<ColumnRange> {
    let mut it = values.iter().copied().filter(|v| v.is_finite());
    let first = it.next()?;
    let (min, max) = it.fold((first, first), |(lo, hi), v| (lo.min(v), hi.max(v)));
    Some(ColumnRange { min, max })
}

fn resolve_range(
    selection: RangeSelection,
    values: &[f64],
    recorded: Option<ColumnRange>,
    what: &str,
) -> Result<ColumnRange> {
    match selection {
        RangeSelection::Manual(r) => Ok(r),
        RangeSelection::Recorded => recorded.ok_or_else(|| Error::InvalidParameter {
            name: format!("{what} range"),
            reason: "no range was recorded in the file".into(),
        }),
        RangeSelection::Computed => finite_range(values).ok_or_else(|| Error::InvalidParameter {
            name: format!("{what} range"),
            reason: "there are no finite values to take a range from".into(),
        }),
    }
}

/// Evaluate an overlay over aligned rows.
///
/// All input slices must have the same length as `inputs.intensity`. See the
/// module docs for the order of operations.
pub fn evaluate_rows(spec: &OverlaySpec, inputs: &OverlayInputs<'_>) -> Result<RowEvaluation> {
    spec.validate()?;
    let n = inputs.intensity.len();
    let check = |what: &str, len: Option<usize>| match len {
        Some(l) if l != n => Err(Error::LengthMismatch {
            what: what.into(),
            expected: n,
            found: l,
        }),
        _ => Ok(()),
    };
    check("threshold values", inputs.threshold.map(<[f64]>::len))?;
    check("brightness values", inputs.brightness.map(<[f64]>::len))?;
    check(
        "cluster survivors",
        inputs.cluster_survivors.map(<[bool]>::len),
    )?;

    // Ranges.
    let intensity_range = if spec.colors.needs_range() {
        let mut r = resolve_range(
            spec.intensity_range,
            inputs.intensity,
            inputs.recorded_intensity_range,
            "intensity",
        )?;
        if spec.symmetric_range {
            let e = r.min.abs().max(r.max.abs());
            r = ColumnRange { min: -e, max: e };
        }
        Some(r)
    } else {
        None
    };
    let brightness_range =
        match (&spec.brightness, inputs.brightness) {
            (Some(b), Some(values)) => Some(resolve_range(
                b.range,
                values,
                inputs.recorded_brightness_range,
                "brightness",
            )?),
            (Some(_), None) => return Err(Error::InvalidParameter {
                name: "brightness".into(),
                reason:
                    "the spec asks for brightness modulation but no brightness values were given"
                        .into(),
            }),
            _ => None,
        };

    let mut diag = OverlayDiagnostics {
        rows: n,
        intensity_range,
        brightness_range,
        ..Default::default()
    };
    let mut colors = vec![Rgba::TRANSPARENT; n];
    let mut passed = vec![false; n];
    let threshold_values = inputs.threshold.unwrap_or(inputs.intensity);
    let dummy_range = ColumnRange { min: 0.0, max: 1.0 };

    for row in 0..n {
        let value = inputs.intensity[row];
        // 1. Missing intensity.
        if !value.is_finite() {
            colors[row] = spec
                .missing_color
                .with_alpha(spec.missing_color.a * spec.opacity);
            diag.missing_intensity += 1;
            continue;
        }
        // 2-3. Color, with the map brightness factor.
        let Some(mut color) = spec.colors.color_for(
            value,
            intensity_range.unwrap_or(dummy_range),
            spec.colormap_brightness,
        )?
        else {
            colors[row] = spec
                .missing_color
                .with_alpha(spec.missing_color.a * spec.opacity);
            diag.missing_intensity += 1;
            continue;
        };
        // 4. Brightness modulation.
        if let (Some(b), Some(values), Some(range)) =
            (&spec.brightness, inputs.brightness, brightness_range)
        {
            let mut bv = values[row];
            if bv.is_finite() {
                if matches!(b.range, RangeSelection::Manual(_)) {
                    bv = bv.clamp(range.min, range.max); // SUMA clips B to a manual range
                }
                let position = range.normalized(bv) as f32;
                let factor = b.scale[0] + position * (b.scale[1] - b.scale[0]);
                color = Rgba {
                    r: color.r * factor,
                    g: color.g * factor,
                    b: color.b * factor,
                    a: color.a,
                };
                // Channels can leave [0, 1] with a factor above 1 or below 0.
                color = Rgba {
                    r: color.r.clamp(0.0, 1.0),
                    g: color.g.clamp(0.0, 1.0),
                    b: color.b.clamp(0.0, 1.0),
                    a: color.a,
                };
            }
        }
        let mut visible = true;
        // 5. Zero and out-of-range hiding.
        if !spec.show_zero && value == 0.0 {
            visible = false;
            diag.hidden_zero += 1;
        }
        if spec.clip == ClipMode::Hide {
            if let Some(r) = intensity_range {
                if !r.contains(value) {
                    visible = false;
                    diag.hidden_out_of_range += 1;
                }
            }
        }
        if !visible {
            colors[row] = color.with_alpha(0.0);
            continue;
        }
        // 6. Threshold.
        let tv = threshold_values[row];
        let passes_threshold = if tv.is_finite() || matches!(spec.threshold, Threshold::Off) {
            spec.threshold.passes(tv)
        } else {
            spec.missing_threshold == MissingThreshold::Show
        };
        let mut row_passed = passes_threshold;
        if passes_threshold {
            if let FailedThreshold::Fade { style, .. } = spec.failed {
                if style.boost > 0.0 {
                    color = boost_color(color, style.boost);
                }
            }
        } else {
            diag.failed_threshold += 1;
            match spec.failed {
                FailedThreshold::Keep => {}
                FailedThreshold::Hide => color = color.with_alpha(0.0),
                FailedThreshold::Dim(f) => {
                    color = Rgba {
                        r: color.r * f,
                        g: color.g * f,
                        b: color.b * f,
                        a: color.a,
                    }
                }
                FailedThreshold::Fade { model, style } => {
                    let mut factor = if tv.is_finite() {
                        model.factor(&spec.threshold, tv)
                    } else {
                        0.0
                    };
                    if let Some(cap) = style.max_alpha {
                        factor = factor.min(cap.clamp(0.0, 1.0));
                    }
                    color = style_faded(color, factor, style);
                }
            }
        }
        // 7. Clusters.
        if let Some(survivors) = inputs.cluster_survivors {
            if !survivors[row] {
                color = color.with_alpha(0.0);
                // Count only rows the cluster rule itself removed: ones that passed the
                // threshold. A row that already failed it was hidden for that reason.
                if passes_threshold {
                    diag.rejected_by_cluster += 1;
                }
                row_passed = false;
            }
        }
        // 8. Opacity.
        color = color.with_alpha(color.a * spec.opacity);
        colors[row] = color;
        passed[row] = row_passed;
    }
    Ok(RowEvaluation {
        colors,
        passed,
        diagnostics: diag,
    })
}

/// Rec. 709 luma of a color's RGB.
fn luma(c: Rgba) -> f32 {
    0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b
}

/// Apply a fade factor and the optional styling to a failing color.
fn style_faded(mut color: Rgba, factor: f32, style: FadeStyle) -> Rgba {
    let factor = factor.clamp(0.0, 1.0);
    color.a *= factor;
    let shortfall = 1.0 - factor;
    if style.desaturate > 0.0 {
        let y = luma(color);
        let amount = shortfall * style.desaturate;
        color.r += (y - color.r) * amount;
        color.g += (y - color.g) * amount;
        color.b += (y - color.b) * amount;
    }
    if style.darken > 0.0 {
        let scale = 1.0 - shortfall * style.darken;
        color.r *= scale;
        color.g *= scale;
        color.b *= scale;
    }
    color
}

/// Push a color away from its own luma, clamped to the gamut.
fn boost_color(c: Rgba, boost: f32) -> Rgba {
    let y = luma(c);
    let push = |v: f32| (y + (v - y) * (1.0 + boost)).clamp(0.0, 1.0);
    Rgba {
        r: push(c.r),
        g: push(c.g),
        b: push(c.b),
        a: c.a,
    }
}

// ---------------------------------------------------------------------------
// Datasets
// ---------------------------------------------------------------------------

/// Which columns of a dataset feed an overlay.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OverlayColumns {
    /// Column holding the values to color.
    pub intensity: usize,
    /// Column to threshold on, if different from the intensity.
    pub threshold: Option<usize>,
    /// Column that modulates brightness, if any.
    pub brightness: Option<usize>,
}

impl OverlayColumns {
    /// Color by column `intensity`, thresholding on the same column.
    pub fn new(intensity: usize) -> Self {
        Self {
            intensity,
            threshold: None,
            brightness: None,
        }
    }
}

/// A finished overlay over a whole domain.
#[derive(Debug, Clone, PartialEq)]
pub struct OverlayEvaluation {
    /// One color per DOMAIN sample (node or voxel); transparent where the dataset
    /// has no row or nothing is drawn. Dense, ready for a viewer or a GPU buffer.
    pub colors: Vec<Rgba>,
    /// `true` per domain sample that is drawn at full strength (see
    /// [`RowEvaluation::passed`]).
    pub passed: Vec<bool>,
    /// Counts and ranges (over the dataset's rows).
    pub diagnostics: OverlayDiagnostics,
}

fn column_values(column: &DataColumn) -> Result<Vec<f64>> {
    if !column.values().is_numeric() {
        return Err(Error::InvalidParameter {
            name: format!("column '{}'", column.label()),
            reason: "is not numeric".into(),
        });
    }
    Ok((0..column.len())
        .map(|r| column.values().get_f64(r).unwrap_or(f64::NAN))
        .collect())
}

/// Evaluate an overlay on a [`Dataset`], expanding sparse rows to the whole domain.
///
/// `cluster_survivors`, if given, has one entry per DOMAIN sample (`true` = the
/// sample is in a surviving cluster).
pub fn evaluate_dataset(
    spec: &OverlaySpec,
    dataset: &Dataset,
    columns: &OverlayColumns,
    cluster_survivors: Option<&[bool]>,
) -> Result<OverlayEvaluation> {
    let pick = |index: usize| -> Result<&DataColumn> {
        dataset
            .columns()
            .get(index)
            .ok_or_else(|| Error::IndexOutOfRange {
                index: index as i64,
                len: dataset.columns().len(),
            })
    };
    let intensity_col = pick(columns.intensity)?;
    let intensity = column_values(intensity_col)?;
    let threshold = columns
        .threshold
        .map(|i| pick(i).and_then(column_values))
        .transpose()?;
    let brightness_col = columns.brightness.map(pick).transpose()?;
    let brightness = brightness_col.map(column_values).transpose()?;

    let domain_len = dataset.domain().sample_count();
    if let Some(cs) = cluster_survivors {
        if cs.len() != domain_len {
            return Err(Error::LengthMismatch {
                what: "cluster survivors (per domain sample)".into(),
                expected: domain_len,
                found: cs.len(),
            });
        }
    }
    // Per-row cluster flags from the per-sample list.
    let row_survivors: Option<Vec<bool>> = cluster_survivors.map(|cs| {
        (0..dataset.row_count())
            .map(|row| dataset.sample_for_row(row).is_some_and(|s| cs[s as usize]))
            .collect()
    });

    let inputs = OverlayInputs {
        intensity: &intensity,
        threshold: threshold.as_deref(),
        brightness: brightness.as_deref(),
        cluster_survivors: row_survivors.as_deref(),
        recorded_intensity_range: intensity_col.range_report().recorded.map(|r| r.range),
        recorded_brightness_range: brightness_col
            .and_then(|c| c.range_report().recorded.map(|r| r.range)),
    };
    let rows = evaluate_rows(spec, &inputs)?;

    let mut colors = vec![Rgba::TRANSPARENT; domain_len];
    let mut passed = vec![false; domain_len];
    for row in 0..dataset.row_count() {
        if let Some(sample) = dataset.sample_for_row(row) {
            colors[sample as usize] = rows.colors[row];
            passed[sample as usize] = rows.passed[row];
        }
    }
    Ok(OverlayEvaluation {
        colors,
        passed,
        diagnostics: rows.diagnostics,
    })
}

// ---------------------------------------------------------------------------
// GPU-ready parameters
// ---------------------------------------------------------------------------

/// The threshold in the flat form a shader uniform wants.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GpuThreshold {
    /// 0 off, 1 above, 2 below, 3 between, 4 outside, 5 absolute-above.
    pub mode: u32,
    /// First bound (`t`, or `lo`).
    pub a: f32,
    /// Second bound (`hi`), or 0.
    pub b: f32,
}

/// Everything a viewer needs to reproduce an overlay on the GPU: a lookup table
/// and a handful of uniforms. Shader code and buffers stay in the viewer; this
/// only fixes the numbers, so the CPU evaluation remains the reference.
#[derive(Debug, Clone, PartialEq)]
pub struct GpuOverlayParams {
    /// Colors sampled at `i / (len - 1)` across the display range. Sample it with
    /// nearest or linear filtering; for pane maps use nearest and a table whose
    /// length is a multiple of the pane count for banding to survive.
    pub lookup: Vec<[f32; 4]>,
    /// Display range `[min, max]` mapping data to texture coordinates.
    pub range: [f32; 2],
    /// The threshold.
    pub threshold: GpuThreshold,
    /// Overall opacity.
    pub opacity: f32,
    /// Whether zero is hidden.
    pub hide_zero: bool,
    /// Whether values outside the range are hidden (rather than clamped).
    pub hide_out_of_range: bool,
}

impl OverlaySpec {
    /// Flatten this spec for a shader given the resolved display `range`.
    /// Label maps are not table-lookups over a range and are rejected.
    pub fn gpu_params(
        &self,
        range: ColumnRange,
        lookup_entries: usize,
    ) -> Result<GpuOverlayParams> {
        self.validate()?;
        if lookup_entries < 2 {
            return Err(Error::InvalidParameter {
                name: "lookup entries".into(),
                reason: "need at least 2".into(),
            });
        }
        let last = (lookup_entries - 1) as f64;
        let lookup = (0..lookup_entries)
            .map(|i| {
                let position = i as f64 / last;
                let value = range.min + position * (range.max - range.min);
                self.colors
                    .color_for(value, range, self.colormap_brightness)
                    .map(|c| c.unwrap_or(self.missing_color).to_array())
            })
            .collect::<Result<Vec<_>>>()
            .and_then(|v| {
                if matches!(self.colors, OverlayColors::Labels(_)) {
                    Err(Error::Unsupported(
                        "label maps cannot be baked into a range lookup table".into(),
                    ))
                } else {
                    Ok(v)
                }
            })?;
        let threshold = match self.threshold {
            Threshold::Off => GpuThreshold {
                mode: 0,
                a: 0.0,
                b: 0.0,
            },
            Threshold::Above(t) => GpuThreshold {
                mode: 1,
                a: t as f32,
                b: 0.0,
            },
            Threshold::Below(t) => GpuThreshold {
                mode: 2,
                a: t as f32,
                b: 0.0,
            },
            Threshold::Between { lo, hi } => GpuThreshold {
                mode: 3,
                a: lo as f32,
                b: hi as f32,
            },
            Threshold::Outside { lo, hi } => GpuThreshold {
                mode: 4,
                a: lo as f32,
                b: hi as f32,
            },
            Threshold::AbsoluteAbove(t) => GpuThreshold {
                mode: 5,
                a: t as f32,
                b: 0.0,
            },
        };
        Ok(GpuOverlayParams {
            lookup,
            range: [range.min as f32, range.max as f32],
            threshold,
            opacity: self.opacity,
            hide_zero: !self.show_zero,
            hide_out_of_range: self.clip == ClipMode::Hide,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::afni_colors::AfniColorScale;
    use crate::column::{ColumnData, ColumnRole};
    use crate::domain::{Domain, SurfaceDomain};
    use crate::labels::{LabelEntry, LabelTable};
    use crate::threshold::{FadeCurve, FadeWidth};

    fn gray_map() -> OverlayColors {
        OverlayColors::Continuous(ContinuousColorMap::grayscale())
    }

    fn spec() -> OverlaySpec {
        let mut s = OverlaySpec::new(gray_map());
        s.intensity_range = RangeSelection::Manual(ColumnRange::new(0.0, 10.0).unwrap());
        s
    }

    fn run(spec: &OverlaySpec, intensity: &[f64]) -> RowEvaluation {
        evaluate_rows(
            spec,
            &OverlayInputs {
                intensity,
                ..Default::default()
            },
        )
        .unwrap()
    }

    #[test]
    fn continuous_colors_follow_the_range_and_clamp() {
        let r = run(&spec(), &[0.0, 5.0, 10.0, 20.0, -5.0]);
        assert_eq!(r.colors[0].r, 0.0);
        assert_eq!(r.colors[1].r, 0.5);
        assert_eq!(r.colors[2].r, 1.0);
        assert_eq!(r.colors[3].r, 1.0, "above the range takes the end color");
        assert_eq!(r.colors[4].r, 0.0);
        assert!(r.passed.iter().all(|&p| p));
        assert_eq!(
            r.diagnostics.intensity_range,
            Some(ColumnRange {
                min: 0.0,
                max: 10.0
            })
        );
    }

    #[test]
    fn computed_recorded_and_symmetric_ranges() {
        let mut s = OverlaySpec::new(gray_map());
        let r = run(&s, &[2.0, 4.0, f64::NAN, 6.0]);
        assert_eq!(
            r.diagnostics.intensity_range,
            Some(ColumnRange { min: 2.0, max: 6.0 }),
            "NaN ignored"
        );
        s.symmetric_range = true;
        let r = run(&s, &[-2.0, 6.0]);
        assert_eq!(
            r.diagnostics.intensity_range,
            Some(ColumnRange {
                min: -6.0,
                max: 6.0
            })
        );
        // Recorded needs a recorded range to exist, and is distinct from computed.
        s.symmetric_range = false;
        s.intensity_range = RangeSelection::Recorded;
        assert!(evaluate_rows(
            &s,
            &OverlayInputs {
                intensity: &[1.0, 2.0],
                ..Default::default()
            }
        )
        .is_err());
        let rec = ColumnRange::new(0.0, 100.0).unwrap();
        let r = evaluate_rows(
            &s,
            &OverlayInputs {
                intensity: &[1.0, 2.0],
                recorded_intensity_range: Some(rec),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(r.diagnostics.intensity_range, Some(rec));
        // Computed with nothing finite is an error, not a silent default.
        s.intensity_range = RangeSelection::Computed;
        assert!(evaluate_rows(
            &s,
            &OverlayInputs {
                intensity: &[f64::NAN],
                ..Default::default()
            }
        )
        .is_err());
    }

    #[test]
    fn missing_intensity_gets_the_missing_color_and_is_not_passed() {
        let mut s = spec();
        s.missing_color = Rgba::from_u8(90, 90, 90, 255);
        let r = run(&s, &[5.0, f64::NAN, f64::INFINITY]);
        assert_eq!(r.colors[1], Rgba::from_u8(90, 90, 90, 255));
        assert_eq!(r.colors[2], r.colors[1]);
        assert_eq!(r.passed, vec![true, false, false]);
        assert_eq!(r.diagnostics.missing_intensity, 2);
        s.missing_color = Rgba::TRANSPARENT;
        assert_eq!(run(&s, &[f64::NAN]).colors[0].a, 0.0);
    }

    #[test]
    fn show_zero_and_out_of_range_hiding() {
        let mut s = spec();
        s.show_zero = false;
        s.clip = ClipMode::Hide;
        let r = run(&s, &[0.0, 5.0, 11.0, -1.0]);
        assert_eq!(
            r.colors.iter().map(|c| c.a).collect::<Vec<_>>(),
            vec![0.0, 1.0, 0.0, 0.0]
        );
        assert_eq!(r.passed, vec![false, true, false, false]);
        assert_eq!(
            (r.diagnostics.hidden_zero, r.diagnostics.hidden_out_of_range),
            (1, 2)
        );
        // A hidden sample keeps its color channels (alpha only), for debugging.
        assert_eq!(r.colors[2].r, 1.0);
    }

    #[test]
    fn threshold_hide_dim_and_keep() {
        let mut s = spec();
        s.threshold = Threshold::AbsoluteAbove(4.0);
        let v = [1.0, 4.0, -5.0, f64::NAN];
        s.failed = FailedThreshold::Hide;
        let r = run(&s, &v);
        assert_eq!(r.passed, vec![false, true, true, false]);
        assert_eq!(r.colors[0].a, 0.0);
        assert_eq!(r.colors[1].a, 1.0);
        s.failed = FailedThreshold::Dim(0.5);
        let r = run(&s, &v);
        assert_eq!(
            (r.colors[0].r, r.colors[0].a),
            (0.05, 1.0),
            "dimmed, still opaque"
        );
        assert_eq!(r.colors[1].r, 0.4, "passing colors untouched");
        s.failed = FailedThreshold::Keep;
        let r = run(&s, &v);
        assert_eq!(r.colors[0].a, 1.0);
        assert!(!r.passed[0], "the pass mask still reports the failure");
        // The NaN intensity is "missing" (step 1) and never reaches the threshold.
        assert_eq!(
            (
                r.diagnostics.failed_threshold,
                r.diagnostics.missing_intensity
            ),
            (1, 1)
        );
    }

    #[test]
    fn a_separate_threshold_column_and_missing_threshold_policy() {
        let mut s = spec();
        s.threshold = Threshold::Above(2.0);
        let intensity = [5.0, 5.0, 5.0];
        let thr = [3.0, 1.0, f64::NAN];
        let inputs = OverlayInputs {
            intensity: &intensity,
            threshold: Some(&thr),
            ..Default::default()
        };
        let r = evaluate_rows(&s, &inputs).unwrap();
        assert_eq!(
            r.passed,
            vec![true, false, false],
            "NaN threshold hides by default"
        );
        s.missing_threshold = MissingThreshold::Show; // AFNI/SUMA behaviour
        let r = evaluate_rows(&s, &inputs).unwrap();
        assert_eq!(r.passed, vec![true, false, true]);
        // Length mismatch is reported.
        let short = [1.0];
        assert!(evaluate_rows(
            &s,
            &OverlayInputs {
                intensity: &intensity,
                threshold: Some(&short),
                ..Default::default()
            }
        )
        .is_err());
    }

    #[test]
    fn fades_scale_alpha_and_styles_are_optional() {
        let mut s = spec();
        s.threshold = Threshold::AbsoluteAbove(4.0);
        s.failed = FailedThreshold::Fade {
            model: FadeModel::Afni {
                curve: FadeCurve::Linear,
                floor: 0.0,
            },
            style: FadeStyle::default(),
        };
        let r = run(&s, &[2.0, 3.99, 4.0, 0.0]);
        assert_eq!(r.colors[0].a, 128.0 / 255.0);
        assert_eq!(r.colors[1].a, 222.0 / 255.0, "the 222 ceiling");
        assert_eq!(r.colors[2].a, 1.0);
        assert_eq!(r.colors[3].a, 0.0, "exactly zero is rejected");
        // The AFNI fade changes opacity only: color channels are untouched.
        assert_eq!(r.colors[0].r, 0.2);
        // Styling (non-AFNI) darkens and desaturates failing colors only.
        let colored = OverlayColors::Panes {
            table: ColorTable::from_rgb_bytes(&[[255, 0, 0], [255, 0, 0]]).unwrap(),
            rule: PaneRule::SumaBanded,
        };
        let mut t = OverlaySpec::new(colored);
        t.intensity_range = RangeSelection::Manual(ColumnRange::new(0.0, 10.0).unwrap());
        t.threshold = Threshold::AbsoluteAbove(4.0);
        t.failed = FailedThreshold::Fade {
            model: FadeModel::Boundary {
                curve: FadeCurve::Linear,
                width: FadeWidth::BoundaryMagnitude,
            },
            style: FadeStyle {
                desaturate: 1.0,
                darken: 0.5,
                boost: 0.0,
                max_alpha: Some(0.6),
            },
        };
        let r = run(&t, &[2.0, 6.0]);
        assert!(r.colors[0].a <= 0.5 + 1e-6 && r.colors[0].a > 0.0);
        assert!(
            r.colors[0].r < 1.0 && r.colors[0].g > 0.0,
            "darker and less saturated"
        );
        assert_eq!(
            r.colors[1],
            Rgba::from_u8(255, 0, 0, 255),
            "passing untouched without boost"
        );
        t.failed = FailedThreshold::Fade {
            model: FadeModel::Boundary {
                curve: FadeCurve::Linear,
                width: FadeWidth::BoundaryMagnitude,
            },
            style: FadeStyle {
                boost: 0.5,
                ..FadeStyle::default()
            },
        };
        // A pure red is already at the edge of the gamut: boost cannot change it.
        assert_eq!(run(&t, &[6.0]).colors[0], Rgba::from_u8(255, 0, 0, 255));
    }

    #[test]
    fn cluster_survivors_hide_unconditionally_even_when_fading() {
        let mut s = spec();
        s.threshold = Threshold::AbsoluteAbove(4.0);
        s.failed = FailedThreshold::Fade {
            model: FadeModel::Suma {
                curve: FadeCurve::Linear,
            },
            style: FadeStyle::default(),
        };
        let v = [8.0, 8.0, 2.0];
        let keep = [true, false, true];
        let r = evaluate_rows(
            &s,
            &OverlayInputs {
                intensity: &v,
                cluster_survivors: Some(&keep),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(r.colors[0].a, 1.0);
        assert_eq!(
            r.colors[1].a, 0.0,
            "passed the threshold but is in a rejected cluster"
        );
        assert!(
            r.colors[2].a > 0.0,
            "a faded sample in a surviving cluster stays faded"
        );
        assert_eq!(r.passed, vec![true, false, false]);
        assert_eq!(r.diagnostics.rejected_by_cluster, 1);
    }

    #[test]
    fn opacity_multiplies_every_alpha_and_is_validated() {
        let mut s = spec();
        s.opacity = 0.5;
        s.missing_color = Rgba::from_u8(1, 2, 3, 255);
        let r = run(&s, &[5.0, f64::NAN]);
        assert_eq!(r.colors[0].a, 0.5);
        assert_eq!(r.colors[1].a, 0.5);
        s.opacity = 1.5;
        assert!(s.validate().is_err());
        s.opacity = f32::NAN;
        assert!(s.validate().is_err());
    }

    #[test]
    fn brightness_modulation_by_a_second_column() {
        let mut s = spec();
        s.brightness = Some(BrightnessModulation {
            range: RangeSelection::Manual(ColumnRange::new(0.0, 1.0).unwrap()),
            scale: [0.0, 1.0],
        });
        let intensity = [10.0, 10.0, 10.0, 10.0];
        let brightness = [0.0, 0.5, 5.0, f64::NAN];
        let r = evaluate_rows(
            &s,
            &OverlayInputs {
                intensity: &intensity,
                brightness: Some(&brightness),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(r.colors[0].r, 0.0, "lowest brightness is black");
        assert_eq!(r.colors[1].r, 0.5);
        assert_eq!(
            r.colors[2].r, 1.0,
            "a manual range clips the brightness column"
        );
        assert_eq!(
            r.colors[3].r, 1.0,
            "a missing brightness leaves the color alone"
        );
        // Brightness values must be supplied when the spec asks for them.
        assert!(evaluate_rows(
            &s,
            &OverlayInputs {
                intensity: &intensity,
                ..Default::default()
            }
        )
        .is_err());
    }

    #[test]
    fn colormap_brightness_factor_scales_and_clamps_the_map() {
        let table = ColorTable::from_rgb_bytes(&[[100, 200, 50], [100, 200, 50]]).unwrap();
        let mut s = OverlaySpec::new(OverlayColors::Panes {
            table,
            rule: PaneRule::SumaBanded,
        });
        s.intensity_range = RangeSelection::Manual(ColumnRange::new(0.0, 1.0).unwrap());
        s.colormap_brightness = 1.5;
        let c = run(&s, &[0.5]).colors[0];
        assert!(
            (c.r - 150.0 / 255.0).abs() < 1e-6 && c.g == 1.0,
            "clamped at 1: {c:?}"
        );
        s.colormap_brightness = 0.0;
        assert!(s.validate().is_err());
        s.colormap_brightness = 2.5;
        assert!(s.validate().is_err());
        s.colormap_brightness = 2.0;
        assert!(s.validate().is_ok());
    }

    #[test]
    fn afni_banded_panes_count_from_the_top_and_differ_from_suma_on_boundaries() {
        // Four panes, lowest value first: A B C D.
        let colors: Vec<[u8; 3]> = vec![[10, 0, 0], [20, 0, 0], [30, 0, 0], [40, 0, 0]];
        let table = ColorTable::from_rgb_bytes(&colors).unwrap();
        let range = ColumnRange::new(0.0, 4.0).unwrap();
        let pick =
            |rule: PaneRule, v: f64| (rule.lookup(&table, v, range).unwrap().to_u8()[0]) as u32;
        // Away from boundaries the two rules agree: A for the lowest quarter...
        for rule in [PaneRule::AfniBanded, PaneRule::SumaBanded] {
            assert_eq!(pick(rule, 0.5), 10, "{rule:?}");
            assert_eq!(pick(rule, 3.5), 40, "{rule:?}");
            assert_eq!(pick(rule, 99.0), 40, "{rule:?} clamps high");
            assert_eq!(pick(rule, -9.0), 10, "{rule:?} clamps low");
        }
        // ...but exactly on a boundary AFNI counts from the top (floor(N/(top-bot)*(top-v))),
        // so v = 1 lands in the pane below SUMA's.
        assert_eq!(pick(PaneRule::SumaBanded, 1.0), 20);
        assert_eq!(pick(PaneRule::AfniBanded, 1.0), 10);
        assert!(PaneRule::AfniBanded
            .lookup(&table, 1.0, ColumnRange { min: 1.0, max: 1.0 })
            .is_err());
    }

    #[test]
    fn suma_panes_hold_the_last_color_over_the_top_pane_and_interpolate_below() {
        let table = ColorTable::from_rgb_bytes(&[[0, 0, 0], [100, 0, 0], [200, 0, 0], [250, 0, 0]])
            .unwrap();
        let range = ColumnRange::new(0.0, 1.0).unwrap();
        let red = |v: f64| {
            PaneRule::SumaInterpolated
                .lookup(&table, v, range)
                .unwrap()
                .to_u8()[0]
        };
        assert_eq!(red(0.0), 0);
        assert_eq!(red(0.125), 50, "halfway between entries 0 and 1");
        assert_eq!(red(0.25), 100, "entry i sits at i/N");
        assert_eq!(red(0.75), 250, "the last entry is reached at (N-1)/N");
        assert_eq!(red(0.9), 250, "and held to the top");
        assert_eq!(red(1.0), 250);
        // A zero-width range gives the middle entry.
        let flat = ColumnRange { min: 3.0, max: 3.0 };
        assert_eq!(
            PaneRule::SumaInterpolated
                .lookup(&table, 3.0, flat)
                .unwrap()
                .to_u8()[0],
            100
        );
    }

    #[test]
    fn direct_mapping_uses_the_truncated_value_as_an_index() {
        let table = ColorTable::from_rgb_bytes(&[[1, 0, 0], [2, 0, 0], [3, 0, 0]]).unwrap();
        let range = ColumnRange { min: 0.0, max: 1.0 };
        let at = |v: f64| PaneRule::Direct.lookup(&table, v, range).unwrap().to_u8()[0];
        assert_eq!(
            (at(0.0), at(1.0), at(1.9), at(2.0), at(50.0), at(-3.0)),
            (1, 2, 2, 3, 3, 1)
        );
        // Direct mapping needs no range, so an all-NaN-free dataset with none still works.
        let mut s = OverlaySpec::new(OverlayColors::Panes {
            table,
            rule: PaneRule::Direct,
        });
        s.intensity_range = RangeSelection::Recorded; // would fail if a range were required
        assert!(evaluate_rows(
            &s,
            &OverlayInputs {
                intensity: &[0.0, 1.0],
                ..Default::default()
            }
        )
        .is_ok());
    }

    #[test]
    fn afni_scales_become_pane_tables_in_the_right_order() {
        let table = AfniColorScale::SpectrumRedToBlue
            .to_color_table(256)
            .unwrap();
        let bytes = AfniColorScale::SpectrumRedToBlue.table(256).unwrap();
        assert_eq!(table.len(), 256);
        // AFNI index 0 (the top of the bar, red) is the LAST entry here.
        assert_eq!(table.colors()[255].to_u8()[..3], bytes[0]);
        assert_eq!(table.colors()[0].to_u8()[..3], bytes[255]);
        // So the highest value is red and the lowest is blue.
        let mut s = OverlaySpec::new(OverlayColors::Panes {
            table,
            rule: PaneRule::AfniBanded,
        });
        s.intensity_range = RangeSelection::Manual(ColumnRange::new(0.0, 1.0).unwrap());
        let r = run(&s, &[1.0, 0.0]);
        assert_eq!(
            r.colors[0].to_u8()[..3],
            bytes[0],
            "highest value: AFNI's top color"
        );
        assert_eq!(
            r.colors[1].to_u8()[..3],
            bytes[255],
            "lowest value: AFNI's bottom color"
        );
    }

    #[test]
    fn label_overlays_color_by_integer_key() {
        let entry = |key: i64, c: [u8; 3]| {
            let [r, g, b] = c.map(|v| f32::from(v) / 255.0);
            LabelEntry {
                key,
                name: format!("k{key}"),
                rgba: Some([r, g, b, 1.0]),
            }
        };
        let table = LabelTable::new(vec![entry(1, [255, 0, 0]), entry(7, [0, 255, 0])]).unwrap();
        let mut s = OverlaySpec::new(OverlayColors::Labels(LabelColorMap::with_defaults(table)));
        s.missing_color = Rgba::from_u8(50, 50, 50, 255);
        let r = run(&s, &[1.0, 7.0, 3.0, 1.5, f64::NAN]);
        assert_eq!(r.colors[0], Rgba::from_u8(255, 0, 0, 255));
        assert_eq!(r.colors[1], Rgba::from_u8(0, 255, 0, 255));
        assert_eq!(r.colors[2], Rgba::TRANSPARENT, "an unlabeled key");
        assert_eq!(
            r.colors[3],
            Rgba::from_u8(50, 50, 50, 255),
            "1.5 is not a label"
        );
        assert_eq!(r.diagnostics.intensity_range, None, "labels need no range");
    }

    fn surface_dataset(indices: Option<Vec<u32>>, values: Vec<f32>, nodes: usize) -> Dataset {
        let domain = Domain::Surface(SurfaceDomain::new(None, nodes).unwrap());
        let col = DataColumn::new("v", ColumnRole::Intensity, ColumnData::Float32(values)).unwrap();
        match indices {
            Some(i) => {
                Dataset::indexed(crate::dataset::DatasetKind::Scalar, domain, i, vec![col]).unwrap()
            }
            None => Dataset::dense(crate::dataset::DatasetKind::Scalar, domain, vec![col]).unwrap(),
        }
    }

    #[test]
    fn dataset_evaluation_expands_sparse_rows_to_the_whole_domain() {
        let ds = surface_dataset(Some(vec![4, 1]), vec![10.0, 0.0], 6);
        let r = evaluate_dataset(&spec(), &ds, &OverlayColumns::new(0), None).unwrap();
        assert_eq!(r.colors.len(), 6);
        assert_eq!(r.colors[4].r, 1.0);
        assert_eq!(r.colors[1].r, 0.0);
        for empty in [0, 2, 3, 5] {
            assert_eq!(
                r.colors[empty],
                Rgba::TRANSPARENT,
                "no row for sample {empty}"
            );
            assert!(!r.passed[empty]);
        }
        assert_eq!(r.diagnostics.rows, 2);
    }

    #[test]
    fn dataset_cluster_survivors_are_per_domain_sample() {
        let ds = surface_dataset(Some(vec![4, 1]), vec![10.0, 10.0], 6);
        let survivors = [false, true, false, false, false, false]; // only sample 1
        let r = evaluate_dataset(&spec(), &ds, &OverlayColumns::new(0), Some(&survivors)).unwrap();
        assert_eq!(r.colors[4].a, 0.0);
        assert_eq!(r.colors[1].a, 1.0);
        assert!(
            evaluate_dataset(&spec(), &ds, &OverlayColumns::new(0), Some(&survivors[..3])).is_err()
        );
    }

    #[test]
    fn dataset_uses_the_recorded_range_and_validates_columns() {
        let domain = Domain::Surface(SurfaceDomain::new(None, 3).unwrap());
        let recorded = crate::column::RecordedRange {
            range: ColumnRange::new(0.0, 100.0).unwrap(),
            min_sample: None,
            max_sample: None,
        };
        let col = DataColumn::new(
            "v",
            ColumnRole::Intensity,
            ColumnData::Float64(vec![0.0, 50.0, 100.0]),
        )
        .unwrap()
        .with_recorded_range(Some(recorded));
        let text = DataColumn::new(
            "t",
            ColumnRole::Generic,
            ColumnData::Text(vec!["a".into(); 3]),
        )
        .unwrap();
        let ds =
            Dataset::dense(crate::dataset::DatasetKind::Scalar, domain, vec![col, text]).unwrap();
        let mut s = OverlaySpec::new(gray_map());
        s.intensity_range = RangeSelection::Recorded;
        let r = evaluate_dataset(&s, &ds, &OverlayColumns::new(0), None).unwrap();
        assert_eq!(r.colors[1].r, 0.5);
        assert!(
            evaluate_dataset(&s, &ds, &OverlayColumns::new(1), None).is_err(),
            "text column"
        );
        assert!(
            evaluate_dataset(&s, &ds, &OverlayColumns::new(9), None).is_err(),
            "no such column"
        );
        let both = OverlayColumns {
            intensity: 0,
            threshold: Some(1),
            brightness: None,
        };
        assert!(
            evaluate_dataset(&s, &ds, &both, None).is_err(),
            "text threshold column"
        );
    }

    #[test]
    fn gpu_parameters_flatten_the_spec() {
        let mut s = spec();
        s.threshold = Threshold::Outside { lo: -2.0, hi: 2.0 };
        s.opacity = 0.75;
        s.show_zero = false;
        s.clip = ClipMode::Hide;
        let range = ColumnRange::new(0.0, 10.0).unwrap();
        let g = s.gpu_params(range, 5).unwrap();
        assert_eq!(g.lookup.len(), 5);
        assert_eq!(g.lookup[0], [0.0, 0.0, 0.0, 1.0]);
        assert_eq!(g.lookup[2][0], 0.5);
        assert_eq!(g.range, [0.0, 10.0]);
        assert_eq!(
            (g.threshold.mode, g.threshold.a, g.threshold.b),
            (4, -2.0, 2.0)
        );
        assert_eq!(
            (g.opacity, g.hide_zero, g.hide_out_of_range),
            (0.75, true, true)
        );
        assert!(s.gpu_params(range, 1).is_err());
        let labels = OverlaySpec::new(OverlayColors::Labels(LabelColorMap::with_defaults(
            LabelTable::new(vec![]).unwrap(),
        )));
        assert!(labels.gpu_params(range, 4).is_err());
    }

    #[test]
    fn invalid_specs_fail_before_touching_data() {
        let mut s = spec();
        s.threshold = Threshold::Between { lo: 3.0, hi: 1.0 };
        assert!(evaluate_rows(
            &s,
            &OverlayInputs {
                intensity: &[1.0],
                ..Default::default()
            }
        )
        .is_err());
        let mut s = spec();
        s.threshold = Threshold::Between { lo: 0.0, hi: 1.0 };
        s.failed = FailedThreshold::Fade {
            model: FadeModel::Afni {
                curve: FadeCurve::Linear,
                floor: 0.0,
            },
            style: FadeStyle::default(),
        };
        assert!(
            s.validate().is_err(),
            "AFNI fade is not defined for Between"
        );
        let mut s = spec();
        s.failed = FailedThreshold::Dim(2.0);
        assert!(s.validate().is_err());
        assert!(ColorTable::new(vec![Rgba::BLACK]).is_err());
    }
}
