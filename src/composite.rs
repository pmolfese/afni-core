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
// Alpha compositing: layering colors over each other. The anatomical underlay is
// the bottom layer; overlay planes (colored statistics, ROI annotations, live
// RGBA overlays from AFNI) are laid over it in order.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `overlay.rs` produces the colors for one overlay plane (`OverlayEvaluation`);
//   this file combines planes. A viewer can instead do the same arithmetic on the
//   GPU; the functions here are the reference.
// * `color.rs` supplies `Rgba`.
//
// STRAIGHT VERSUS PREMULTIPLIED ALPHA (the decision, stated once)
//
// Everything in this crate, including `Rgba`, is STRAIGHT alpha: the stored RGB is
// the color of the object and alpha says how much of it is there; a 50%-opaque red
// is (1, 0, 0, 0.5). Compositing converts to premultiplied form internally, because
// that is where the "over" operator is a plain weighted sum, and converts back.
// `to_premultiplied` / `from_premultiplied` exist for GPU buffers that want
// premultiplied data. Mixing the two conventions is the classic source of dark
// fringes around translucent regions, so the type system keeps them apart: a
// premultiplied color is a different type, `Premultiplied`.
// ---------------------------------------------------------------------------

//! Alpha compositing with straight-alpha colors.

use crate::color::Rgba;
use crate::error::{Error, Result};

/// A color whose RGB has already been multiplied by its alpha. Kept distinct from
/// [`Rgba`] (straight alpha) so the two cannot be confused.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Premultiplied {
    /// Red times alpha.
    pub r: f32,
    /// Green times alpha.
    pub g: f32,
    /// Blue times alpha.
    pub b: f32,
    /// Alpha.
    pub a: f32,
}

/// Convert a straight-alpha color to premultiplied form.
pub fn to_premultiplied(c: Rgba) -> Premultiplied {
    Premultiplied {
        r: c.r * c.a,
        g: c.g * c.a,
        b: c.b * c.a,
        a: c.a,
    }
}

/// Convert back to straight alpha. A fully transparent color has no defined RGB;
/// it becomes transparent black.
pub fn from_premultiplied(p: Premultiplied) -> Rgba {
    if p.a <= 0.0 {
        return Rgba::TRANSPARENT;
    }
    Rgba {
        r: (p.r / p.a).min(1.0),
        g: (p.g / p.a).min(1.0),
        b: (p.b / p.a).min(1.0),
        a: p.a,
    }
}

/// The Porter-Duff "over" operator: `top` laid over `bottom`, both straight alpha.
///
/// In premultiplied terms `out = top + bottom * (1 - top.a)`. Over an opaque
/// bottom this is the familiar `top.a * top + (1 - top.a) * bottom` per channel.
///
/// ```
/// use afni_core::color::Rgba;
/// use afni_core::composite::over;
///
/// let red_half = Rgba::new(1.0, 0.0, 0.0, 0.5)?;
/// let result = over(red_half, Rgba::WHITE);
/// assert_eq!(result.to_u8(), [255, 128, 128, 255]); // half red over white is pink
/// # Ok::<(), afni_core::Error>(())
/// ```
pub fn over(top: Rgba, bottom: Rgba) -> Rgba {
    let (t, b) = (to_premultiplied(top), to_premultiplied(bottom));
    let k = 1.0 - t.a;
    from_premultiplied(Premultiplied {
        r: t.r + b.r * k,
        g: t.g + b.g * k,
        b: t.b + b.b * k,
        a: t.a + b.a * k,
    })
}

/// One plane to composite: its colors (one per sample) and an extra opacity that
/// scales every alpha in the plane.
#[derive(Debug, Clone, Copy)]
pub struct Layer<'a> {
    /// Straight-alpha colors, one per sample.
    pub colors: &'a [Rgba],
    /// Opacity of the whole plane in `[0, 1]`; 1 leaves the colors' own alpha.
    pub opacity: f32,
}

impl<'a> Layer<'a> {
    /// A layer at full opacity.
    pub fn new(colors: &'a [Rgba]) -> Self {
        Self {
            colors,
            opacity: 1.0,
        }
    }
}

/// Composite `layers` over `underlay`, in order: the first layer is laid directly
/// over the underlay, the second over that, and so on (so later layers are on
/// top). All planes must have the underlay's length.
pub fn composite_layers(underlay: &[Rgba], layers: &[Layer<'_>]) -> Result<Vec<Rgba>> {
    for (i, layer) in layers.iter().enumerate() {
        if layer.colors.len() != underlay.len() {
            return Err(Error::LengthMismatch {
                what: format!("layer {i} samples"),
                expected: underlay.len(),
                found: layer.colors.len(),
            });
        }
        if !(0.0..=1.0).contains(&layer.opacity) {
            return Err(Error::InvalidParameter {
                name: format!("layer {i} opacity"),
                reason: format!("{} is outside [0, 1]", layer.opacity),
            });
        }
    }
    let mut out = underlay.to_vec();
    for layer in layers {
        for (pixel, &c) in out.iter_mut().zip(layer.colors) {
            *pixel = over(c.with_alpha(c.a * layer.opacity), *pixel);
        }
    }
    Ok(out)
}

/// An opaque grayscale underlay from anatomical values: each value is placed in
/// `[lo, hi]` (clamped) and becomes gray. Non-finite values are transparent black.
pub fn gray_underlay(values: &[f64], lo: f64, hi: f64) -> Result<Vec<Rgba>> {
    crate::numeric::ensure_finite("underlay low", lo)?;
    crate::numeric::ensure_finite("underlay high", hi)?;
    if hi <= lo {
        return Err(Error::InvalidParameter {
            name: "underlay range".into(),
            reason: format!("high {hi} must exceed low {lo}"),
        });
    }
    Ok(values
        .iter()
        .map(|&v| {
            if v.is_finite() {
                let g = ((v - lo) / (hi - lo)).clamp(0.0, 1.0) as f32;
                Rgba {
                    r: g,
                    g,
                    b: g,
                    a: 1.0,
                }
            } else {
                Rgba::TRANSPARENT
            }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(r: f32, g: f32, b: f32, a: f32) -> Rgba {
        Rgba::new(r, g, b, a).unwrap()
    }

    #[test]
    fn over_an_opaque_bottom_is_the_weighted_sum() {
        let top = c(1.0, 0.0, 0.0, 0.25);
        let bottom = c(0.0, 0.0, 1.0, 1.0);
        let out = over(top, bottom);
        assert_eq!(out.a, 1.0);
        assert!((out.r - 0.25).abs() < 1e-6 && (out.b - 0.75).abs() < 1e-6 && out.g == 0.0);
    }

    #[test]
    fn opaque_top_hides_the_bottom_and_transparent_top_shows_it() {
        let bottom = c(0.2, 0.4, 0.6, 1.0);
        assert_eq!(over(c(1.0, 1.0, 0.0, 1.0), bottom), c(1.0, 1.0, 0.0, 1.0));
        assert_eq!(over(Rgba::TRANSPARENT, bottom), bottom);
        // Over nothing at all: the top stays as it is.
        assert_eq!(
            over(c(1.0, 0.0, 0.0, 0.5), Rgba::TRANSPARENT),
            c(1.0, 0.0, 0.0, 0.5)
        );
        assert_eq!(
            over(Rgba::TRANSPARENT, Rgba::TRANSPARENT),
            Rgba::TRANSPARENT
        );
    }

    #[test]
    fn straight_alpha_has_no_dark_fringe_when_two_translucent_colors_meet() {
        // Two 50% reds: the result is 75% opaque and still pure red (a premultiplied
        // mix-up would darken it).
        let red = c(1.0, 0.0, 0.0, 0.5);
        let out = over(red, red);
        assert!((out.a - 0.75).abs() < 1e-6);
        assert!((out.r - 1.0).abs() < 1e-6 && out.g == 0.0 && out.b == 0.0);
    }

    #[test]
    fn premultiplied_conversion_round_trips_and_handles_transparent() {
        let straight = c(0.8, 0.4, 0.2, 0.5);
        let pre = to_premultiplied(straight);
        assert_eq!((pre.r, pre.g, pre.b, pre.a), (0.4, 0.2, 0.1, 0.5));
        let back = from_premultiplied(pre);
        for (a, b) in back.to_array().iter().zip(straight.to_array()) {
            assert!((a - b).abs() < 1e-6);
        }
        assert_eq!(
            from_premultiplied(to_premultiplied(Rgba::TRANSPARENT)),
            Rgba::TRANSPARENT
        );
    }

    #[test]
    fn over_is_associative_for_straight_alpha_colors() {
        let (a, b, d) = (
            c(1.0, 0.0, 0.0, 0.3),
            c(0.0, 1.0, 0.0, 0.6),
            c(0.0, 0.0, 1.0, 0.9),
        );
        let left = over(over(a, b), d);
        let right = over(a, over(b, d));
        for (x, y) in left.to_array().iter().zip(right.to_array()) {
            assert!((x - y).abs() < 1e-5, "{left:?} vs {right:?}");
        }
    }

    #[test]
    fn layers_stack_in_order_with_per_layer_opacity() {
        let under = vec![Rgba::BLACK; 2];
        let red = vec![c(1.0, 0.0, 0.0, 1.0), Rgba::TRANSPARENT];
        let green = vec![c(0.0, 1.0, 0.0, 1.0); 2];
        // Red first, then green on top at half opacity.
        let out = composite_layers(
            &under,
            &[
                Layer::new(&red),
                Layer {
                    colors: &green,
                    opacity: 0.5,
                },
            ],
        )
        .unwrap();
        assert!(
            (out[0].r - 0.5).abs() < 1e-6 && (out[0].g - 0.5).abs() < 1e-6,
            "{:?}",
            out[0]
        );
        assert!(
            (out[1].g - 0.5).abs() < 1e-6 && out[1].r == 0.0,
            "the transparent red shows the black"
        );
        // No layers: the underlay unchanged.
        assert_eq!(composite_layers(&under, &[]).unwrap(), under);
        // Bad inputs.
        assert!(composite_layers(&under, &[Layer::new(&red[..1])]).is_err());
        assert!(composite_layers(
            &under,
            &[Layer {
                colors: &red,
                opacity: 1.5
            }]
        )
        .is_err());
    }

    #[test]
    fn gray_underlay_clamps_and_marks_missing() {
        let u = gray_underlay(&[0.0, 5.0, 10.0, 20.0, f64::NAN], 0.0, 10.0).unwrap();
        assert_eq!(u[0], Rgba::BLACK);
        assert_eq!(u[1].r, 0.5);
        assert_eq!(u[3], Rgba::WHITE, "clamped");
        assert_eq!(u[4], Rgba::TRANSPARENT);
        assert!(gray_underlay(&[1.0], 1.0, 1.0).is_err());
    }
}
