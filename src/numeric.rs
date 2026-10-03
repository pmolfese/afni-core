// PUBLIC DOMAIN NOTICE
//
// This file is part of afni-core, which was written by employees of the United
// States Government (National Institutes of Health) as part of their official
// duties. Under Title 17, Section 105 of the United States Code it is a
// "United States Government Work" and is not subject to copyright protection
// in the United States. Outside the United States, rights are waived under
// CC0 1.0 Universal. See the LICENSE file at the crate root.
//
// ---------------------------------------------------------------------------
// WHAT THIS FILE IS
//
// The numeric conventions every other module in `afni-core` follows:
//
//   1. Precision:   `f64` for math, `f32` for stored display/mesh buffers.
//   2. Indices:     checked conversion between file integers and `usize`/`u32`.
//   3. Non-finite:  one explicit policy for NaN and +/-infinity.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * Phase 1 (datasets/domains) uses `checked_index` to validate sparse node
//   lists against a domain, and `usize_to_u32` when building index buffers.
// * Phase 2-3 (statistics, FDR) use `ensure_finite` and `NonFinitePolicy` so a
//   NaN statistic is an error or a skipped sample, never a made-up p-value.
// * Phase 4-5 (colors, overlays) use `narrow_to_f32` when producing the f32
//   buffers a viewer uploads to the GPU.
// * Errors come from `error.rs`.
// ---------------------------------------------------------------------------

//! Shared numeric conventions: precision, checked indices, and NaN/Inf policy.
//!
//! # The three conventions
//!
//! **Precision.** Statistical math (p-values, critical values, interpolation of
//! FDR curves) is done in [`f64`]. Data that is merely *displayed* (mesh
//! coordinates, per-vertex colors, uploaded intensity buffers) is stored as
//! [`f32`], halving memory and matching what GPUs consume. Convert at the
//! boundary, in one place, with [`widen_to_f64`] and [`narrow_to_f32`], so the
//! precision loss is visible in the code.
//!
//! **Indices.** Files store indices as signed or unsigned integers of various
//! widths; Rust slices are indexed by [`usize`]; GPU index buffers are
//! usually `u32`. Never use a bare `as` cast between them: `-1i64 as usize`
//! silently becomes a gigantic number. Use [`checked_index`] and
//! [`usize_to_u32`], which return an error instead.
//!
//! **NaN and infinity.** `NaN` compares unequal to everything (even itself) and
//! poisons arithmetic, so it must be handled deliberately. Every operation
//! that can meet a non-finite value takes or documents a [`NonFinitePolicy`];
//! there is no hidden default that replaces such values with zero.

use crate::error::{Error, Result};

// ---------------------------------------------------------------------------
// Precision conversion
// ---------------------------------------------------------------------------

/// Widen a stored `f32` to `f64` for calculation.
///
/// This is lossless: every `f32` is exactly representable as an `f64`. The
/// function exists mainly so call sites read as an intentional step
/// ("data is now in math precision") rather than an incidental `as` cast.
///
/// ```
/// assert_eq!(afni_core::numeric::widen_to_f64(0.5_f32), 0.5_f64);
/// ```
#[inline]
pub fn widen_to_f64(x: f32) -> f64 {
    f64::from(x)
}

/// Narrow a calculated `f64` to `f32` for a display/mesh buffer.
///
/// This **rounds** to the nearest representable `f32` and is therefore lossy.
/// Values too large for `f32` become infinity, which is the standard IEEE
/// behavior of an `as` cast; NaN stays NaN. Callers that cannot tolerate
/// either should check with [`ensure_finite`] first.
///
/// ```
/// assert_eq!(afni_core::numeric::narrow_to_f32(0.5_f64), 0.5_f32);
/// ```
#[inline]
pub fn narrow_to_f32(x: f64) -> f32 {
    // `as` is exactly the rounding conversion we want here; this is the one
    // sanctioned place in the crate where a float narrowing cast appears.
    x as f32
}

// ---------------------------------------------------------------------------
// Checked index conversion
// ---------------------------------------------------------------------------

/// Convert a signed file index into a valid `usize` index into `len` samples.
///
/// Returns [`Error::IndexOutOfRange`] if `index` is negative or `>= len`.
/// The check is done on the signed value *before* converting, so negative
/// numbers cannot wrap around.
///
/// ```
/// use afni_core::numeric::checked_index;
///
/// assert_eq!(checked_index(2, 5).unwrap(), 2);
/// assert!(checked_index(5, 5).is_err());  // one past the end
/// assert!(checked_index(-1, 5).is_err()); // negative
/// ```
pub fn checked_index(index: i64, len: usize) -> Result<usize> {
    // `usize::try_from` fails for negative numbers (and, on 32-bit targets,
    // for values above `u32::MAX`), which is exactly the first half of the
    // check we need.
    match usize::try_from(index) {
        Ok(i) if i < len => Ok(i),
        // Either negative/too large for usize, or simply >= len.
        _ => Err(Error::IndexOutOfRange { index, len }),
    }
}

/// Convert a `usize` to a `u32`, as needed for GPU index buffers and many
/// on-disk formats, failing instead of truncating.
///
/// ```
/// use afni_core::numeric::usize_to_u32;
///
/// assert_eq!(usize_to_u32(7).unwrap(), 7);
/// # #[cfg(target_pointer_width = "64")]
/// assert!(usize_to_u32(u32::MAX as usize + 1).is_err());
/// ```
pub fn usize_to_u32(value: usize) -> Result<u32> {
    u32::try_from(value).map_err(|_| Error::IndexOverflow { value })
}

// ---------------------------------------------------------------------------
// Non-finite handling
// ---------------------------------------------------------------------------

/// Return `value` unchanged if it is finite, otherwise an [`Error::NonFinite`]
/// naming `what` was wrong.
///
/// Use this for *parameters and statistics that must be real numbers*, such as
/// a degrees-of-freedom value or a threshold. For *data samples*, where a
/// non-finite entry may legitimately mean "missing", use [`NonFinitePolicy`].
///
/// ```
/// use afni_core::numeric::ensure_finite;
///
/// assert_eq!(ensure_finite("dof", 10.0).unwrap(), 10.0);
/// assert!(ensure_finite("dof", f64::NAN).is_err());
/// assert!(ensure_finite("dof", f64::INFINITY).is_err());
/// ```
pub fn ensure_finite(what: &str, value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::NonFinite {
            what: what.to_owned(),
            value,
        })
    }
}

/// What to do with a NaN or infinite *data sample*.
///
/// AFNI datasets contain NaN for "no data" (for example outside a mask), so a
/// blanket error would be unusable, and a blanket replacement with zero would
/// be scientifically wrong. Callers therefore choose, explicitly:
///
/// | Policy | Finite sample | Non-finite sample |
/// |--------|---------------|-------------------|
/// | [`Reject`](Self::Reject) | kept | error |
/// | [`Skip`](Self::Skip) | kept | dropped (reported as "missing") |
/// | [`Propagate`](Self::Propagate) | kept | kept as-is |
///
/// Apply a policy to one sample with [`NonFinitePolicy::check`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NonFinitePolicy {
    /// Fail with [`Error::NonFinite`]. The safe choice for inputs that should
    /// never contain missing values.
    Reject,
    /// Treat the sample as missing: it takes no part in the calculation and the
    /// output for it is "no value".
    Skip,
    /// Pass the value through untouched; downstream arithmetic then follows
    /// ordinary IEEE rules (NaN in, NaN out).
    Propagate,
}

impl NonFinitePolicy {
    /// Apply this policy to one sample.
    ///
    /// The return type encodes the three outcomes:
    ///
    /// * `Ok(Some(x))` - use `x`.
    /// * `Ok(None)` - the sample is missing; skip it.
    /// * `Err(_)` - the policy is [`Reject`](Self::Reject) and `x` is not finite.
    ///
    /// `what` names the quantity for the error message.
    ///
    /// ```
    /// use afni_core::numeric::NonFinitePolicy::{Propagate, Reject, Skip};
    ///
    /// assert_eq!(Reject.check("t", 1.5).unwrap(), Some(1.5));
    /// assert!(Reject.check("t", f64::NAN).is_err());
    /// assert_eq!(Skip.check("t", f64::NAN).unwrap(), None);
    /// assert!(Propagate.check("t", f64::NAN).unwrap().unwrap().is_nan());
    /// ```
    pub fn check(self, what: &str, x: f64) -> Result<Option<f64>> {
        // Finite values are always usable, whatever the policy.
        if x.is_finite() {
            return Ok(Some(x));
        }
        match self {
            NonFinitePolicy::Reject => ensure_finite(what, x).map(Some),
            NonFinitePolicy::Skip => Ok(None),
            NonFinitePolicy::Propagate => Ok(Some(x)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checked_index_accepts_only_valid_range() {
        assert_eq!(checked_index(0, 3).unwrap(), 0);
        assert_eq!(checked_index(2, 3).unwrap(), 2);
        // Boundaries: first invalid on each side.
        assert_eq!(
            checked_index(3, 3),
            Err(Error::IndexOutOfRange { index: 3, len: 3 })
        );
        assert_eq!(
            checked_index(-1, 3),
            Err(Error::IndexOutOfRange { index: -1, len: 3 })
        );
        // An empty domain has no valid index at all.
        assert!(checked_index(0, 0).is_err());
        // Extreme values must not wrap.
        assert!(checked_index(i64::MIN, 10).is_err());
        assert!(checked_index(i64::MAX, 10).is_err());
    }

    #[test]
    #[cfg(target_pointer_width = "64")]
    fn usize_to_u32_detects_overflow() {
        assert_eq!(usize_to_u32(u32::MAX as usize).unwrap(), u32::MAX);
        let too_big = u32::MAX as usize + 1;
        assert_eq!(
            usize_to_u32(too_big),
            Err(Error::IndexOverflow { value: too_big })
        );
    }

    #[test]
    fn ensure_finite_rejects_nan_and_both_infinities() {
        assert!(ensure_finite("x", 0.0).is_ok());
        assert!(ensure_finite("x", -0.0).is_ok());
        assert!(ensure_finite("x", f64::NAN).is_err());
        assert!(ensure_finite("x", f64::INFINITY).is_err());
        assert!(ensure_finite("x", f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn policies_treat_every_non_finite_kind_alike() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(NonFinitePolicy::Reject.check("x", bad).is_err());
            assert_eq!(NonFinitePolicy::Skip.check("x", bad).unwrap(), None);
            let kept = NonFinitePolicy::Propagate.check("x", bad).unwrap();
            // NaN != NaN, so compare bit patterns to prove "unchanged".
            assert_eq!(kept.map(f64::to_bits), Some(bad.to_bits()));
        }
        // Finite values pass under every policy.
        for p in [
            NonFinitePolicy::Reject,
            NonFinitePolicy::Skip,
            NonFinitePolicy::Propagate,
        ] {
            assert_eq!(p.check("x", 2.5).unwrap(), Some(2.5));
        }
    }

    #[test]
    fn precision_round_trip_is_exact_for_f32_values() {
        // f32 -> f64 -> f32 never changes an f32 value (including NaN-ness).
        for x in [0.0_f32, -1.5, 1e-30, 3.4e38, f32::MIN_POSITIVE] {
            assert_eq!(narrow_to_f32(widen_to_f64(x)), x);
        }
        assert!(narrow_to_f32(widen_to_f64(f32::NAN)).is_nan());
        // Narrowing an out-of-range f64 overflows to infinity, as documented.
        assert!(narrow_to_f32(1e300).is_infinite());
    }
}
