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
// The generic signal operations behind AFNI's resting-state correlation tools,
// with no knowledge of surfaces, seeds or files:
//
// * removing a mean, linear or quadratic trend (`Detrend`),
// * Legendre polynomial regressors (`legendre_basis`, AFNI's `THD_build_polyref`),
// * projecting regressors ("orts") out of a series (`OrtProjector`),
// * L2 normalization (`normalize_l2`, AFNI's `THD_normalize`),
// * an FFT bandpass with orts, exactly as `THD_bandpass_vectors` does it
//   (`bandpass_vectors`), including the number of degrees of freedom it removes.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `instacorr.rs` builds the seed-correlation pipeline from these pieces; this
//   file knows nothing about correlation.
// * There is no dependency: the FFT here is a small radix-2 / Bluestein
//   implementation (any length), written for clarity and checked against AFNI's
//   own `1dBandpass` in tests/signal_conformance.rs.
//
// WHAT IS EXACT, AND WHAT IS NOT
//
// The arithmetic follows AFNI step by step: the same FFT length (the series length
// rounded up to an even number), the same band indices (computed in 32-bit float
// like AFNI, ties to even), the same edge taper (0.5, or 0.05 when orts are removed),
// the 0 and Nyquist frequencies always dropped, filtered orts projected out, and the
// same count of removed degrees of freedom. AFNI computes in `f32` and this crate in
// `f64`, so results agree to about 1e-5 relative, not bit for bit.
//
// Differences on purpose (see the roadmap log):
// * Bad input is an error. AFNI quietly substitutes a time step of 1.0 for a
//   non-positive one and, for `ftop <= fbot`, assigns the value to `fbot` (a typo in
//   thd_bandpass.c that turns the call into a low-pass).
// * Orts are projected out by orthonormalizing them (modified Gram-Schmidt, dropping
//   directions smaller than 1e-8 of the column) instead of AFNI's SVD pseudo-inverse
//   with the same 1e-8 cutoff; the two agree for any sensible set of regressors.
// ---------------------------------------------------------------------------

//! Generic signal operations: detrending, orts, power spectra, and FFT bandpass.

use crate::error::{Error, Result};
use crate::numeric::{ensure_finite, NonFinitePolicy};

/// AFNI's `ICOR_MAX_FTOP`: a top frequency at or above this means "no upper limit"
/// (a high-pass filter).
pub const MAX_FTOP_HZ: f64 = 99_999.0;

/// The shortest series AFNI will bandpass.
pub const MIN_BANDPASS_SAMPLES: usize = 9;

// ---------------------------------------------------------------------------
// Legendre polynomials and detrending
// ---------------------------------------------------------------------------

/// The Legendre polynomial of order `m` at `x` (the recurrence
/// `k P_k = (2k-1) x P_{k-1} - (k-1) P_{k-2}`).
pub fn legendre(x: f64, m: usize) -> f64 {
    match m {
        0 => 1.0,
        1 => x,
        _ => {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=m {
                let kf = k as f64;
                let next = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                p0 = p1;
                p1 = next;
            }
            p1
        }
    }
}

/// `count` Legendre polynomial regressors of orders `0..count` sampled at `len`
/// equally spaced points over `[-1, 1]` (AFNI's `THD_build_polyref`). Errors unless
/// `count >= 1` and `len > count`.
pub fn legendre_basis(count: usize, len: usize) -> Result<Vec<Vec<f64>>> {
    if count < 1 || len <= count {
        return Err(Error::InvalidParameter {
            name: "Legendre basis".into(),
            reason: format!("{count} regressors need more than {count} samples, got {len}"),
        });
    }
    let step = 2.0 / (len as f64 - 1.0);
    Ok((0..count)
        .map(|order| {
            (0..len)
                .map(|k| legendre(step * k as f64 - 1.0, order))
                .collect()
        })
        .collect())
}

/// What trend to remove before filtering (AFNI's `qdet`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Detrend {
    /// Nothing.
    None,
    /// Subtract the mean.
    Mean,
    /// Subtract the least-squares line.
    Linear,
    /// Subtract the least-squares quadratic.
    #[default]
    Quadratic,
}

impl Detrend {
    /// The polynomial degree removed, if any.
    fn degree(self) -> Option<usize> {
        match self {
            Self::None => Option::None,
            Self::Mean => Some(0),
            Self::Linear => Some(1),
            Self::Quadratic => Some(2),
        }
    }

    /// The shortest series AFNI detrends: below it the series is left alone (a mean
    /// needs 2 points, a line 3, a quadratic 4).
    fn min_len(self) -> usize {
        match self {
            Self::None => usize::MAX,
            Self::Mean => 2,
            Self::Linear => 3,
            Self::Quadratic => 4,
        }
    }

    /// Remove the trend from `values` in place, as AFNI's `THD_const_detrend`,
    /// `THD_linear_detrend` and `THD_quadratic_detrend` do (a least-squares
    /// polynomial fit in the sample index). Series shorter than the minimum are
    /// left untouched, like AFNI.
    pub fn apply(self, values: &mut [f64]) {
        let Some(degree) = self.degree() else { return };
        if values.len() < self.min_len() {
            return;
        }
        // Projecting out an orthonormal basis of the polynomials is the same least
        // squares fit (and stays accurate for long series).
        let basis = orthonormal_basis(
            &legendre_basis(degree + 1, values.len()).expect("length checked above"),
            1e-8,
        );
        project_out(values, &basis);
    }
}

// ---------------------------------------------------------------------------
// Orts
// ---------------------------------------------------------------------------

/// Modified Gram-Schmidt: an orthonormal basis for the span of `columns`. A column
/// whose part outside the span so far is smaller than `relative_tolerance` times its
/// own length adds nothing and is dropped.
fn orthonormal_basis(columns: &[Vec<f64>], relative_tolerance: f64) -> Vec<Vec<f64>> {
    let mut basis: Vec<Vec<f64>> = Vec::new();
    // Like a pseudo-inverse, judge "negligible" against the LARGEST column, so a
    // column that a filter has reduced to rounding noise is dropped instead of being
    // blown up to unit length.
    let scale = columns.iter().map(|c| norm(c)).fold(0.0_f64, f64::max);
    for column in columns {
        let original = scale;
        if original == 0.0 {
            continue;
        }
        let mut v = column.clone();
        // Subtract the parts along the basis so far; do it twice, which keeps the
        // basis orthogonal to machine precision even for awkward columns.
        for _ in 0..2 {
            project_out(&mut v, &basis);
        }
        let n = norm(&v);
        if n > relative_tolerance * original {
            v.iter_mut().for_each(|x| *x /= n);
            basis.push(v);
        }
    }
    basis
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// Subtract from `values` its components along each (orthonormal) basis vector.
fn project_out(values: &mut [f64], basis: &[Vec<f64>]) {
    for b in basis {
        let w: f64 = values.iter().zip(b).map(|(x, y)| x * y).sum();
        values.iter_mut().zip(b).for_each(|(x, y)| *x -= w * y);
    }
}

/// Removes the part of a series that lies in the span of some regressors (orts):
/// ordinary least squares, `y - X (X'X)^-1 X' y`.
#[derive(Debug, Clone, PartialEq)]
pub struct OrtProjector {
    basis: Vec<Vec<f64>>,
    len: usize,
}

impl OrtProjector {
    /// Build from regressors that all have `len` samples. Regressors that are zero or
    /// (nearly) combinations of earlier ones add nothing and are ignored.
    pub fn new(orts: &[Vec<f64>], len: usize) -> Result<Self> {
        for (i, ort) in orts.iter().enumerate() {
            if ort.len() != len {
                return Err(Error::LengthMismatch {
                    what: format!("ort {i}"),
                    expected: len,
                    found: ort.len(),
                });
            }
            if ort.iter().any(|v| !v.is_finite()) {
                return Err(Error::NonFinite {
                    what: format!("ort {i}"),
                    value: f64::NAN,
                });
            }
        }
        Ok(Self {
            basis: orthonormal_basis(orts, 1e-8),
            len,
        })
    }

    /// Number of independent directions that will be removed.
    pub fn rank(&self) -> usize {
        self.basis.len()
    }

    /// Remove the regressors' part from `values` (which must have the same length).
    pub fn apply(&self, values: &mut [f64]) {
        debug_assert_eq!(values.len(), self.len);
        project_out(values, &self.basis);
    }
}

// ---------------------------------------------------------------------------
// Normalization
// ---------------------------------------------------------------------------

/// Scale `values` to unit L2 length and return the scale factor, as AFNI's
/// `THD_normalize` does. A vector whose squared length is `<= 1e-20` is left
/// unchanged and `None` is returned.
pub fn normalize_l2(values: &mut [f64]) -> Option<f64> {
    let sumsq: f64 = values.iter().map(|x| x * x).sum();
    // Too small, or not a number (NaN or infinite): nothing sensible to scale.
    if !sumsq.is_finite() || sumsq <= 1e-20 {
        return None;
    }
    let factor = 1.0 / sumsq.sqrt();
    values.iter_mut().for_each(|x| *x *= factor);
    Some(factor)
}

// ---------------------------------------------------------------------------
// FFT (radix-2 and Bluestein, any length)
// ---------------------------------------------------------------------------

/// A complex number (just what the FFT needs).
#[derive(Debug, Clone, Copy, PartialEq)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    const ZERO: Self = Self { re: 0.0, im: 0.0 };

    fn mul(self, o: Self) -> Self {
        Self {
            re: self.re * o.re - self.im * o.im,
            im: self.re * o.im + self.im * o.re,
        }
    }

    fn conj(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    fn scale(self, s: f64) -> Self {
        Self {
            re: self.re * s,
            im: self.im * s,
        }
    }
}

/// In-place iterative radix-2 FFT; `data.len()` must be a power of two. `inverse`
/// flips the sign of the exponent (no 1/n scaling).
fn fft_pow2(data: &mut [Complex], inverse: bool) {
    let n = data.len();
    debug_assert!(n.is_power_of_two());
    // Reorder into bit-reversed order.
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            data.swap(i, j);
        }
    }
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut len = 2;
    while len <= n {
        let angle = sign * 2.0 * std::f64::consts::PI / len as f64;
        for start in (0..n).step_by(len) {
            for k in 0..len / 2 {
                let w = Complex {
                    re: (angle * k as f64).cos(),
                    im: (angle * k as f64).sin(),
                };
                let a = data[start + k];
                let b = data[start + k + len / 2].mul(w);
                data[start + k] = Complex {
                    re: a.re + b.re,
                    im: a.im + b.im,
                };
                data[start + k + len / 2] = Complex {
                    re: a.re - b.re,
                    im: a.im - b.im,
                };
            }
        }
        len <<= 1;
    }
}

/// FFT of any length: radix-2 when the length is a power of two, otherwise
/// Bluestein's algorithm (a convolution done with power-of-two FFTs).
fn fft(data: &mut [Complex], inverse: bool) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    if n.is_power_of_two() {
        fft_pow2(data, inverse);
        return;
    }
    let sign = if inverse { 1.0 } else { -1.0 };
    // chirp[k] = exp(sign * i * pi * k^2 / n); k^2 is reduced mod 2n to keep the
    // angle small and the cosines accurate.
    let chirp: Vec<Complex> = (0..n)
        .map(|k| {
            let k2 = ((k as u128 * k as u128) % (2 * n as u128)) as f64;
            let a = sign * std::f64::consts::PI * k2 / n as f64;
            Complex {
                re: a.cos(),
                im: a.sin(),
            }
        })
        .collect();
    let m = (2 * n - 1).next_power_of_two();
    let mut a = vec![Complex::ZERO; m];
    let mut b = vec![Complex::ZERO; m];
    for k in 0..n {
        a[k] = data[k].mul(chirp[k]);
    }
    b[0] = chirp[0].conj();
    for k in 1..n {
        b[k] = chirp[k].conj();
        b[m - k] = chirp[k].conj();
    }
    fft_pow2(&mut a, false);
    fft_pow2(&mut b, false);
    for k in 0..m {
        a[k] = a[k].mul(b[k]);
    }
    fft_pow2(&mut a, true);
    let inv_m = 1.0 / m as f64;
    for k in 0..n {
        data[k] = a[k].scale(inv_m).mul(chirp[k]);
    }
}

/// Which bins a real-input [`power_spectrum`] returns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpectrumLayout {
    /// All `nfft` bins, including the redundant negative-frequency half.
    Full,
    /// Bins `0..=floor(nfft / 2)`.
    ///
    /// Interior bins are not doubled. This is exactly the corresponding prefix
    /// of [`Full`](Self::Full), which keeps the operation unambiguous for
    /// callers that want to combine particular bins themselves.
    OneSided,
}

/// Scaling applied to squared FFT magnitudes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SpectrumNormalization {
    /// Raw squared magnitude of the unnormalized forward transform.
    Raw,
    /// Divide every squared magnitude by the FFT length.
    ///
    /// For [`SpectrumLayout::Full`], the sum of the returned powers equals the
    /// sum of squares of the zero-padded input (Parseval scaling).
    DivideByFftLength,
}

/// Explicit choices for [`power_spectrum`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PowerSpectrumOptions {
    /// Transform length. `None` uses the input length; `Some(n)` zero-pads to
    /// `n` and must not be shorter than the input.
    pub fft_len: Option<usize>,
    /// Whether to return the full transform or its nonnegative-frequency half.
    pub layout: SpectrumLayout,
    /// Scaling applied after squared magnitudes are computed.
    pub normalization: SpectrumNormalization,
    /// Policy for NaN and infinity in the ordered input series.
    ///
    /// `Reject` reports the first bad sample and `Propagate` lets ordinary FFT
    /// arithmetic carry it into the result. `Skip` cannot preserve temporal
    /// positions and returns an error if a non-finite sample is encountered;
    /// replace such samples explicitly before calling when that is intended.
    pub non_finite: NonFinitePolicy,
}

impl PowerSpectrumOptions {
    /// Construct a spectrum request while making every convention explicit.
    pub const fn new(
        fft_len: Option<usize>,
        layout: SpectrumLayout,
        normalization: SpectrumNormalization,
        non_finite: NonFinitePolicy,
    ) -> Self {
        Self {
            fft_len,
            layout,
            normalization,
            non_finite,
        }
    }
}

/// Full or one-sided power spectrum of a real-valued series.
///
/// The input is copied into a zero-filled transform buffer and passed through
/// the same arbitrary-length FFT used by [`bandpass_vectors`]. The forward FFT
/// itself is unnormalized. [`PowerSpectrumOptions`] makes padding, returned
/// bins, output scaling, and non-finite handling visible at the call site.
///
/// This function computes only squared magnitudes. The private complex FFT and
/// its inverse remain implementation details so callers cannot accidentally
/// depend on their internal representation or scaling convention.
pub fn power_spectrum(values: &[f64], options: PowerSpectrumOptions) -> Result<Vec<f64>> {
    if values.is_empty() {
        return Err(Error::Empty("power spectrum input".into()));
    }
    let fft_len = options.fft_len.unwrap_or(values.len());
    if fft_len < values.len() {
        return Err(Error::InvalidParameter {
            name: "FFT length".into(),
            reason: format!(
                "must be at least the {} input samples, got {fft_len}",
                values.len()
            ),
        });
    }

    let mut transformed = vec![Complex::ZERO; fft_len];
    for (sample, &value) in values.iter().enumerate() {
        transformed[sample].re = if value.is_finite() {
            // Keep the hot path allocation-free; the diagnostic string below
            // is needed only for an actual bad sample.
            value
        } else {
            match options.non_finite {
                NonFinitePolicy::Reject => {
                    return Err(Error::NonFinite {
                        what: format!("power spectrum sample {sample}"),
                        value,
                    });
                }
                NonFinitePolicy::Skip => {
                    return Err(Error::InvalidParameter {
                    name: "power spectrum non-finite policy".into(),
                    reason: "Skip cannot remove a sample from an ordered time series; replace it explicitly or choose Reject/Propagate".into(),
                });
                }
                NonFinitePolicy::Propagate => value,
            }
        };
    }
    fft(&mut transformed, false);

    let bins = match options.layout {
        SpectrumLayout::Full => fft_len,
        SpectrumLayout::OneSided => fft_len / 2 + 1,
    };
    let scale = match options.normalization {
        SpectrumNormalization::Raw => 1.0,
        SpectrumNormalization::DivideByFftLength => 1.0 / fft_len as f64,
    };
    Ok(transformed
        .into_iter()
        .take(bins)
        .map(|value| (value.re * value.re + value.im * value.im) * scale)
        .collect())
}

// ---------------------------------------------------------------------------
// Bandpass
// ---------------------------------------------------------------------------

/// A pass band: the time step between samples and the lowest and highest frequency
/// kept, in Hz.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BandSpec {
    /// Seconds between samples. Must be positive.
    pub dt: f64,
    /// Lowest frequency kept (`0` for a low-pass). Must be `>= 0`.
    pub fbot: f64,
    /// Highest frequency kept. Must be above `fbot`; a value of
    /// [`MAX_FTOP_HZ`] or more means no upper limit (a high-pass), and one above the
    /// Nyquist frequency is clipped to it, as in AFNI.
    pub ftop: f64,
}

/// How a pass band maps onto the FFT of a series of a given length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BandPlan {
    /// FFT length: the series length rounded up to an even number.
    pub nfft: usize,
    /// First FFT bin kept (0 means no low cut).
    pub jbot: usize,
    /// Last FFT bin kept (`nfft / 2 - 1` means no high cut).
    pub jtop: usize,
}

impl BandPlan {
    /// Whether the band actually removes frequencies (otherwise only the mean and
    /// Nyquist bins would go, and AFNI skips the FFT entirely).
    pub fn filters(&self) -> bool {
        self.jbot > 0 || self.jtop < self.nfft / 2 - 1
    }
}

/// Round half to even, in `f32` (C's `rintf`). Stable Rust only has this for `f64`
/// from 1.77, and AFNI does the band-index arithmetic in `f32`.
fn rint_f32(x: f32) -> f32 {
    let r = x.round();
    if (x - x.trunc()).abs() == 0.5 {
        2.0 * (x / 2.0).round()
    } else {
        r
    }
}

/// Round half to even (C's `rint`).
fn rint_f64(x: f64) -> f64 {
    let r = x.round();
    if (x - x.trunc()).abs() == 0.5 {
        2.0 * (x / 2.0).round()
    } else {
        r
    }
}

/// Work out the FFT length and the kept bins for a series of `len` samples, the way
/// `THD_bandpass_vectors` does. Errors for a series shorter than 9, a non-positive
/// time step or a band with `fbot < 0` or `ftop <= fbot`.
pub fn band_plan(len: usize, band: &BandSpec) -> Result<BandPlan> {
    if len < MIN_BANDPASS_SAMPLES {
        return Err(Error::InvalidParameter {
            name: "series length".into(),
            reason: format!("bandpass needs at least {MIN_BANDPASS_SAMPLES} samples, got {len}"),
        });
    }
    ensure_finite("time step", band.dt)?;
    ensure_finite("low frequency", band.fbot)?;
    if band.dt <= 0.0 {
        return Err(Error::InvalidParameter {
            name: "time step".into(),
            reason: format!("{} is not positive", band.dt),
        });
    }
    if band.fbot < 0.0 || band.ftop.is_nan() || band.ftop <= band.fbot {
        return Err(Error::InvalidParameter {
            name: "pass band".into(),
            reason: format!(
                "needs 0 <= fbot < ftop, got fbot {} and ftop {}",
                band.fbot, band.ftop
            ),
        });
    }
    let nfft = len + len % 2;
    let nby2 = nfft / 2;
    let nhalf = nby2 - 1;
    // AFNI does these three lines in 32-bit float.
    let df = 1.0_f32 / (nfft as f32 * band.dt as f32);
    let qbot = rint_f32(band.fbot as f32 / df);
    let qtop = rint_f32(band.ftop as f32 / df);
    let mut jbot = if qbot < nhalf as f32 {
        qbot as usize
    } else {
        0
    };
    let mut jtop = if qtop < nhalf as f32 {
        qtop as usize
    } else {
        nhalf
    };
    if band.ftop >= MAX_FTOP_HZ {
        jtop = nhalf;
    }
    if jtop >= nby2 {
        jtop = nhalf;
    }
    if jbot > jtop {
        // The band collapsed (both cutoffs round to the same bin): AFNI keeps it all.
        jbot = 0;
        jtop = nhalf;
    }
    Ok(BandPlan { nfft, jbot, jtop })
}

/// The number of dimensions left after bandpassing (twice the number of kept bins,
/// at most the series length): `THD_bandpass_remain_dim`.
pub fn bandpass_remaining_dimension(len: usize, band: &BandSpec) -> Result<usize> {
    let plan = band_plan(len, band)?;
    // This function (unlike the filter) treats a band of one or two bins as "ignore
    // the filter" and keeps everything.
    let (jbot, jtop) = if plan.jbot + 1 >= plan.jtop {
        (0, plan.nfft / 2 - 1)
    } else {
        (plan.jbot, plan.jtop)
    };
    Ok((2 * (jtop - jbot + 1)).min(len))
}

/// What [`bandpass_vectors`] did.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BandpassReport {
    /// The number of linear dimensions projected out: AFNI's return value. This is
    /// what to subtract from the degrees of freedom of anything computed afterwards.
    pub removed_dof: usize,
    /// The FFT length used (0 if no filtering was done).
    pub nfft: usize,
}

/// Detrend, bandpass and remove orts from a set of equal-length series, in place;
/// AFNI's `THD_bandpass_vectors`.
///
/// In order: (1) the `detrend` is removed from each series; (2) if `band` is given
/// and removes anything, each series is FFT-filtered, with the 0 and Nyquist bins
/// dropped and the band edges tapered (0.5, or 0.05 when there are orts); (3) the
/// orts are filtered the same way and then projected out of every series, so the
/// filter cannot reintroduce them.
///
/// Returns how many dimensions were removed (see [`BandpassReport`]), which includes
/// the detrend (a line counts 1, a quadratic 2, a mean 0, as in AFNI), the removed
/// frequency bins, and the orts, all scaled by `len / nfft` when the FFT was
/// zero-padded.
///
/// Errors: any series of the wrong length, fewer than 9 samples when a band or
/// detrend is requested, as many orts as samples, or a bad band (see [`band_plan`]).
/// With no band, no detrend and no orts, nothing happens and the report says so.
pub fn bandpass_vectors(
    vectors: &mut [Vec<f64>],
    band: Option<&BandSpec>,
    detrend: Detrend,
    orts: &[Vec<f64>],
) -> Result<BandpassReport> {
    let len = vectors.first().map_or(0, Vec::len);
    if band.is_none() && detrend == Detrend::None && orts.is_empty() {
        return Ok(BandpassReport {
            removed_dof: 0,
            nfft: 0,
        });
    }
    if len < MIN_BANDPASS_SAMPLES {
        return Err(Error::InvalidParameter {
            name: "series length".into(),
            reason: format!("filtering needs at least {MIN_BANDPASS_SAMPLES} samples, got {len}"),
        });
    }
    if let Some(bad) = vectors.iter().position(|v| v.len() != len) {
        return Err(Error::LengthMismatch {
            what: format!("series {bad}"),
            expected: len,
            found: vectors[bad].len(),
        });
    }
    if orts.len() >= len {
        return Err(Error::InvalidParameter {
            name: "orts".into(),
            reason: format!("{} regressors for {len} samples", orts.len()),
        });
    }
    let plan = match band {
        Some(b) => Some(band_plan(len, b)?),
        None => None,
    };
    let nfft = plan.map_or(len + len % 2, |p| p.nfft);
    let mut ndof = 0_usize;

    // 1. Detrend. (AFNI counts nothing for a plain mean here; the FFT step below
    //    counts the zero-frequency bin.)
    match detrend {
        Detrend::Linear => ndof += 1,
        Detrend::Quadratic => ndof += 2,
        _ => {}
    }
    vectors.iter_mut().for_each(|v| detrend.apply(v));

    // 2. Bandpass by FFT, two real series per complex FFT like AFNI.
    let mut filtered = false;
    if let Some(plan) = plan.filter(|p| p.filters()) {
        filtered = true;
        let nby2 = plan.nfft / 2;
        ndof += 2; // the zero and Nyquist bins
        if plan.jbot >= 1 {
            ndof += 2 * plan.jbot - 1;
        }
        ndof += 2 * (nby2 - plan.jtop) - 1;
        let taper = if orts.is_empty() { 0.5 } else { 0.05 };
        filter_series(vectors, &plan, taper);
    }

    // 3. Orts: filter a copy of each the same way, then project them out.
    if !orts.is_empty() {
        let mut q: Vec<Vec<f64>> = orts.to_vec();
        // The recursion in AFNI has no orts of its own, so the taper is 0.5 here.
        bandpass_vectors(&mut q, band, detrend, &[])?;
        let projector = OrtProjector::new(&q, len)?;
        vectors.iter_mut().for_each(|v| projector.apply(v));
        ndof += orts.len();
    }

    if nfft > len {
        ndof = rint_f64(len as f64 / nfft as f64 * ndof as f64) as usize;
    }
    Ok(BandpassReport {
        removed_dof: ndof,
        nfft: if filtered { nfft } else { 0 },
    })
}

/// FFT-filter every series in place (pairs share one complex FFT).
fn filter_series(vectors: &mut [Vec<f64>], plan: &BandPlan, taper: f64) {
    let (nfft, nby2) = (plan.nfft, plan.nfft / 2);
    let len = vectors[0].len();
    let mut buffer = vec![Complex::ZERO; nfft];
    let mut i = 0;
    while i < vectors.len() {
        let paired = i + 1 < vectors.len();
        // Load one or two series: x in the real part, y in the imaginary part.
        buffer.iter_mut().for_each(|c| *c = Complex::ZERO);
        for k in 0..len {
            buffer[k].re = vectors[i][k];
            if paired {
                buffer[k].im = vectors[i + 1][k];
            }
        }
        fft(&mut buffer, false);
        // Zero the mean and Nyquist bins; taper the band edges; zero the rest.
        buffer[0] = Complex::ZERO;
        buffer[nby2] = Complex::ZERO;
        if plan.jbot >= 1 {
            buffer[plan.jbot] = buffer[plan.jbot].scale(taper);
            buffer[nfft - plan.jbot] = buffer[nfft - plan.jbot].scale(taper);
            for j in 1..plan.jbot {
                buffer[j] = Complex::ZERO;
                buffer[nfft - j] = Complex::ZERO;
            }
        }
        buffer[plan.jtop] = buffer[plan.jtop].scale(taper);
        buffer[nfft - plan.jtop] = buffer[nfft - plan.jtop].scale(taper);
        for j in (plan.jtop + 1)..nby2 {
            buffer[j] = Complex::ZERO;
            buffer[nfft - j] = Complex::ZERO;
        }
        fft(&mut buffer, true);
        let inv = 1.0 / nfft as f64;
        for k in 0..len {
            vectors[i][k] = buffer[k].re * inv;
            if paired {
                vectors[i + 1][k] = buffer[k].im * inv;
            }
        }
        i += 2;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ramp(n: usize) -> Vec<f64> {
        (0..n).map(|i| i as f64).collect()
    }

    #[test]
    fn legendre_polynomials_and_basis() {
        assert_eq!(legendre(0.5, 0), 1.0);
        assert_eq!(legendre(0.5, 1), 0.5);
        assert!((legendre(0.5, 2) - (3.0 * 0.25 - 1.0) / 2.0).abs() < 1e-15);
        assert!((legendre(0.3, 3) - (5.0 * 0.027 - 3.0 * 0.3) / 2.0).abs() < 1e-14);
        // P_k(1) = 1 for every order.
        for m in 0..30 {
            assert!((legendre(1.0, m) - 1.0).abs() < 1e-12);
        }
        let basis = legendre_basis(3, 5).unwrap();
        assert_eq!(basis[0], vec![1.0; 5]);
        assert_eq!(basis[1], vec![-1.0, -0.5, 0.0, 0.5, 1.0]);
        assert!(legendre_basis(0, 5).is_err() && legendre_basis(5, 5).is_err());
    }

    #[test]
    fn detrending_removes_exactly_what_it_should() {
        let n = 20;
        // A line plus a constant disappears under a linear detrend...
        let mut v: Vec<f64> = ramp(n).iter().map(|t| 3.0 + 0.5 * t).collect();
        Detrend::Linear.apply(&mut v);
        assert!(v.iter().all(|x| x.abs() < 1e-12));
        // ...a quadratic under the quadratic one, but not under the linear one.
        let quad: Vec<f64> = ramp(n).iter().map(|t| 1.0 + t + 0.1 * t * t).collect();
        let mut q = quad.clone();
        Detrend::Quadratic.apply(&mut q);
        assert!(q.iter().all(|x| x.abs() < 1e-10));
        let mut l = quad.clone();
        Detrend::Linear.apply(&mut l);
        assert!(l.iter().any(|x| x.abs() > 0.1));
        // The mean detrend changes only the level.
        let mut m = vec![1.0, 2.0, 3.0, 6.0];
        Detrend::Mean.apply(&mut m);
        assert_eq!(m, vec![-2.0, -1.0, 0.0, 3.0]);
        // AFNI leaves short series alone: line < 3 samples, quadratic < 4, mean < 2.
        let mut two = vec![1.0, 5.0];
        Detrend::Linear.apply(&mut two);
        assert_eq!(two, vec![1.0, 5.0]);
        let mut three = vec![1.0, 5.0, 2.0];
        Detrend::Quadratic.apply(&mut three);
        assert_eq!(three, vec![1.0, 5.0, 2.0]);
        let mut one = vec![4.0];
        Detrend::Mean.apply(&mut one);
        assert_eq!(one, vec![4.0]);
        let mut none = vec![1.0, 2.0, 3.0];
        Detrend::None.apply(&mut none);
        assert_eq!(none, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn ort_projection_is_least_squares() {
        let n = 12;
        let x: Vec<f64> = ramp(n);
        let ones = vec![1.0; n];
        // y = 2 x + 5 + wiggle; removing {1, x} leaves only the wiggle's residual.
        let wiggle: Vec<f64> = (0..n)
            .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
            .collect();
        let mut y: Vec<f64> = (0..n).map(|i| 2.0 * x[i] + 5.0 + wiggle[i]).collect();
        let p = OrtProjector::new(&[ones.clone(), x.clone()], n).unwrap();
        assert_eq!(p.rank(), 2);
        p.apply(&mut y);
        // What is left is orthogonal to both regressors.
        let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(p, q)| p * q).sum::<f64>();
        assert!(dot(&y, &ones).abs() < 1e-10 && dot(&y, &x).abs() < 1e-10);
        // A duplicate or zero regressor adds no dimension.
        let p = OrtProjector::new(&[ones.clone(), ones.clone(), vec![0.0; n]], n).unwrap();
        assert_eq!(p.rank(), 1);
        assert!(OrtProjector::new(&[vec![1.0; 3]], n).is_err());
        assert!(OrtProjector::new(&[vec![f64::NAN; n]], n).is_err());
    }

    #[test]
    fn normalization_matches_afni_including_the_tiny_vector_case() {
        let mut v = vec![3.0, 4.0];
        assert_eq!(normalize_l2(&mut v), Some(0.2));
        assert!((v[0] - 0.6).abs() < 1e-15 && (v[1] - 0.8).abs() < 1e-15);
        // Squared length <= 1e-20 is left alone, like THD_normalize.
        let mut tiny = vec![1e-11, 0.0];
        assert_eq!(normalize_l2(&mut tiny), None);
        assert_eq!(tiny, vec![1e-11, 0.0]);
        let mut empty: Vec<f64> = vec![];
        assert_eq!(normalize_l2(&mut empty), None);
    }

    #[test]
    fn fft_matches_a_direct_dft_for_every_kind_of_length() {
        for n in [1, 2, 3, 5, 8, 12, 30, 53, 64, 100, 138] {
            let signal: Vec<Complex> = (0..n)
                .map(|k| Complex {
                    re: (k as f64 * 0.7).sin() + 0.1 * k as f64,
                    im: (k as f64 * 1.3).cos(),
                })
                .collect();
            let mut fast = signal.clone();
            fft(&mut fast, false);
            for (k, got) in fast.iter().enumerate() {
                let mut want = Complex::ZERO;
                for (t, x) in signal.iter().enumerate() {
                    let a = -2.0 * std::f64::consts::PI * (k * t) as f64 / n as f64;
                    let w = Complex {
                        re: a.cos(),
                        im: a.sin(),
                    };
                    let p = x.mul(w);
                    want.re += p.re;
                    want.im += p.im;
                }
                assert!(
                    (got.re - want.re).abs() < 1e-9 && (got.im - want.im).abs() < 1e-9,
                    "n={n} k={k}"
                );
            }
            // And the inverse undoes it (after scaling by 1/n).
            fft(&mut fast, true);
            for (a, b) in fast.iter().zip(&signal) {
                assert!((a.re / n as f64 - b.re).abs() < 1e-10);
                assert!((a.im / n as f64 - b.im).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn public_power_spectrum_makes_padding_layout_and_scaling_explicit() {
        let full_raw = PowerSpectrumOptions::new(
            Some(4),
            SpectrumLayout::Full,
            SpectrumNormalization::Raw,
            NonFinitePolicy::Reject,
        );
        let constant = power_spectrum(&[1.0; 4], full_raw).unwrap();
        assert_eq!(constant.len(), 4);
        assert!((constant[0] - 16.0).abs() < 1.0e-12);
        assert!(constant[1..].iter().all(|power| power.abs() < 1.0e-12));

        // Zero padding a unit impulse gives unit power in every full-spectrum
        // bin, which also fixes the forward-transform normalization convention.
        assert_eq!(power_spectrum(&[1.0, 0.0], full_raw).unwrap(), vec![1.0; 4]);

        let one_sided = PowerSpectrumOptions::new(
            Some(4),
            SpectrumLayout::OneSided,
            SpectrumNormalization::Raw,
            NonFinitePolicy::Reject,
        );
        assert_eq!(power_spectrum(&[1.0, 0.0], one_sided).unwrap().len(), 3);

        let parseval = PowerSpectrumOptions::new(
            None,
            SpectrumLayout::Full,
            SpectrumNormalization::DivideByFftLength,
            NonFinitePolicy::Reject,
        );
        let signal = [1.0, -2.0, 3.0, -4.0, 5.0];
        let spectral_energy: f64 = power_spectrum(&signal, parseval).unwrap().iter().sum();
        let signal_energy: f64 = signal.iter().map(|value| value * value).sum();
        assert!((spectral_energy - signal_energy).abs() < 1.0e-10);

        assert!(power_spectrum(
            &signal,
            PowerSpectrumOptions {
                fft_len: Some(4),
                ..parseval
            }
        )
        .is_err());
        assert!(power_spectrum(&[], parseval).is_err());
    }

    #[test]
    fn public_power_spectrum_non_finite_policy_is_order_preserving() {
        let options = |non_finite| {
            PowerSpectrumOptions::new(
                None,
                SpectrumLayout::Full,
                SpectrumNormalization::Raw,
                non_finite,
            )
        };
        let values = [1.0, f64::NAN, 2.0];
        assert!(matches!(
            power_spectrum(&values, options(NonFinitePolicy::Reject)),
            Err(Error::NonFinite { ref what, value })
                if what == "power spectrum sample 1" && value.is_nan()
        ));
        assert!(matches!(
            power_spectrum(&values, options(NonFinitePolicy::Skip)),
            Err(Error::InvalidParameter { ref name, .. })
                if name == "power spectrum non-finite policy"
        ));
        let propagated = power_spectrum(&values, options(NonFinitePolicy::Propagate)).unwrap();
        assert!(propagated.iter().all(|power| power.is_nan()));
    }

    #[test]
    fn band_plans_follow_afni() {
        // 100 samples at 2 s: nfft 100, df 0.005, band 0.01..0.1 is bins 2..20
        // (the same numbers 1dBandpass prints).
        let plan = band_plan(
            100,
            &BandSpec {
                dt: 2.0,
                fbot: 0.01,
                ftop: 0.1,
            },
        )
        .unwrap();
        assert_eq!(
            plan,
            BandPlan {
                nfft: 100,
                jbot: 2,
                jtop: 20
            }
        );
        assert!(plan.filters());
        // Odd lengths round up to even: 127 -> 128, bins 3..26.
        let plan = band_plan(
            127,
            &BandSpec {
                dt: 2.0,
                fbot: 0.01,
                ftop: 0.1,
            },
        )
        .unwrap();
        assert_eq!(
            plan,
            BandPlan {
                nfft: 128,
                jbot: 3,
                jtop: 26
            }
        );
        // A top above Nyquist is clipped, a bottom of 0 means low-pass: nothing cut
        // at the bottom and nothing at the top means no filtering at all.
        let all = band_plan(
            100,
            &BandSpec {
                dt: 1.0,
                fbot: 0.0,
                ftop: MAX_FTOP_HZ,
            },
        )
        .unwrap();
        assert!(!all.filters());
        // Bad input is an error.
        let ok = BandSpec {
            dt: 1.0,
            fbot: 0.01,
            ftop: 0.1,
        };
        assert!(band_plan(8, &ok).is_err());
        assert!(band_plan(100, &BandSpec { dt: 0.0, ..ok }).is_err());
        assert!(band_plan(100, &BandSpec { fbot: -1.0, ..ok }).is_err());
        assert!(band_plan(100, &BandSpec { ftop: 0.01, ..ok }).is_err());
        assert!(band_plan(
            100,
            &BandSpec {
                ftop: f64::NAN,
                ..ok
            }
        )
        .is_err());
    }

    #[test]
    fn bandpass_keeps_in_band_tones_and_removes_the_rest() {
        let n = 200;
        let dt = 1.0;
        // 0.05 Hz is inside 0.02..0.1; 0.3 Hz and a slow drift are outside.
        let tone = |hz: f64| -> Vec<f64> {
            (0..n)
                .map(|t| (2.0 * std::f64::consts::PI * hz * t as f64 * dt).sin())
                .collect()
        };
        let mixed: Vec<f64> = (0..n)
            .map(|t| tone(0.05)[t] + tone(0.3)[t] + 0.01 * t as f64 + 4.0)
            .collect();
        let mut v = vec![mixed];
        let band = BandSpec {
            dt,
            fbot: 0.02,
            ftop: 0.1,
        };
        let report = bandpass_vectors(&mut v, Some(&band), Detrend::Quadratic, &[]).unwrap();
        // Close to the in-band tone (edge effects of a 200-sample window remain).
        let want = tone(0.05);
        let err: f64 = (20..180)
            .map(|t| (v[0][t] - want[t]).abs())
            .fold(0.0, f64::max);
        assert!(err < 0.2, "max error {err}");
        // Dimensions removed: quadratic 2 + (0, Nyquist) 2 + low (2*4-1) + high (2*(100-20)-1).
        let plan = band_plan(n, &band).unwrap();
        let expected = 2 + 2 + (2 * plan.jbot - 1) + (2 * (n / 2 - plan.jtop) - 1);
        assert_eq!(report.removed_dof, expected);
        assert_eq!(report.nfft, n);
    }

    #[test]
    fn removed_dof_scales_when_the_fft_is_padded() {
        // 127 samples -> nfft 128: AFNI scales the count by 127/128 and rounds.
        let band = BandSpec {
            dt: 2.0,
            fbot: 0.01,
            ftop: 0.1,
        };
        let mut v = vec![(0..127)
            .map(|i| (i as f64 * 0.37).sin())
            .collect::<Vec<_>>()];
        let r = bandpass_vectors(&mut v, Some(&band), Detrend::Linear, &[]).unwrap();
        // 1 + 2 + (2*3-1) + (2*(64-26)-1) = 83; 83 * 127/128 = 82.35 -> 82.
        assert_eq!((r.removed_dof, r.nfft), (82, 128));
    }

    #[test]
    fn orts_are_filtered_then_projected_and_counted() {
        let n = 100;
        let band = BandSpec {
            dt: 1.0,
            fbot: 0.02,
            ftop: 0.12,
        };
        let ort: Vec<f64> = (0..n)
            .map(|t| (2.0 * std::f64::consts::PI * 0.06 * t as f64).sin())
            .collect();
        let signal: Vec<f64> = (0..n)
            .map(|t| 2.0 * ort[t] + (2.0 * std::f64::consts::PI * 0.09 * t as f64).cos())
            .collect();
        let mut with = vec![signal.clone()];
        let r = bandpass_vectors(
            &mut with,
            Some(&band),
            Detrend::Quadratic,
            std::slice::from_ref(&ort),
        )
        .unwrap();
        let mut without = vec![signal.clone()];
        let r0 = bandpass_vectors(&mut without, Some(&band), Detrend::Quadratic, &[]).unwrap();
        assert_eq!(r.removed_dof, r0.removed_dof + 1);
        // What remains is orthogonal to the FILTERED ort.
        let mut q = vec![ort];
        bandpass_vectors(&mut q, Some(&band), Detrend::Quadratic, &[]).unwrap();
        let dot: f64 = with[0].iter().zip(&q[0]).map(|(a, b)| a * b).sum();
        assert!(dot.abs() < 1e-8, "dot {dot}");
        // The ort's contribution really was removed from the series.
        let energy = |v: &[f64]| v.iter().map(|x| x * x).sum::<f64>();
        assert!(energy(&with[0]) < energy(&without[0]));
    }

    #[test]
    fn degenerate_calls() {
        // Nothing requested: nothing done.
        let mut v = vec![vec![1.0, 2.0]];
        let r = bandpass_vectors(&mut v, None, Detrend::None, &[]).unwrap();
        assert_eq!(r.removed_dof, 0);
        assert_eq!(v, vec![vec![1.0, 2.0]]);
        // Too short, ragged, or too many orts are errors.
        let mut short = vec![vec![0.0; 5]];
        assert!(bandpass_vectors(&mut short, None, Detrend::Linear, &[]).is_err());
        let mut ragged = vec![vec![0.0; 12], vec![0.0; 11]];
        assert!(bandpass_vectors(&mut ragged, None, Detrend::Linear, &[]).is_err());
        let mut ok = vec![vec![0.0; 10]];
        let orts = vec![vec![1.0; 10]; 10];
        assert!(bandpass_vectors(&mut ok, None, Detrend::None, &orts).is_err());
        // The remaining dimension: 2 * (kept bins), capped by the length.
        let band = BandSpec {
            dt: 2.0,
            fbot: 0.01,
            ftop: 0.1,
        };
        assert_eq!(
            bandpass_remaining_dimension(100, &band).unwrap(),
            2 * (20 - 2 + 1)
        );
    }
}
