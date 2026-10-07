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
// The single error type (`Error`) and result alias (`Result<T>`) used by every
// fallible function in `afni-core`.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// * `numeric.rs` is the first user: its checked conversions return these errors.
// * Later phases (statistics, datasets, colors, ...) will add variants here
//   rather than defining their own error types, so that callers only ever need
//   to handle one `Error`.
// * `afni-io` has its own, separate `Error` for file problems (bad bytes, missing
//   files). Core errors are about *meaning* (a non-finite statistic, an index
//   outside a domain), never about files. When an `afni-io` adapter fails, it
//   will wrap or convert into whichever of the two errors fits.
// ---------------------------------------------------------------------------

//! Error and result types shared across the crate.

/// Convenience alias for `Result<T, afni_core::Error>`.
///
/// Writing `Result<f64>` instead of `std::result::Result<f64, Error>` keeps
/// signatures short. Because it shadows the standard `Result` only inside
/// modules that import it, there is no ambiguity elsewhere.
pub type Result<T> = std::result::Result<T, Error>;

/// Everything that can go wrong inside `afni-core`.
///
/// # Design notes for readers new to Rust
///
/// * `#[derive(thiserror::Error)]` writes the boilerplate that makes this a
///   proper `std::error::Error`; the `#[error("...")]` text on each variant is
///   what `Display` (and therefore `println!("{err}")`) prints.
/// * `#[non_exhaustive]` promises callers that we may add variants later.
///   Downstream `match` statements must include a `_ =>` arm, so adding a
///   variant in a later phase is not a breaking change.
/// * Every variant carries enough context (names, values, lengths) to explain
///   the failure without a debugger. The roadmap requires that invalid metadata
///   produce an explicit error, never a plausible-looking guessed value.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// A value that must be finite was NaN or infinite.
    ///
    /// Statistical code never silently replaces such values (see
    /// [`crate::numeric::NonFinitePolicy`]).
    #[error("{what} must be finite, got {value}")]
    NonFinite {
        /// Human-readable name of the offending quantity, e.g. `"t statistic"`.
        what: String,
        /// The value that was rejected.
        value: f64,
    },

    /// An index was outside `0..len`.
    ///
    /// `index` is an `i64` so that negative indices from files (which are
    /// signed in several AFNI formats) can be reported faithfully instead of
    /// wrapping around to a huge unsigned number.
    #[error("index {index} is out of range for a domain of {len} samples")]
    IndexOutOfRange {
        /// The index as it was supplied.
        index: i64,
        /// Number of valid samples; valid indices are `0..len`.
        len: usize,
    },

    /// A `usize` did not fit in the narrower integer type a buffer requires.
    #[error("value {value} does not fit in a 32-bit index")]
    IndexOverflow {
        /// The value that was too large.
        value: usize,
    },

    /// A dataset column index was outside `0..len`.
    #[error("column index {index} is out of range for {len} columns")]
    ColumnIndexOutOfRange {
        /// The requested zero-based column index.
        index: usize,
        /// Number of columns in the dataset.
        len: usize,
    },

    /// A parameter was finite but outside the range its operation accepts,
    /// for example a negative degrees-of-freedom count.
    #[error("invalid parameter {name}: {reason}")]
    InvalidParameter {
        /// Parameter name, e.g. `"dof"`.
        name: String,
        /// Why it was rejected.
        reason: String,
    },

    /// The same sample index appears twice in an indexed mapping.
    ///
    /// Datasets store one row per sample, so a repeated index would make
    /// "the value at node N" ambiguous. It is rejected rather than resolved by
    /// guessing which row wins.
    #[error("sample index {index} appears more than once")]
    DuplicateIndex {
        /// The repeated index.
        index: u32,
    },

    /// Two things that must have the same length did not, such as a column
    /// and the number of rows the dataset's mapping expects.
    #[error("{what}: expected {expected}, found {found}")]
    LengthMismatch {
        /// What was being compared, e.g. `"column 'beta' rows"`.
        what: String,
        /// The length required.
        expected: usize,
        /// The length supplied.
        found: usize,
    },

    /// Two objects refer to different spatial sample domains.
    ///
    /// This is distinct from a length mismatch: a two-voxel volume and a
    /// two-node surface have the same number of samples but cannot safely share
    /// a mask or per-sample result.
    #[error("domain mismatch: {what}")]
    DomainMismatch {
        /// What incompatible objects were being combined.
        what: String,
    },

    /// Something that must contain data was empty (a dataset with no columns,
    /// a column with no rows, an empty label name).
    #[error("{0} is empty")]
    Empty(String),

    /// No answer exists or one could not be found: a probability with no
    /// matching statistic, or a numerical routine that failed to converge.
    #[error("no solution: {0}")]
    NoSolution(String),

    /// The request is well-formed but this crate does not implement it yet.
    ///
    /// Returned instead of an approximation, so an unsupported statistic can
    /// never masquerade as a real p-value.
    #[error("unsupported: {0}")]
    Unsupported(String),
}
