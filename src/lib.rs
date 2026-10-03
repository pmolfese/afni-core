// PUBLIC DOMAIN NOTICE
//
// afni-core was written by employees of the United States Government (National
// Institutes of Health) as part of their official duties. Under Title 17,
// Section 105 of the United States Code it is a "United States Government
// Work" and is not subject to copyright protection in the United States. It is
// therefore in the public domain in the United States and may be used, copied,
// modified and redistributed without restriction. Outside the United States,
// rights are waived under the Creative Commons CC0 1.0 Universal Public Domain
// Dedication. See the LICENSE file at the crate root.
//
// The software is provided "as is", without warranty of any kind. See LICENSE.
//
// ---------------------------------------------------------------------------
// WHAT THIS FILE IS
//
// The crate root. It holds the crate-level documentation (shown on docs.rs and
// by `cargo doc`), declares the modules, and sets crate-wide lint policy.
//
// HOW IT RELATES TO THE REST OF THE CRATE
//
// Each `pub mod` below is one area of functionality. Right now only the
// foundation modules exist (Phase 0 of afni-core_ROADMAP.md); the table in the
// documentation lists what later phases will add.
// ---------------------------------------------------------------------------

//! File-neutral data models and algorithms for AFNI/SUMA tools.
//!
//! `afni-core` decides what decoded neuroimaging data *means*: which
//! statistical test a column represents and what p-value a value implies, how
//! numbers become colors, which samples pass a threshold, and how clusters and
//! regions of interest behave. It is shared by the `sumaru` surface viewer, the
//! `afniru` slice viewer, command-line tools, and tests.
//!
//! # What this crate does **not** do
//!
//! This is enforced by a test (`tests/dependency_direction.rs`), not just
//! convention. `afni-core` contains no:
//!
//! * file parsing or writing (that is the job of `afni-io`),
//! * GUI toolkits, GPU APIs, or shader code,
//! * sockets or AFNI/SUMA "talk" protocol,
//! * process launching or dependence on AFNI executables at runtime.
//!
//! Every public operation works on caller-supplied slices and small value
//! types, so the same code can run in a viewer, a CLI, or a test, and a GPU
//! path can share its semantics.
//!
//! # Where this crate sits
//!
//! ```text
//!   sumaru ───────▶ afni-core ◀─────── afniru
//!     │                 ▲                 │
//!     └───────────▶ afni-io ◀─────────────┘
//!                       │
//!                       └──▶ afni-core   (adapters only)
//! ```
//!
//! The arrow `afni-io → afni-core` is allowed (so `afni-io` can offer
//! "convert this file's data into core types"). The reverse is **forbidden**:
//! it would make the meaning layer depend on file formats and create a cycle.
//! See `docs/ARCHITECTURE.md` for the reasoning and the type-ownership plan.
//!
//! # Modules
//!
//! | Module | Purpose | Roadmap phase |
//! |--------|---------|---------------|
//! | [`error`] | The crate-wide [`Error`] and [`Result`] | 0 |
//! | [`numeric`] | `f32`/`f64` rules, checked indices, NaN/Inf policy | 0 |
//! | [`domain`] | Surface node sets and volume voxel grids | 1 |
//! | [`mapping`] | Dense and indexed (sparse) row-to-sample maps | 1 |
//! | [`column`](mod@column) | Typed columns and column metadata | 1 |
//! | [`dataset`] | The validated [`Dataset`](dataset::Dataset) model | 1 |
//! | [`stat`] | `StatKind` / `StatSpec` (what a statistic *is*) | 1 (math in 2) |
//! | [`special`] | log-gamma/beta, incomplete beta/gamma, normal tails (in log space) | 2 |
//! | [`stats`] | Tail-aware p-values and critical values | 2 |
//! | [`color`] | `Rgba`, color stops, continuous maps, interpolation modes | 4 |
//! | [`afni_colors`] | AFNI's built-in color scales, exactly | 4 |
//! | [`suma_colormaps`] | SUMA's standard named color maps, exactly | 4 |
//! | [`composite`] | Alpha compositing (straight alpha) of underlay and overlay planes | 5 |
//! | [`overlay`] | Overlay evaluation: colors, pass mask, diagnostics | 5 |
//! | [`cluster`] | Connected-cluster labeling on a surface (SurfClust-compatible) | 6 |
//! | [`mesh`] | Mesh geometry (normals, areas, volume) and distance searches | 6 |
//! | [`topology`] | Validated triangle-mesh connectivity, diagnostics, rings | 6 |
//! | [`threshold`] | Thresholds, transparent thresholding, matched-p transfer | 5 |
//! | [`fdr`] | FDR/MDF curves, q-values, Benjamini-Hochberg | 3 |
//! | [`curve`] | Validated FDR/MDF curve tables | 1 (lookups in 3) |
//! | [`labels`] | Label tables (keys, names, colors) | 1 (colors in 4) |
//! | *fdr* | FDR/MDF curves and multiple-comparison helpers | 3 |
//! | *overlay* | Thresholding, overlay evaluation, compositing | 5 |
//! | *surface, cluster* | Topology, geometry, surface clustering | 6 |
//! | *volume* | Voxel neighborhoods and volume clustering | 7 |
//! | *roi* | ROI model and set operations | 8 |
//! | *timeseries* | Preprocessing and seed correlation | 9 |
//!
//! # Numeric conventions in one paragraph
//!
//! Statistical math uses `f64`; stored display buffers use `f32`; index
//! conversions are checked rather than cast; and NaN/infinity are always
//! handled by an explicit [`numeric::NonFinitePolicy`]. Invalid metadata yields
//! an [`Error`], never a plausible-looking guess. Details are in [`numeric`].
//!
//! # Statistics at a glance
//!
//! ```
//! use afni_core::stat::{StatKind, StatSpec};
//! use afni_core::stats::{critical_value, p_value, Tail};
//!
//! // A t statistic with 23 degrees of freedom, as AFNI writes `Ttest(23)`.
//! let spec = StatSpec::parse("Ttest(23)").unwrap();
//! // The tail is always your choice; there is no hidden default.
//! let p = p_value(&spec, 2.5, Tail::TwoSided)?;
//! println!("p = {:.4} (ln p = {:.3})", p.p(), p.ln_p());
//! // And back: the |t| that a two-sided p of 0.01 requires.
//! let t = critical_value(&spec, 0.01, Tail::TwoSided)?;
//! assert!((t - 2.8073).abs() < 1e-4);
//! # Ok::<(), afni_core::Error>(())
//! ```
//!
//! # Quick start
//!
//! ```
//! use afni_core::numeric::{checked_index, NonFinitePolicy};
//!
//! // A sparse surface dataset says "this row belongs to node 41".
//! // Before using it as a slice index, check it against the domain size.
//! let node_count = 100;
//! let node = checked_index(41, node_count)?;
//! assert_eq!(node, 41);
//!
//! // A NaN sample is skipped rather than silently turned into a number.
//! assert_eq!(NonFinitePolicy::Skip.check("sample", f64::NAN)?, None);
//! # Ok::<(), afni_core::Error>(())
//! ```

// Missing docs are warnings so the "plenty of comments" goal is machine-checked:
// every public item must be documented.
#![warn(missing_docs)]
#![warn(missing_debug_implementations)]

pub mod afni_colors;
pub mod cluster;
pub mod color;
pub mod column;
pub mod composite;
pub mod curve;
pub mod dataset;
pub mod domain;
pub mod error;
pub mod fdr;
pub mod labels;
pub mod mapping;
pub mod mesh;
pub mod numeric;
pub mod overlay;
pub mod special;
pub mod stat;
pub mod stats;
pub mod suma_colormaps;
pub mod threshold;
pub mod topology;

// Re-export the two most-used names at the crate root so callers can write
// `afni_core::Error` instead of `afni_core::error::Error`.
pub use error::{Error, Result};

// Compile (and run) the Rust example in README.md as a doctest, so the README
// cannot drift away from the API. Only built under `cargo test`.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeDoctests;
