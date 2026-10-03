# afni-core architecture

*Public domain: United States Government work (17 U.S.C. § 105); see `LICENSE`.*

This document records the Phase 0 decisions. Code comments and the roadmap
(`afni-core_ROADMAP.md`) refer back to it.

## 1. Separate projects, one-way dependency

`afni-core`, `afni-io`, and the viewers (`sumaru`, `afniru`) are **separate
Cargo projects**, each with its own `Cargo.toml`, `target/`, and (for
`afni-io`) git history. There is deliberately no shared Cargo workspace.

```text
   sumaru ───────▶ afni-core ◀─────── afniru
     │                 ▲                 │
     └───────────▶ afni-io ◀─────────────┘
                       │
                       └──▶ afni-core   (adapters only, from Phase 1)
```

* `afni-io → afni-core` is allowed: `afni-io` will gain adapters that turn
  decoded files into core types, via `afni-core = { path = "../afni-core" }`.
* `afni-core → afni-io` is **forbidden** (it would create a cycle and make the
  meaning layer depend on file formats).

**Enforcement:** `tests/dependency_direction.rs` runs `cargo tree` on the real
resolved dependency graph and fails if `afni-io`, GUI/GPU/windowing/socket
crates, or `flate2` appear. It also scans `src/` for `std::fs`, `std::path`,
`std::process`, `std::net`, `use afni_io`, and `unsafe`.

## 2. What belongs where

| Concern | Crate |
|---|---|
| Bytes, syntax, compression, raw attributes, round-trip fidelity | `afni-io` |
| Domains, typed datasets, p-values, FDR, colors, thresholds, clusters, ROIs | `afni-core` |
| Windows, GPU, cameras, picking, sockets, preferences, CLI strings | viewers/tools |

Boundary test: a public `afni-core` function must work on caller-supplied
slices and small value types, knowing no path, file format, GUI, GPU, socket,
or AFNI executable.

## 3. Type-ownership decision (roadmap Phase 0, item 3)

**Decision: move, do not mirror — and move each type only in the phase that
first needs it.**

* Semantic types (`StatKind`, `StatSpec`, `ThresholdCurve`, label semantics,
  column ranges) are owned by `afni-core` once their phase lands (2, 3, 4, 1).
* At that time `afni-io` depends on `afni-core` and **re-exports** the type
  from its old path, so `afni_io::stat::StatSpec` keeps compiling.
* Parsing a header string such as `COLMS_STATSYM` into a core type stays in
  `afni-io` as an adapter function; core never sees NIML/HEAD syntax.
* Nothing moved in Phase 0. In Phase 1, `StatKind` and `StatSpec` moved (re-exported
  by `afni-io`); the raw `ThresholdCurve` and `LabelTable` stay in `afni-io` as
  unvalidated file types, with validated core equivalents and adapters
  (`afni_io::adapt`). There is no duplicate *semantic* model. If a
  phase finds a type cannot move cleanly, record it in the roadmap discovery
  log and fall back to explicit `From`/`TryFrom` adapters with a dated
  deprecation plan.

Rationale: avoiding duplicates is easier than reconciling them later, and
moving types alongside the code that gives them meaning keeps each change
reviewable and testable against AFNI.

Raw unknown metadata that must round-trip lives in an `afni-io` wrapper
`{ core_dataset, raw_extra_metadata }`, never in `afni-core`.

## 4. Numeric and error conventions

Implemented in `src/numeric.rs` and `src/error.rs`:

* `f64` for statistical math; `f32` for stored display/mesh buffers
  (`widen_to_f64`, `narrow_to_f32`).
* Checked index conversion (`checked_index`, `usize_to_u32`); no bare `as`.
* NaN/Inf handled by an explicit `NonFinitePolicy` (`Reject`, `Skip`,
  `Propagate`) or `ensure_finite`; never silently replaced.
* One crate-wide `Error`; invalid metadata is an error or
  `Error::Unsupported`, never a guessed result.
* `unsafe_code = "forbid"`, Rust 1.74 MSRV, `#![warn(missing_docs)]`.

## 5. Conformance harness

* Fixtures are committed in `tests/data/conformance/` with `# afni_version:`
  provenance, so ordinary `cargo test` never needs AFNI.
* `tests/data/regenerate_conformance.sh` regenerates them from live AFNI.
* `AFNI_CORE_LIVE=1 cargo test` replays every case against the installed AFNI
  and checks the committed values still match.
* `tests/common/mod.rs` provides the single comparison path for later phases:
  `Fixture`, `Tolerance`, `assert_close`, `assert_monotonic`,
  `assert_round_trip`.

Current fixture: `cdf.ref`, generated with AFNI_26.2.08.

## 6. Datasets and adapters (Phase 1)

`afni_core::dataset::Dataset` = `Domain` (surface nodes or volume voxels) +
`SampleMap` (dense, or indexed/sparse) + typed `DataColumn`s (values, stat,
FDR/MDF curves, label table, recorded range). Policies: out-of-domain and
duplicate indices are errors; NaN/Inf are kept as "no data"; the node count of a
sparse surface dataset is supplied by the caller, never guessed.

`afni_io::adapt` converts `NimlDataset` (with a round-trip `NimlEnvelope`),
`Gifti`, `Volume`, and `OneD` into core datasets. Raw attributes the core model
does not understand live in `NimlExtras` on the I/O side.

## 7. Statistics (Phase 2)

`afni_core::stats` turns a `StatSpec` plus a statistic into a `Probability`
(`p` and `ln p`) for an explicit `Tail` (`Lower`, `Upper`, `TwoSided`), and back
with `critical_value`. Built on `afni_core::special` (log-space incomplete
beta/gamma and normal tails; no external math dependency).

Decisions: the tail is never inferred; two-sided is defined only for
distributions with a centre of symmetry; failures are errors, never AFNI's
sentinel values; discrete distributions use the true step function with
documented integer quantiles; gamma's second parameter is a rate; noncentral t,
F and chi-square are computed exactly (never approximated by the central forms). Correlation metadata keeps AFNI's
`Correl(samples, nfit, nort)` form; the NIfTI one-parameter form is converted
to the equivalent `Correl(dof + 1, 1, 0)` at ingest (`StatSpec::from_intent`,
`IntentOrigin`). `nfit > 1` is a multiple correlation (`R >= 0`, one-sided).
Conformance: `tests/stats_conformance.rs` against AFNI's `nifticdf` (412 cases at
full precision), `cdf`, `ccalc`, and `p2dsetstat`.

## 8. FDR and q-values (Phase 3)

`afni_core::fdr` uses stored FDR/MDF curves (`q_value_for_threshold`,
`threshold_for_q`, `missed_detection_fraction`, `minimum_q`) and builds them from
data (`fdrize`, `fdr_curves`, a port of `mri_fdrize`/`mri_fdr_curve`). The
interpolation is AFNI's four-point clamped cubic and inverse
(`ThresholdCurve::interpolate`, `inverse_interpolate`). p-values
(`stats::Probability`) and q-values (`fdr::QValue`) are distinct types. A plain
Benjamini-Hochberg/Yekutieli helper (`benjamini_hochberg`) is kept separate from
AFNI's algorithm, which floors p, ignores p >= 0.9999, estimates the number of true
positives, and works in single precision. Conformance: `afni-io/tests/fdr_conformance.rs`
against `3drefit -addFDR` curves, `3dFDR`, and `fdrval`. The roadmap discovery log
lists every AFNI quirk that was copied or guarded.

## 9. Colors (Phase 4)

`afni_core::color` holds `Rgba` (f32, straight alpha), `ColorStop`,
`ContinuousColorMap` (explicit `Interpolation` and `InterpolationSpace`, hard edges
at duplicate stops, GPU-ready `lookup_table`), and `ColorMap` (continuous or
labels). `afni_core::afni_colors` reproduces AFNI's nine 256-entry scales exactly
(`AfniColorScale`), including `DC_spectrum_AJJ`/`ZSS`; golden tables are generated
from AFNI's own C functions. `afni_core::labels` adds label colors: a
`LabelColorPolicy` (unlabeled keys, uncolored entries, key 0), the stable fallback
palette, and `LabelColorMap`. File label tables convert through
`afni_io::adapt::label_table_to_core` with keys, order and colors preserved.

## 10. Thresholds, overlays and compositing (Phase 5)

`threshold` defines `Threshold` (exact boundary behavior per mode; `AbsoluteAbove`
and SUMA's `Outside` are distinct because SUMA's own modes disagree at the
boundary), `FadeModel` (AFNI volume fade with byte rounding and the 222/255
ceiling; SUMA's surface falloff; a generalized boundary fade), and
`transfer_threshold` (match p-values across statistics, preserving the tail).
`overlay` turns data plus an immutable `OverlaySpec` into colors, a pass mask and
diagnostics (`evaluate_rows`, `evaluate_dataset`); it has no cache or viewer
state. Colors come from a continuous map, a pane table with SUMA/AFNI lookup rules
(`ColorTable`, `PaneRule`), or label keys; `GpuOverlayParams` flattens a spec for a
shader. `composite` layers straight-alpha colors with the Porter-Duff over
operator. SUMA parity is tested against `ScaleToMap` (54 cases); the roadmap
discovery log lists what that oracle cannot reach.

## 11. Surface topology, geometry and clustering (Phase 6)

`topology` validates a triangle list and derives neighbors, edges, face neighbors,
components, ring layers and a diagnostic `TopologyReport` (boundary, non-manifold,
winding, bow-ties, Euler characteristic); only an out-of-range index or empty node
set is an error. `mesh` adds coordinates: triangle and node normals (SUMA's
convention: the normalized sum of unit triangle normals), areas (a third of the
adjacent triangles per node), bounds, edge lengths, signed volume, and a reusable
shortest-path `NeighborhoodSearcher`. `cluster` groups active nodes
(`Connectivity::{EdgeRings, GraphDistance}`, `Tails::{Merged, Separate}`, size limits,
three sort orders) and returns labels, per-cluster statistics and a `survivor_mask`
that plugs into `overlay`. `ClusterNeighborhoods` precomputes searches for repeated
clustering. Conformance is against `SurfaceMetrics`, `SurfMeasures` and `SurfClust`
(`tests/mesh_conformance.rs`, `tests/cluster_conformance.rs`); the one known
disagreement (SUMA's millimetre radius search) and its direction are pinned by a test.
