# afni-core roadmap

Goal: build the file-neutral data models and algorithms shared by `sumaru`, a
future 2D AFNI slice viewer, command-line tools, and tests. `afni-io` should
decode and encode files; `afni-core` should decide what statistical metadata
means, how values become colors, which samples pass a threshold, and how
clusters, ROIs, and derived datasets behave.

This roadmap follows the existing long-term split described in
`UPDATE_AFNI_CRATE.md` and `afni-io_ROADMAP.md`. It expands the originally
named `stats`, `color`, `overlay`, `cluster`, and `dataset` areas after
reviewing the corresponding `sumaru` modules and AFNI's reference C source.

Status: ✅ done · 🚧 in progress · ⬜ not started

| # | Phase | Status |
|---|-------|--------|
| 0 | Workspace, dependency direction, and conformance harness | ✅ |
| 1 | File-neutral domains and datasets | ✅ |
| 2 | Statistical interpretation and p-values | ✅ |
| 3 | FDR/MDF curves and multiple-comparison helpers | ✅ |
| 4 | Colors, continuous maps, and label tables | ✅ |
| 5 | Thresholding, overlay evaluation, and compositing | ✅ |
| 6 | Surface topology, geometry metrics, and clustering | ✅ |
| 7 | Volume neighborhoods and clustering | ⬜ |
| 8 | File-neutral ROI model and operations | ⬜ |
| 9 | Time-series preprocessing and seed correlation | ⬜ |
| 10 | Graph, tract, and other derived semantic models | ⬜ |
| 11 | Consumer migration, performance, and API stabilization | ⬜ |

Open questions awaiting a decision are numbered in [Open decisions](#open-decisions) (none open; D1–D20 are resolved as R6–R25).

---

## Design boundaries

Target dependency direction:

```text
afni-core  <-  afni-io
    ^              ^
    |              |
 sumaru        afni-view
```

- `afni-core` owns semantic types and pure algorithms. It must not depend on
  file paths, NIML/XML/header syntax, GUI toolkits, GPU APIs, sockets, or AFNI
  executables.
- `afni-io` keeps raw format representations where round-trip fidelity needs
  them, and adds adapters into `afni-core` types. It may depend on
  `afni-core`; the reverse dependency would make the semantic layer depend on
  file formats and should be avoided in the final workspace.
- `sumaru` and `afni-view` own interaction state, rendering resources,
  background jobs, command construction, and user preferences. They consume
  deterministic results from `afni-core`.
- Algorithms should operate on slices and small value types where practical,
  so a CPU viewer, a CLI, and a GPU upload path can share the same semantics.
- Preserve `unsafe_code = "forbid"`, Rust 1.74 until the workspace explicitly
  raises its MSRV, and a small dependency surface. Heavy or specialized
  algorithms such as FFT preprocessing should be feature-gated if their
  dependencies are not needed by every consumer.
- Invalid metadata should return typed errors or an explicit unsupported
  result. Do not silently invent degrees of freedom, replace non-finite data,
  or guess a statistical tail.

### What stays outside afni-core

- File parsing and writing: `afni-io`.
- AFNI/SUMA talk sockets and message framing: `sumaru` or an I/O/protocol
  crate.
- `egui`, `wgpu`, shaders, camera controls, picking, screenshots, and windows:
  viewers.
- CLI string construction such as `SurfClust ...`: consumer tooling. The
  underlying cluster parameters and results belong here.
- Persistent preferences and project/session state: consumers.

---

## Source map

These are the starting points, not instructions to copy behavior blindly.
Differences found during implementation should be recorded in the discovery
log with a deliberate AFNI-compatible or corrected decision.

| Area | sumaru source | AFNI reference source |
|------|---------------|-----------------------|
| Statistics | `src/stats.rs` | `thd_statpval.c`, `mri_stats.c`, `nifti/nifticdf/nifticdf.c`, `p2dsetstat.c`, `dsetstat2p.c` |
| FDR/MDF | `src/dataset.rs::AfniFdrCurve` | `mri_fdrize.c`, `thd_fdrcurve.c`, `mri_floatvec.c`, `fdrval.c` |
| Dataset model | `src/dataset.rs` | `suma_datasets.c`, `suma_datasets.h`, `3ddata.h` |
| Colors | `src/color.rs` | `display.c::DC_spectrum_AJJ`, `display.c` bigmaps, `SUMA/SUMA_Color.c` |
| Overlays | `src/overlay.rs`, `viewer/mesh.rs` | `SUMA/SUMA_Color.c`, `SUMA/SUMA_Load_Surface_Object.c`, `afni_func.c` |
| Surface/domain | `src/surface.rs` | `SUMA/SUMA_Load_Surface_Object.c`, `SUMA/SUMA_GeomComp.c`, `SUMA/SUMA_Macros.h` |
| Surface clusters | `src/cluster.rs` | `SUMA/SUMA_SurfClust.c`, `SUMA/SUMA_SurfClust.h` |
| Volume clusters | none yet | `ptaylor/3dClusterize.c`, `mri_clusterize.c`, `afni_cluster.c` |
| ROIs | `src/roi.rs`, `viewer/roi.rs` | `SUMA/SUMA_input.c`, `SUMA/SUMA_niml.c`, `SUMA/SUMA_ROI2dataset.c` |
| Seed correlation | `src/instacorr.rs` | `thd_instacorr.c`, `thd_bandpass.c`, `SUMA/SUMA_dot.c` |
| Graphs/tracts | `src/graph_dataset.rs`, `src/tractography.rs` | FATCAT `Graph_Bucket` handling, `TrackIO.h`, relevant SUMA dataset code |

`afni-io` already supplies the format-side inputs: `StatKind`, `StatSpec`,
`ThresholdCurve`, typed arrays, `NimlDataset`, `LabelTable`, `NodeRoi`, volume
geometry, and triangular `Surface` data. Phase 0 must decide which semantic
types move into `afni-core` and which remain raw I/O types with conversion
adapters.

---

## Phase 0 — Workspace, dependency direction, and conformance harness ✅

- [x] ~~Create a Cargo workspace~~ **Decided against:** the crates stay separate
  projects (own `Cargo.toml`/`target`, `afni-io` its own git repo). Dependency
  direction is enforced by a test instead; see `docs/ARCHITECTURE.md`.
- [x] Write crate-level scope documentation and enforce the boundaries above.
- [x] Choose the migration plan (move each type in the phase that needs it;
  `afni-io` re-exports; see `docs/ARCHITECTURE.md` §3) for types currently public from `afni-io`:
  - move generally semantic types into `afni-core` and re-export them from
    `afni-io`, or
  - temporarily keep conversion types in both crates with explicit `From` /
    `TryFrom` adapters and a deprecation plan.
- [x] Avoid an `afni-core -> afni-io -> afni-core` cycle. The final dependency
  direction must be one-way.
- [x] Establish shared numeric/error conventions: `f64` for statistical math,
  `f32` for stored display/mesh buffers, checked index conversion, and explicit
  handling of NaN/Inf.
- [x] Add a reference-test harness that can run AFNI commands when available
  but retains committed expected values so normal tests do not require AFNI on
  `PATH`.
- [x] Record the AFNI version/commit used to generate conformance fixtures.
- [x] Add property tests or dense table tests (helpers now; sparse/dense and
  CPU/GPU equivalence cases arrive with Phases 1 and 5) for monotonicity, round trips,
  sparse/dense equivalence, and CPU/GPU-ready output equivalence.

Definition of done: both crates build independently, the dependency direction
is documented and enforced, and later phases have one standard way to compare
against AFNI.

---

## Phase 1 — File-neutral domains and datasets ✅

Start from `sumaru/src/dataset.rs`, but make the model useful for both surfaces
and volumes rather than encoding one viewer's state.

- [x] Define stable domain identifiers and explicit sample domains:
  `SurfaceDomain` (node count/topology identity) and `VolumeDomain`
  (dimensions/affine or grid identity).
- [x] Define dense and indexed row-to-sample mappings. Preserve sparse surface
  node lists without expanding them until a consumer asks for a dense view.
- [x] Define typed column storage for integer, float, and text data without
  coercing label keys or double-precision values to `f32`.
- [x] Port `Dataset`, `DataColumn`, `ColumnRange`, `ColumnRole`, time step/start,
  parent IDs, and validation from `sumaru/src/dataset.rs`.
- [x] Store `StatSpec`, FDR/MDF curves, units, and label tables as typed column
  metadata rather than reparsing strings in consumers.
- [x] Separate recorded range metadata from a freshly computed finite-data
  range; preserve both when they disagree.
- [x] Define policies for NaN/Inf, empty columns, duplicate sparse indices,
  out-of-domain indices, and missing samples.
- [x] Add adapters (`afni_io::adapt`) from:
  - `afni_io::NimlDataset`, including `COLMS_TYPE`, `COLMS_STATSYM`, time step,
    sparse indices, parents, FDR/MDF attributes, and label tables;
  - GIfTI data arrays and their intents/metadata;
  - AFNI/NIfTI volume sub-bricks through `afni_io::Volume`;
  - `.1D` tables where a domain is supplied by the caller.
- [x] Add the reverse adapter needed to write a core surface dataset through
  `afni-io`, while keeping unknown raw attributes in an I/O-side envelope.
- [x] Test dense/sparse equivalence and conversion of every committed
  statistical, label, and time-series fixture.

Open design decision: a raw file can contain unknown metadata that must
round-trip, while a core dataset should stay format-neutral. Prefer an
I/O-side wrapper containing `{ core_dataset, raw_extra_metadata }` over putting
NIML/HEAD attributes into `afni-core`.

---

## Phase 2 — Statistical interpretation and p-values ✅

`afni-io` currently reads statistical metadata but intentionally does not
interpret it. This phase supplies that missing layer.

### API and semantics

- [x] Add explicit tail selection:

  ```rust
  pub enum Tail { Lower, Upper, TwoSided }
  ```

  Do not use a global `two_sided_p_value` API: t/correlation/z are commonly
  two-sided, while F, chi-square, beta, binomial, gamma, and Poisson use an
  upper or lower tail according to the question.
- [x] Return a probability representation that retains `ln(p)` as well as a
  best-effort `p`, so very small probabilities can be displayed without
  underflow.
- [x] Provide checked forward and inverse operations, conceptually:
  `p_value(spec, statistic, tail)` and
  `critical_value(spec, probability, tail)`.
- [x] Distinguish a threshold probability (`p <= ...` for all samples that
  pass) from the exact probability of one sample (`p = ...`) at the UI edge.
- [x] Validate parameter counts, finiteness, domains, scale parameters, and
  degrees of freedom before calculation.
- [x] Never take `abs(value)` indiscriminately. It is appropriate for a
  two-sided symmetric null, but wrong for one-sided, shifted, noncentral, and
  direct-p distributions.

### Distribution coverage

- [x] First port and harden sumaru's correlation, t, F, z, and chi-square
  calculations.
- [x] Complete AFNI classic stat codes 2–10 by adding beta, binomial, gamma,
  and Poisson, matching `THD_stat_to_pval` / `THD_pval_to_stat`.
- [x] Implement direct probability intents:
  - `Pval`: the stored value is `p` and must be in `[0,1]`;
  - `LogPval`: `p = exp(-abs(value))`;
  - `Log10Pval`: `p = 10^(-abs(value))`.
- [x] Correct the semantic documentation for `LogPval` and `Log10Pval`:
  NIfTI permits a signed log representation but its reference library emits
  positive `-ln(p)` and `-log10(p)` values.
- [x] Add central Normal, Logistic, Laplace, Uniform, Weibull, Chi, inverse
  Gaussian, and extreme-value distributions.
- [x] Add noncentral t, F, and chi-square only with a numerically reliable
  implementation and AFNI/NIfTI conformance data; do not approximate them with
  their central forms. Chi-square and F are Poisson mixtures of the central
  forms; t is Lenth's AS 243 with a direct-integral fallback for the far tail
  where the series cancels. Both tails stay accurate in log space except as
  noted in the discovery log.
- [x] Give discrete inverse distributions a documented quantile convention.
  Generic floating-point bisection over `[0, +infinity)` is not sufficient.

### Correlation metadata ambiguity

- [x] Resolve AFNI versus NIfTI correlation parameters before exposing the
  general API:
  - AFNI `Correl(samples, nfit, nort)` uses three parameters;
  - NIfTI `NIFTI_INTENT_CORREL` defines one parameter, degrees of freedom;
  - AFNI's NIfTI/GIfTI paths often copy `intent_p1..3` as though they were AFNI
    parameters.
- [x] Represent those forms distinctly or carry provenance. The present
  `StatSpec::from_nifti_intent` zero-fills the AFNI three-parameter shape and
  cannot reliably interpret a standards-compliant third-party correlation
  intent.
- [x] Add fixtures for an AFNI correlation sub-brick, AFNI-written GIfTI, and a
  standards-compliant one-parameter NIfTI correlation intent.

### Validation

- [x] Compare common forward/inverse values against `cdf`, `ccalc`,
  `p2dsetstat` (`dsetstat2p` shares its code path), including explicit one- and two-sided cases.
- [x] Compare all NIfTI intent codes against `nifticdf.c` reference tables.
- [x] Test boundary probabilities, fractional degrees of freedom, invalid
  parameters, negative statistics, tiny tails, and forward/inverse round trips.
- [x] Document intentional departures from legacy AFNI sentinel behavior such
  as returning `0`, `-1`, or `99.99` on errors; Rust APIs should return errors.

---

## Phase 3 — FDR/MDF curves and multiple-comparison helpers ✅

- [x] Move the semantic `ThresholdCurve` into `afni-core` or add a lossless
  core equivalent with validated `x0`, nonzero `dx`, and finite samples.
  (Core has the validated `ThresholdCurve`; `afni-io` keeps its raw, unvalidated
  one for round trips, with adapters between them.)
- [x] Port AFNI's clamped four-point interpolation from `mri_floatvec.c`
  exactly. Sumaru currently uses a clamped Catmull–Rom polynomial, which is
  similar but not the same interpolation.
- [x] Implement statistic-threshold to `z(q)` and `q`, plus inverse `q` to
  threshold, matching `thd_fdrcurve.c` edge behavior.
- [x] Implement MDF lookup and document that its x-axis is `log10(p)`, unlike
  the statistic-threshold x-axis of an FDR curve.
- [x] Port and verify sumaru's automatic FDR-curve construction against
  `mri_fdrize.c`, including ignored zeros/non-finite values, minimum sample
  counts, true-positive estimation, q floors, and the 101-sample curve.
- [x] Provide a simple Benjamini-Hochberg helper for callers that have raw
  p-values but no stored curve; keep it distinct from AFNI's stored-curve
  algorithm.
- [x] Keep p-values and q-values separate in all names and display-facing
  results.
- [x] Compare interpolation and generated curves point-for-point with AFNI's
  `fdrval`, `3dFDR`, and committed `FDRCURVE_*` fixtures. Curves match to 1e-5,
  `3dFDR` z-scores match exactly (difference 0), `fdrval` to its 5 printed digits.

---

## Phase 4 — Colors, continuous maps, and label tables ✅

- [x] Port `Rgba`, `ColorStop`, `ContinuousColorMap`, and label-color lookup
  from `sumaru/src/color.rs` without `anyhow` or viewer dependencies.
- [x] Reconcile `afni_io::labels::LabelTable` with the display label table:
  preserve file order and attributes in I/O, while core provides validated
  keys, names, optional colors, lookup, and fallback policies.
- [x] Port AFNI's built-in color maps exactly, including byte endpoints and
  central gaps (all nine `display.c` scales at five sizes, byte for byte): `DC_spectrum_AJJ`, `DC_spectrum_ZSS`, and `display.c` bigmaps.
- [x] Add explicit interpolation modes (continuous, stepped/direct, nearest)
  and define duplicate-stop behavior.
- [x] Define whether interpolation occurs in encoded RGB or linear-light RGB;
  use AFNI/SUMA encoded-RGB behavior for parity and name any perceptual option
  separately.
- [x] Port stable fallback label colors, key-zero/unlabeled policy, alpha
  handling, brightness scaling, and NaN/missing colors.
- [x] Support building a color map from AFNI/SUMA label tables without losing
  integer keys or silently quantizing colors.
- [x] Generate compact golden tables from AFNI for every built-in map and test
  all 256 entries where AFNI uses a 256-entry map. (Tables come from AFNI's own
  C functions, extracted verbatim; see the discovery log for why not the binary.)

---

## Phase 5 — Thresholding, overlay evaluation, and compositing ✅

Port the pure parts of `sumaru/src/overlay.rs`; do not port viewer caches or
GPU resources as core state.

- [x] Define a validated `Threshold` with Off, Above, Below, Between, and
  Outside modes, including exact inclusive-boundary behavior.
- [x] Represent symmetric/absolute thresholding explicitly instead of relying
  on callers to manufacture `[-T, T]` ranges.
- [x] Define range selection (recorded/computed auto range versus manual),
  symmetric intensity ranges, clipping, zero display, opacity, brightness
  columns, and missing/non-finite policies.
- [x] Define `OverlaySpec` as immutable display intent and return an
  `OverlayEvaluation` containing colors, pass/fail mask, and diagnostics.
  Consumers may cache it; the core model should not hide invalidation inside a
  mutable viewer cache.
- [x] Port continuous and discrete-label mapping, sparse row-to-sample
  expansion, brightness modulation, threshold mask modes, and cluster masking.
- [x] Separate AFNI-compatible fade behavior from sumaru enhancements such as
  desaturation, darkening, and saturation boost.
- [x] Add pure alpha compositing for anatomical underlay, ordered overlay
  planes, live RGBA overlays, and ROI annotations. Specify straight versus
  premultiplied alpha.
- [x] Add threshold transfer by matched p-value as a core helper: preserve the
  source tail choice and convert through the destination `StatSpec`.
- [x] Produce GPU-ready lookup tables/uniform values, but keep shader code and
  buffers in the viewer.
- [ ] Test CPU evaluation against SUMA for threshold modes, `shw_0`, intensity
  clipping, brightness, direct color mapping, and alpha. Then use the same
  cases as shader conformance tests in each viewer. **Partly done:** 54 cases
  replayed against SUMA's own mapping code (`ScaleToMap`) cover interpolated,
  banded and direct mapping, clipping, auto range, mask ranges, `shw_0`-style zero
  masking and the brightness factor, node by node. Threshold-mode boundaries and
  both alpha-fade formulas are checked against the C source by unit tests, NOT
  against a running SUMA/AFNI (no command-line path exercises them). The committed
  `scaletomap.ref` is ready to be reused as the shader conformance set.

---

## Phase 6 — Surface topology, geometry metrics, and clustering ✅

Surface clustering currently depends on substantial reusable work in
`sumaru/src/surface.rs`. Extract that foundation before moving the cluster
algorithm.

- [x] Define a validated triangle topology: node/face counts, node neighbors,
  face neighbors, edges, boundary/non-manifold diagnostics, and stable topology
  identity.
- [x] Compute face normals, vertex normals, face areas, per-node areas, total
  area, bounds, and edge lengths without a graphics dependency.
- [x] Provide connected components, ring neighborhoods, and bounded geodesic
  graph searches reusable by clustering and ROI tools.
- [x] Port `ClusterParams`, `ClusterInput`, `ClusterLabels`, and summaries from
  `sumaru/src/cluster.rs`.
- [x] Support minimum node count and minimum surface area, edge-ring and
  millimeter-radius connectivity, largest-first labels, and deterministic tie
  breaking.
- [x] Keep positive and negative tails separate when requested. Document that
  sumaru's bisided behavior intentionally differs from `SurfClust`, which can
  merge touching opposite signs after masking.
- [x] Add peak node/value, min/max, area, node count, centroid/center of mass,
  and optional coordinate summaries useful to viewers and reports.
- [x] Compare against `SurfClust` for every mode it can express and retain
  separate tests for corrected bisided semantics. (22 runs: every edge-ring case
  matches exactly, rank by rank on the irregular mesh and as a set on the regular
  one; five of six millimetre cases match; the sixth differs in a documented
  direction. Central-node columns are not computed.)
- [x] Benchmark wide ring/radius searches; precompute reusable neighborhoods
  or weighted adjacency when repeated interactive clustering warrants it.

---

## Phase 7 — Volume neighborhoods and clustering ⬜

This phase supports the future slice viewer and avoids making the surface
cluster API masquerade as a volume algorithm.

- [ ] Define checked 3D index/coordinate conversion over a `VolumeDomain`.
- [ ] Implement AFNI NN1/NN2/NN3 voxel connectivity and explicit face/edge/
  corner neighborhood iteration.
- [ ] Support minimum voxel count and physical volume using voxel dimensions
  or affine-derived voxel volume.
- [ ] Support one-sided, two-sided/bisided, within-range, and outside-range
  threshold masks with deterministic cluster ordering.
- [ ] Report peak voxel/value, voxel count, physical volume, index centroid,
  and transformed world-coordinate centroid.
- [ ] Match `3dClusterize` and `mri_clusterize.c` on committed synthetic
  volumes, including oblique geometry where only reported coordinates—not
  connectivity—use the affine.
- [ ] Keep cluster simulation/correction tables out of the first pass; add a
  separate phase later if `3dClustSim` compatibility becomes a requirement.

---

## Phase 8 — File-neutral ROI model and operations ⬜

- [ ] Port the semantic model from `sumaru/src/roi.rs`: identity, parent
  domain/surface, side, label key, appearance, creation/edit status, drawing
  type, provenance, and ordered stroke/fill data.
- [ ] Convert losslessly to/from `afni_io::NodeRoi`, preserving unknown numeric
  codes in the I/O envelope even when core cannot interpret them.
- [ ] Provide canonical unique-node sets, validation against a surface domain,
  node-range summaries, label-table entries, and dataset conversion.
- [ ] Add topology-based open/closed path validation, shortest-path joining,
  boundary extraction, and filled-area construction as pure operations.
- [ ] Define union, intersection, difference, dilation/erosion by rings or
  geodesic distance, and connected-component cleanup.
- [ ] Preserve edit operations as explicit commands suitable for undo/redo;
  keep mouse/picking gestures in the viewer.
- [ ] Compare NIML replay and ROI-to-dataset output against
  `SUMA_NIMLDrawnROI_to_DrawnROI`, `SUMA_ROI2dataset`, and compact committed
  fixtures.

---

## Phase 9 — Time-series preprocessing and seed correlation ⬜

- [ ] Port `InstaCorrOptions`, prepared row storage, invalid-row masks, and
  seed correlation from `sumaru/src/instacorr.rs`.
- [ ] Separate generic signal operations—demean, detrend, Legendre bases,
  projection, normalization, FFT bandpass—from the InstaCorr orchestration.
- [ ] Match `THD_bandpass_vectors`, `thd_instacorr.c`, and
  `SUMA/SUMA_dot.c`, including `normalize_dset`, `polort`, filtered nuisance
  regressors, removed degrees of freedom, Nyquist clipping, and short-series
  rejection.
- [ ] Return derived correlation data with correct statistical metadata rather
  than only a `Vec<f32>`; the correlation degrees of freedom must reflect
  removed regressors/filtering.
- [ ] Support a single node seed first, then averaged ROI/multi-node seeds with
  a clearly defined normalization order.
- [ ] Feature-gate FFT dependencies if the base crate otherwise needs none.
- [ ] Test against AFNI/SUMA prepared vectors and output correlations, not just
  internal Rust reference calculations.

---

## Phase 10 — Graph, tract, and other derived semantic models ⬜

Only models and algorithms shared by more than one consumer belong here; file
syntax remains in `afni-io`.

- [ ] Extract the file-neutral portion of `GraphDataset`: nodes, labels,
  coordinates, full/triangular/sparse edge layouts, measures, ranges, and
  efficient matrix materialization.
- [ ] Put `Graph_Bucket` parsing/writing in `afni-io` and conversion/validation
  in adapters.
- [ ] Extract tract bundles, tract points, bounds, and selection/filtering;
  keep `TAYLOR_TRACT_DATUM` parsing in `afni-io`.
- [ ] Reuse core color/range/threshold types for graph edges and tract scalar
  attributes instead of inventing viewer-local variants.
- [ ] Consider reusable resampling/domain-mapping primitives only after their
  coordinate and interpolation semantics are specified. Do not infer that
  equal node counts imply compatible surfaces except under an explicit
  standard-template policy.
- [ ] Add other semantic models—surface states, annotations, time courses—only
  when both a parser and at least one nontrivial consumer need them.

---

## Phase 11 — Consumer migration, performance, and API stabilization ⬜

- [ ] Migrate `sumaru` in this order: dataset adapters → stats/FDR → colors →
  overlays → topology/clusters → ROI → InstaCorr. Keep compatibility facades
  so UI work is not mixed into algorithm extraction.
- [ ] Bring up `afni-view` against the same dataset/stat/color/overlay APIs and
  add volume clustering only after its volume model is stable.
- [ ] Remove duplicate `sumaru` implementations only after fixture and visual
  behavior tests pass.
- [ ] Benchmark large surface datasets, repeated threshold changes, FDR lookup,
  clustering, and seed correlation. Avoid allocations proportional to all
  samples on every slider movement when a reusable prepared form is possible.
- [ ] Stabilize names around statistical tail, threshold, sample domain, and
  sparse mappings; add serialization only for application state that is truly
  format-neutral.
- [ ] Publish crate documentation with small end-to-end examples using data
  supplied by a caller rather than reading files inside `afni-core`.
- [ ] Establish semantic versioning and a deprecation window for types moved
  from `afni-io` or `sumaru`.

---

## Suggested first usable milestone

The smallest milestone that materially helps both viewers is Phases 0–5 with
limited distribution coverage:

1. A validated surface `Dataset` adapter from `afni_io::NimlDataset`/GIfTI.
2. Tail-aware correlation, t, F, z, chi-square, and direct p-value intents.
3. Exact stored FDR-curve lookup.
4. AFNI continuous maps and label colors.
5. Pure threshold/overlay evaluation producing a dense RGBA buffer and pass
   mask.

That replaces sumaru's duplicated statistical and overlay interpretation while
leaving topology, ROI editing, clustering, and signal processing for later
increments.

## Definition of done for afni-core 1.0

- `afni-core` contains no file parsing, GUI, GPU, socket, or process-launching
  code.
- `afni-io` can adapt every supported statistical volume/surface dataset into
  core semantic types without losing typed values, sparse mappings, stat
  metadata, FDR curves, or label tables.
- `sumaru` and `afni-view` use the same p/q, color, threshold, and clustering
  semantics.
- Supported AFNI behavior is covered by committed conformance data and optional
  live comparisons to the AFNI executables/source implementation.
- Unsupported statistical kinds or malformed metadata produce explicit errors,
  never plausible-looking guessed p-values.
- CPU reference behavior and any viewer GPU implementation share the same
  golden cases.

---

## Open decisions

None open. D1–D20 were all decided on 2026-10-02 and are recorded below as R6–R25 (IDs are
never reused). New questions get the next free `D` number and a short context / pros / cons
entry here. Standing constraint from the owner: sumaru behavior is not to change for now, so
sumaru-facing items were deferred. "Deferred" entries stay valid and can be reopened.

### Resolved

| ID | Phase | Decision | Where recorded |
|----|-------|----------|----------------|
| R1 | 1 | Keep a file's representation: an index list of 0..n-1 stays indexed (not collapsed to dense). | Phase 1 |
| R2 | 1 | Add NaN-aware equality (`eq_nan_aware`) alongside IEEE `==`. | Phase 1 |
| R3 | 1 | Malformed FDR/MDF curves: error by default, `SkipWithWarning` on request. | Phase 1 |
| R4 | 2 | Implement the deferred noncentral t, F and chi-square exactly. | Phase 2 (noncentral) |
| R5 | 2 | Correlation in files: classify by structure, error on malformed headers, allow an explicit writer override. | Phase 2 (correlation in files) |
| R6 | 2 | (was D1) Deprecate, with a warning but keep, the lenient `Option`-returning stat APIs. Done: `StatSpec::from_nifti_intent` and `afni-io`'s `DataArray::stat` carry `#[deprecated]`. | Phase 2 (lenient intent APIs) |
| R7 | 3 | (was D2) New FDR curves default to AFNI's exact tie order. | Phase 3 (tie order) |
| R8 | 3 | (was D3) Below the curve's range, take max abs(stat) from the data, not from stored statistics. | Phase 3 (max abs stat) |
| R9 | 3 | (was D4) Skip `3dFDR -old`. | Phase 3 (extended statistics) |
| R10 | 3 | (was D5) Keep the `fdr` fixture dump files for now. | Phase 3 (fixtures) |
| R11 | 4 | (was D7) Port SUMA's named colormaps, verify them against AFNI, and label sumaru's three non-AFNI maps (`fire`, `afni_p2_spanned`, `amber_monochrome`) as sumaru's own. Implementation pending (Phase 11 or earlier). | Phase 4 (not AFNI definitions) |
| R12 | 4 | (was D8) `afni-io` depends on `afni-core` through a git dependency tracking `main` for now; move to tagged releases later. Implementation pending: edit `afni-io/Cargo.toml` and document a `[patch]` override for local work. | Phase 4 (building afni-io alone) |
| R13 | 5 | (was D10) A NaN threshold value is hidden by default; `MissingThreshold::Show` gives AFNI parity. | Phase 5 (non-finite threshold values) |
| R14 | 6 | (was D18) Add an absolute-value-weighted center of mass as a separate field; keep SurfClust's signed one. Implementation pending. | Phase 6 (center of mass is fragile) |
| R15 | 4 | (was D6) Strict IEEE bytes are "AFNI's" color scale; the one-byte FMA difference stays documented. | Phase 4 (platform-dependent bytes) |
| R16 | 5 | (was D9) Deferred: keep core's `Outside` and `AbsoluteAbove` as they are; sumaru migrates later. | Phase 5 (thresholds disagree about their own boundary) |
| R17 | 5 | (was D11) Deferred: sumaru keeps its continuous stops; core offers both models. | Phase 5 (panes are not N stops at i/(N-1)) |
| R18 | 5 | (was D12) Deferred: sumaru's fade is unchanged; core offers `FadeModel::Afni`/`Suma` separately. | Phase 5 (alpha fades: AFNI and SUMA are different formulas) |
| R19 | 5 | (was D13) Left as is (negative values transparent under `Above`) until the Phase 11 check against AFNI. Revisit then. | Phase 5 (unpinned: AFNI's one-sided fade) |
| R20 | 5 | (was D14) Deferred: no non-linear colormap type until a map needs it. | Phase 5 (what ScaleToMap does and does not cover) |
| R21 | 6 | (was D15) Keep true shortest-path radius; no SUMA-compatible mode. Core is correct. | Phase 6 (SUMA's millimetre radius is not graph distance) |
| R22 | 6 | (was D16) Two-sided clustering keeps positive and negative clusters separate (`Separate`, sumaru's behavior). Core's current default of `Merged` (SurfClust) is unchanged; callers choose `Separate` explicitly. | Phase 6 (sumaru migration) |
| R23 | 6 | (was D17) Keep SurfClust's seed and tie order. | Phase 6 (SurfClust quirks) |
| R24 | 6 | (was D19) Deferred: no central-node columns. | Phase 6 (not covered) |
| R25 | 6 | (was D20) `TopologyId` stays triangle-order sensitive. | Phase 6 (topology decisions) |

---

## Discovery log

- **2026-10-02 · Planning.** `afni-io` already reads all 23 AFNI/NIfTI
  statistical kinds but intentionally leaves p-value math to `afni-core`.
- **2026-10-02 · Planning.** Sumaru's `AfniStatSpec` is a good seed for common
  distributions but supports only correlation, t, F, z, and chi-square and
  calls the mixed two-sided/upper-tail operation `two_sided_p_value`.
- **2026-10-02 · Planning.** AFNI's newer `p2dsetstat`/`dsetstat2p` interface
  requires the caller to choose one- versus two-sided testing; sidedness is not
  safely inferable from a stat code alone.
- **2026-10-02 · Planning.** AFNI correlation metadata and standards-compliant
  NIfTI correlation intents use different parameterizations. A unified stat
  code without provenance is insufficient for reliable p-values.
- **2026-10-02 · Planning.** NIfTI `LOGPVAL`/`LOG10PVAL` are interpreted safely
  as `exp(-abs(value))`/`10^-abs(value)`; the reference library emits positive
  `-log(p)` values.
- **2026-10-02 · Planning.** Sumaru's FDR lookup uses clamped Catmull–Rom;
  AFNI's `interp_floatvec` uses a different clamped four-point cubic. Exact q
  parity requires the AFNI polynomial.
- **2026-10-02 · Planning.** Sumaru's dataset, overlay, color, cluster, ROI,
  surface-domain, and InstaCorr modules contain substantial pure logic that can
  move to a shared crate. Viewer stack selection, GPU caches, commands, and UI
  controls should remain consumer-owned.
- **2026-10-02 · Planning.** Surface and volume clustering share threshold and
  reporting concepts but have different neighborhood/size semantics; they
  should share result vocabulary without being forced through one algorithm.
- **2026-10-02 · Phase 0.** Cargo requires workspace members to live below the
  workspace root, so a combined `afni-core`/`afni-io` workspace would need a
  `Cargo.toml` in the shared parent directory. The project owner prefers three
  independent projects, so no workspace exists.
- **2026-10-02 · Phase 0.** AFNI_26.2.08 `cdf -t2p` is **two-sided** for `fizt`,
  `fitt`, and `fico` (z=0 gives p=1) but **upper-tail** for `fift`/`fict`. This
  confirms the explicit-`Tail` design; Phase 2 must not inherit it implicitly.
- **2026-10-02 · Phase 0.** `cdf` prints only 6 significant digits, so
  fixtures from it support a ~2e-6 relative tolerance. Tighter conformance needs
  `ccalc`/`p2dsetstat` output or AFNI's C source values.
- **2026-10-02 · Phase 1.** AFNI-written "dense" surface datasets still contain an
  `INDEX_LIST` of `0..n-1`, so `afni-io` reports them as sparse and the surface
  node count is not recorded anywhere in the file. Adapters therefore require a
  caller-supplied `node_count` for any dataset with an index list (never
  guessed from `max(index)+1`). `SampleMap::is_identity` detects dense data in
  disguise.
- **2026-10-02 · Phase 1.** Derived `PartialEq` on a dataset holding NaN
  ("missing") samples is never equal to its own clone (IEEE: NaN != NaN).
  Compare values with a defined fill, or add a NaN-aware comparison if needed.
- **2026-10-02 · Phase 1.** `StatKind`/`StatSpec` moved to `afni-core` and are
  re-exported by `afni-io` (no API break). `afni-io`'s raw `ThresholdCurve` stays
  unvalidated for round-trip fidelity; adapters convert it to the validated core
  curve and reject malformed curves with an error.
- **2026-10-02 · Phase 1.** `COLMS_RANGE` is recomputed by the NIML writer, so a
  recorded range read from a file does not survive a write; the core dataset
  keeps it only as read-time provenance (`RangeReport`).
- **2026-10-02 · Phase 1.** Volume sub-bricks become `Float32` columns (scaled
  values), so integer label volumes lose their integer storage type; their label
  table and `Label` role are kept. Revisit if exact integer storage is needed.
- **2026-10-02 · Phase 1.** GIfTI `stat()` copies `intent_p1..3` as AFNI
  parameters; the correlation ambiguity remains for Phase 2.
- **2026-10-02 · Phase 2.** AFNI's gamma "scale" parameter is a **rate**
  (`x * scale` is passed to the incomplete gamma, as in CDFLIB and the NIfTI
  definition). `nifti_stats -q 2 GAMMA 3 2` equals `Q(3, 4)`. afni-core names it
  `rate`.
- **2026-10-02 · Phase 2.** AFNI/CDFLIB evaluate binomial and Poisson tails at
  non-integer statistics using a continuous extension (`P(X > 2.5)` for Poisson(3)
  differs from `P(X > 2)`). afni-core uses the true step function (`floor(x)`)
  and an integer quantile for the inverse (smallest `k` with `P(X > k) <= p`).
  Conformance compares only integer statistics.
- **2026-10-02 · Phase 2.** AFNI `correl_p2t` returns `0` when `nort < 1`, although
  `correl_t2p` accepts it; afni-core accepts `nort >= 0` in both directions.
- **2026-10-02 · Phase 2.** `p2dsetstat -2sided` on a correlation sub-brick with
  `nfit = 2` returns the **upper tail of R** (multiple correlation, R >= 0), not
  a two-sided value; `-1sided` doubles p first. afni-core models `nfit > 1` as a
  one-sided multiple correlation and reproduces both outputs
  (`Correl(30,2,1)`: p = 0.01 gives 0.537614, p = 0.02 gives 0.501569).
- **2026-10-02 · Phase 2.** `3dAFNItoNIFTI -pure` and AFNI's GIfTI writer copy
  `Correl(samples, nfit, nort)` into `intent_p1..3` (verified: 30, 2, 1). AFNI
  always has `nfit >= 1`, the NIfTI standard leaves p2 = p3 = 0, so the two are
  separable by structure. `StatSpec::from_intent` classifies by that rule, or
  takes an explicit `IntentOrigin`, and errors on anything else; the standard
  form `Correl(dof)` is stored as the equivalent `Correl(dof + 1, 1, 0)`.
- **2026-10-02 · Phase 2.** AFNI's CDFLIB inverse is only good to ~1e-8 absolute
  (`-1 0.5 CORREL 18` returns -5.0e-9 where the exact answer is 0), so inverse
  conformance uses an absolute tolerance of 5e-8.
- **2026-10-02 · Phase 2.** AFNI's uniform reports `q = 1 - u` outside the support
  (e.g. 1.25 for a statistic below the range); afni-core returns the clamped
  probability.
- **2026-10-02 · Phase 2.** `ccalc` prints six *decimal places* (not significant
  digits), and `cdf`/`ccalc` are two-sided for `fitt`/`fico`. Use
  `tests/data/regenerate_nifti_stats.sh` (full-precision `nifticdf`) when tight
  comparisons are needed.
- **2026-10-02 · Phase 2.** The inverse Gaussian upper tail is computed as
  `1 - cdf`, so it loses accuracy below about 1e-12 (AFNI does the same). All
  other distributions keep both tails accurate in log space.
- **2026-10-02 · Phase 2 (noncentral).** AFNI's CDFLIB noncentral chi-square and F
  (`cumchn`, `cumfnc`) are only accurate to about 1e-5 absolute: for `F(3,20; 2)`
  at 3 AFNI gives 0.8237480 while two independent implementations (afni-core and
  a separate pure-Python summation) give 0.8237555. Conformance for those two
  uses an absolute tolerance of 2e-5 and relies on committed independent values
  for exactness. Noncentral t agrees with AFNI to about 2e-10 forward and 2e-9 for
  inverses, including negative noncentrality. Inverses are verified by round trip.
- **2026-10-02 · Phase 2 (noncentral t).** AS 243's upper tail is a difference of
  two positive sums when the noncentrality is negative (or, mirrored, a lower
  tail with positive noncentrality), and cancels almost completely in the far
  tail on the wrong side of the shift (e.g. the lower tail of `nct(30, 4)` at t = -1).
  afni-core detects the loss of digits and integrates the definition instead
  (`E[Q(tW - delta)]`, trapezoid rule in `x = ln w`), cross-checked against the
  series wherever both are valid.
- **2026-10-02 · Phase 2 (correlation in files).** NIfTI headers written by AFNI
  carry no marker (`descrip` and `intent_name` are empty), but none is needed:
  AFNI always writes `nfit >= 1` in `intent_p2`, the standard leaves it 0, so
  every well-formed file of either kind is classified correctly. Only malformed
  headers remain, and those are now errors (`Volume::stats`,
  `DataArray::stat_with_origin`), never silently "no statistic". `AdaptOptions`
  can state the writer (`intent_origin`) or skip a bad statistic with a warning
  (`statistics: SkipWithWarning`). A correlation with `nfit > 1`, which AFNI
  does write (for example `Correl(30,2,1)`), is a multiple correlation and is
  one-sided.
- **2026-10-02 · Phase 3.** AFNI's stored-curve algorithms are in single precision
  and depend on details that are easy to lose. Reproduced on purpose: `f32`
  arrays, p floored at `1e-15` and rounded to `f32`, `qsmal`/`m1`/`qfac` logic,
  `0.1666667` (not 1/6) in the cubic weights, and AFNI's own quicksort. With all
  of that ported, `3dFDR` z-scores match with a worst difference of exactly 0.
- **2026-10-02 · Phase 3 (tie order).** Many samples share exactly the same z(q)
  (the flat top of the step-up procedure), and AFNI's unstable quicksort then
  decides which statistic starts the curve (`x0`, and so every sample position,
  shifts by up to ~1e-3 relative). A stable or sorted-by-statistic order gave
  `x0 = 0.000391` where AFNI has `0.000473`. `cs_sort_ff.c` is now ported line for
  line (`afni_qsort_float_float`) and is the default; `deterministic_ties: true`
  gives the reproducible ordering instead. **Decide later** whether new curves
  should default to the deterministic order.
- **2026-10-02 · Phase 3 (AFNI quirks copied, flagged for review).**
  (a) `interp_floatvec` returns the first sample for every `x` when the curve has
  only two samples (`itop <= 1`); copied.
  (b) `estimate_m1` bins `(int)((p - 0.15) * 20)`, which truncates toward zero, so
  p in (0.10, 0.15) is counted in bin 0 although the comment says 0.15..0.95;
  copied.
  (c) `qsmal` uses 0 as "unset", so a first small q found at rank 0 is lost;
  harmless in practice, copied.
  (d) `mri_fdr_curve` divides by `t2 - t1`, which is 0 only for a degenerate
  curve (all statistics equal); guarded here (fraction 0), AFNI would produce NaN.
  (e) `interp_inverse_floatvec` divides by `yp - ym` on a flat segment equal to
  `y`; guarded (left end).
  (f) `student_t2p` returns p = 1 for `dof < 1`; afni-core evaluates it.
- **2026-10-02 · Phase 3 (fdrval edge behaviour).** `fdrval -qinput` replaces
  `q <= 0` by `1e-9` and `q >= 1` by `0.99999` and nothing else, so `q = 1e-12`
  is *not* clamped. afni-core's `threshold_for_q` takes `q` in `(0, 1]` and treats
  `q = 0` as an error (infinite z) and `q = 1` as threshold 0. A first draft of
  the conformance test clamped `q < 1e-9` and failed, which is how this was found.
- **2026-10-02 · Phase 3 (max |stat|).** For a q below the curve's range AFNI uses
  `DSET_BSTAT_MAXABS`, the stored brick maximum. afni-core computes it from the
  column values (`DataColumn::threshold_for_q`), which can differ if a file's
  stored statistics are stale. Revisit if this shows up.
- **2026-10-02 · Phase 3 (extended statistics).** AFNI builds FDR curves only for
  the classic codes 2-10 (`FUNC_IS_STAT`). afni-core does the same by default and
  requires an explicit `FdrOptions::tail` for any other kind rather than
  inheriting an undefined AFNI rule. `3dFDR -old` (`flags & 1`: count p = 1
  voxels, skip the m0 correction) is **not** ported.
- **2026-10-02 · Phase 3 (fixtures).** AFNI's own `p2t` functions return sentinel
  statistics (`99.99` for t, `999.99` for F) for p below about 1e-6 to 1e-4, so the
  random fixture `fdr+orig` contains a few voxels at those values; they exercise
  the `1e-15` p floor and the z cap and are kept. The fixture data come from
  `jRandomDataset`, so regenerating `make_volume_fixtures.sh` changes the data and
  every reference derived from it together (documented in the script). The new
  fixtures add about 0.9 MB to afni-io's `tests/data` (the per-voxel dumps are the
  bulk; drop them if repository size matters).
- **2026-10-02 · Phase 3 (sumaru migration).** `sumaru`'s `AfniFdrCurve::from_statistics`
  follows the same outline but differs in ways that will change numbers when it is
  replaced: `f64` throughout, stable sorting (so different tie order and `x0`), no
  MDF curve, p not rounded to `f32`, and Catmull-Rom interpolation. Re-run its
  tests against `afni_core::fdr` and expect small, explainable differences.
- **2026-10-02 · Phase 3 (performance, for Phase 11).** `fdrize` evaluates one
  incomplete-beta p-value per sample (via `PValueEvaluator`); a 1e6-voxel map was
  not benchmarked. The sort, histogram and z(q) conversion are linear or
  `n log n`. Benchmark before wiring it to an interactive slider.
- **2026-10-02 · Phase 3 (process note).** A scratch cleanup using `rm -f *` was
  refused by the environment's safety check, so the fixture trial used fresh
  directories instead; nothing in the repository was affected.
- **2026-10-02 · Phase 4 (platform-dependent bytes).** AFNI's colorscale code is
  floating-point and platform sensitive. Compiled with clang's default
  `-ffp-contract=on` (fused multiply-add, the default on Apple silicon), six
  entries across 45 scale/size combinations differ by one byte from strict IEEE
  arithmetic, always at the 0/360-degree hue wrap where `60 - ii*(60/(n/2-1))`
  lands on `-7e-15` instead of `0`: `[255,0,11]` versus `[255,11,0]` (for example
  `Reds_and_Blues` at 64 entries, index 31). They occur only at sizes other than
  AFNI's 256 default (and 128). The golden file uses strict IEEE
  (`-ffp-contract=off`), which is what Rust computes and what any compiler without
  FMA produces. **Open question:** which behavior should be the "AFNI" one, given
  that an M-series AFNI build shows the FMA bytes at non-default sizes.
- **2026-10-02 · Phase 4 (golden tables are not from the binary).** `display.c`
  needs X11/Motif and the bigmaps are not exposed by any command-line tool
  (`MakeColorMap`/`ScaleToMap` serve SUMA's own named maps, a different list).
  The golden tables are therefore produced by compiling the exact AFNI functions
  (`mypow`, `DC_spectrum_AJJ`, `DC_spectrum_ZSS`, `NJ_bigmaps_init` and the macros
  they use) extracted from the source tree by `regenerate_afni_colorscales.sh`.
  It proves the port matches the *source*; it does not prove the shipped binary was
  built from that source.
- **2026-10-02 · Phase 4 (AFNI quirks copied).** (a) AFNI's `pow` is
  `exp(y*log(x))` with 0 for `x <= 0`, and channels are `(int)(255*pow + 0.5)`;
  ported exactly, because the math library's `pow` can differ by a byte.
  (b) `Reds_and_Blues` subtracts `NBIG_MTOP + 1` (not the half-way point) in its
  upper half, so its first upper entries use a negative offset and start above 240
  degrees. (c) `Reds_and_Blues_w_Green` paints only two entries (`n/2-1`, `n/2`)
  green; sumaru widened that band to eight so it stays visible, which is a viewer
  choice and is *not* in afni-core. (d) `DC_spectrum_AJJ` uses `s` where its mirror
  branch uses `sb` for the blue ramp; both are 250.
- **2026-10-02 · Phase 4 (orientation).** AFNI stores a color bar with index 0 at
  the TOP (the highest value). `AfniColorScale::table` keeps that order so it can be
  compared byte for byte; `to_color_map` flips it so position 0 is the lowest value.
  sumaru's `spectrum_red_to_blue` already used the flipped orientation.
- **2026-10-02 · Phase 4 (not AFNI definitions).** sumaru's `blue_white_red`,
  `fire`, `grayscale`, `afni_p2_spanned` and `amber_monochrome` were *not* ported:
  none is a `display.c` bigmap, and the provenance of `afni_p2_spanned` and `fire`
  is unclear (the stops look hand-chosen). `amber_monochrome` is attributed to
  `pbardefs.h` and should be checked against it. Only a plain grayscale ramp is in
  core, marked as not AFNI's. SUMA's own named colormaps (`RGYBR20`, `bgyr19`,
  `ngray20`, ... from `SUMA_Color.c`) are a separate list and are **not** yet
  ported; do that if the viewers need them.
- **2026-10-02 · Phase 4 (label colors).** The stable fallback palette (ten colors,
  key 0 gray, `|key|-1 mod 10`) is a sumaru design, not AFNI's; it is kept because
  both viewers need distinguishable regions for datasets with no table, but is
  named and documented as such. Unlabeled keys are transparent by default (as in
  sumaru). A label dataset's table colors come back bit-for-bit (colors read from
  8-bit data are `value/255`; nothing is re-quantized), and duplicate keys in a
  file are an error rather than a silent merge. sumaru's label keys were `i32`;
  core uses `i64`.
- **2026-10-02 · Phase 4 (interpolation decisions to confirm).** Interpolation
  defaults to encoded RGB, as AFNI/SUMA do. A linear-light option exists and is
  named as such. Duplicate stops make a hard edge, and at the exact shared position
  the *earlier* stop wins; sumaru's loop behaves the same, now stated and tested.
  `Stepped` takes the stop at or below the position, which is how an N-pane color
  bar behaves, but AFNI's exact value-to-pane rule (the index arithmetic in the
  pbar code) is **not** ported; Phase 5 must pin it against SUMA before overlays
  claim parity.
- **2026-10-02 · Phase 4 (building afni-io alone).** `afni-io` now depends on
  `afni-core` by relative path (`../afni-core`), and the two repositories are
  published separately. Until `afni-io` points at a git or registry dependency, a
  fresh clone of `afni-io` does not build; both must be checked out side by side.
  The published `afni-io` (at its "phase 6" commit) does not yet have the
  dependency, so this bites when the adapter work is committed. **Decide:** a git
  dependency in `afni-io/Cargo.toml`, publishing `afni-core`, or a workspace.
- **2026-10-02 · Phase 5 (thresholds disagree about their own boundary).** SUMA's
  modes are not consistent: `ABS_LESS_THAN` hides `-t < v < t` strictly, so
  `|v| = t` PASSES, while `OUTSIDE_RANGE` hides `[lo, hi]` inclusively, so the ends
  FAIL. sumaru's single `Outside` mode passes the ends (`v <= min || v >= max`),
  which is SUMA's absolute mode, not SUMA's outside-range mode. Core has both,
  honestly named (`AbsoluteAbove`, `Outside`) and tested at the boundary; sumaru's
  `Outside(-T, T)` should migrate to `AbsoluteAbove(T)`. SUMA's `LESS_THAN` is in
  fact "show `v >= t`" (`Above`); `Below` exists only in sumaru and is kept.
- **2026-10-02 · Phase 5 (non-finite threshold values).** AFNI and SUMA test
  thresholds with `<`/`>`, which are false for NaN, so a sample whose THRESHOLD value
  is NaN is never hidden and is drawn at full strength. Core hides it by default and
  offers `MissingThreshold::Show` to reproduce the reference programs. An
  intensity that is NaN is a different case (AFNI casts it to `int`, which is
  undefined behavior in C; SUMA never gets that far); core gives it a configurable
  `missing_color` (transparent by default; sumaru used gray 0.35).
- **2026-10-02 · Phase 5 (panes are not N stops at i/(N-1)).** SUMA places color `i`
  of an `N`-color map at `i/N` of the range and holds the last color from `(N-1)/N`
  to the top (`Vscl = p * N`, mix `i` with `i + 1`). A continuous map with `N` stops
  at `i/(N-1)` (what sumaru builds from a 256-entry AFNI table) differs by up to
  1/256 of the range. Core keeps both: `ContinuousColorMap` for general gradients and
  `ColorTable` + `PaneRule` for SUMA/AFNI parity. sumaru's overlay should use panes
  for the AFNI scales to match SUMA exactly.
- **2026-10-02 · Phase 5 (AFNI volumes index panes from the top).** AFNI's volume
  overlay uses `j = (int)(N/(top-bot) * (top - v))` from the top of the bar, SUMA's
  `i = (int)((v-bot)/(top-bot) * N)` from the bottom. They agree except exactly on a
  pane boundary, where AFNI lands one pane lower. Both are implemented
  (`PaneRule::AfniBanded`, `SumaBanded`) and the difference is a unit test. All of
  this arithmetic is `f32` because the references are `float`.
- **2026-10-02 · Phase 5 (alpha fades: AFNI and SUMA are different formulas).**
  AFNI (`AFNI_newnewfunc_overlay`): `255*((1-floor)*|v|/t)^k + 255*floor` as a byte,
  `rintf` (ties to even), CLAMPED TO 222, v == 0 rejected. SUMA
  (`alphaOpacitiesForOverlay`): `min(1, |v|/t)`, squared if quadratic; no byte, no
  222 cap, no floor, and it uses only `ThreshRange[0]`. Both are modeled
  (`FadeModel::Afni`, `Suma`); sumaru's cubic/quartic curves, `max_alpha`,
  desaturate/darken/boost are NOT AFNI and live in `FadeCurve::{Cubic,Quartic}` and
  `FadeStyle`, separate from the AFNI behavior. sumaru's default (`max_alpha` 0.85,
  desaturate 0.5, darken 0.35) is therefore a viewer style, not "AFNI-compatible".
- **2026-10-02 · Phase 5 (unpinned: AFNI's one-sided fade).** For `Above(t)` AFNI
  fades `0 < v < t`, rejects `v == 0`, and (with `thb = 0`) neither fades nor rejects a
  NEGATIVE threshold value, which then draws opaque unless the separate positive-only
  flag hides it. Core makes negative values transparent for `Above`. Whether that
  matches the GUI could not be established without running it. Revisit with
  DriveSuma/AFNI GUI automation in Phase 11.
- **2026-10-02 · Phase 5 (what ScaleToMap does and does not cover).** The 54
  reference cases pass exactly, but the CLI itself is narrower than the library:
  `-br` is limited to `(0, 1]` (the library allows `(0, 2]`; factors above 1 are unit
  tested only), zeros are masked by default, the default mask color is black unless
  `-msk` is given (the documented 0.3 applies only then), and `-apr` runs a
  separate function (`SUMA_ScaleToMap_alaAFNI`) that has the same arithmetic. NOT
  covered: color maps with per-color fractions (non-equal panes,
  `SUMA_Linearize_Color_Map`), `-anr`, `-perc_clp`, `top_frac`, "no color" rows
  (`-1 -1 -1`), and SUMA's coordinate-bias and contouring options. Those need a
  non-linear colormap type that core does not have yet.
- **2026-10-02 · Phase 5 (smaller notes).** (a) `Rgba::scaled_clamped` is an inherent
  method defined in `overlay.rs`, away from the other `Rgba` code; move it into
  `color.rs` when convenient. (b) Brightness modulation can push channels outside
  0..1 for factors above 1 or below 0; SUMA leaves that to the GL clamp, core clamps
  immediately. (c) `FadeStyle::max_alpha` is a sumaru ceiling applied to the factor
  for ANY fade model, so combining it with `FadeModel::Afni` gives min(222/255, cap).
  (d) The conformance fixture was easy to get wrong: the first `ScaleToMap` run
  produced empty cases because `-br 1.5` is rejected by the CLI and `set -e` hid it;
  the script now records an empty case instead of aborting, and the test asserts all
  54 cases are populated.
- **2026-10-02 · Phase 6 (SUMA's millimetre radius is not graph distance).**
  `SurfClust -rmm r` (r > 0) calls `SUMA_getoffsets2`, which builds breadth-first
  LAYERS and gives each node a distance through the previous layer. Its precursor
  selection looks buggy: it initializes `n_prec` to the node's first neighbor, then
  compares `OffVect[n_prec] + Seg` where `Seg` is a SQUARED length against a running
  minimum that mixes squared lengths and distances, so the stored distance is not the
  minimum and is usually too large. Effect: SUMA connects FEWER nodes than a true
  shortest-path search. Measured on the irregular 642-node sphere with
  `-athresh 2.2`: at `-rmm 14` SurfClust finds 35 clusters, true graph distance 32
  (and a simple corrected layered search also 32); at 6 and 9 mm and on the regular
  sphere at 8 mm they agree exactly. Exact parity would require reproducing the quirk,
  which also depends on SUMA's (unspecified) neighbor ordering, so afni-core uses true
  Dijkstra distances and the conformance test expects exactly this one disagreement
  and its direction (afni-core never finds MORE clusters). Decide whether a
  "SUMA-compatible radius" mode is worth the effort.
- **2026-10-02 · Phase 6 (edge rings match exactly).** `-rmm -N` (N edge layers)
  has none of that ambiguity: every ring case matches SurfClust in cluster count,
  node counts, areas, means, extremes and their nodes, variance, standard error,
  centroid and center of mass.
- **2026-10-02 · Phase 6 (tie order depends on float noise).** SurfClust sorts
  clusters by area with a stable selection sort, so equal areas keep discovery order,
  but on a REGULAR mesh many clusters have mathematically equal areas whose
  single-precision sums differ in the last bit, and that noise decides their order.
  It cannot be reproduced. afni-core sorts stably in discovery order (f64 sums over
  ascending members), and the regular-mesh conformance cases are compared as sets;
  on the irregular mesh the rank order matches exactly.
- **2026-10-02 · Phase 6 (SurfClust quirks).** (a) A node whose VALUE is exactly 0 is
  inactive whatever the threshold says; core makes that
  `exclude_zero_values` (off by default, on in the conformance cases). (b) Clusters
  start from the HIGHEST node index downward (the seed scan has no `break`), which
  fixes discovery order, tie order and `-sort_none` order; core copies it and
  records `seed_node`. (c) A negative `-amm2` with no `-n` means "at least that many
  nodes" and is mapped. (d) The help text names `-ir_range` but the program accepts
  `-in_range`. (e) Extremes ties go to the first-discovered node in SUMA and to the
  lowest node number here; the test data avoid exact ties.
- **2026-10-02 · Phase 6 (thresholds confirmed).** SurfClust's `-thresh` (`>=`),
  `-athresh` (`|v| >=`, boundary passes), `-in_range` (inclusive) and `-ex_range`
  (`v < lo || v > hi`, boundary fails) match `Threshold::{Above, AbsoluteAbove,
  Between, Outside}` exactly, which independently confirms the Phase 5 boundary
  rules (the live program had not been usable for that before).
- **2026-10-02 · Phase 6 (geometry conventions pinned).** Against `SurfaceMetrics`
  and `SurfMeasures`: triangle normal `(v1-v0)x(v2-v0)` (outward for the usual
  winding); NODE normal is the normalized sum of UNIT triangle normals, which was
  distinguished from area-weighted (differs by up to ~0.7 on the irregular mesh) and
  angle-weighted; node area is one third of the adjacent triangle areas. SUMA works
  in `float`, so areas agree to about 1e-6 relative and the tolerances are those of
  the reference.
- **2026-10-02 · Phase 6 (center of mass is fragile).** SurfClust's center of mass is
  `sum(v xyz)/sum(v)` with signed values, so a cluster whose values nearly cancel gives
  an absurd position (coordinates near -900 on a 50 mm sphere in the test data) and an
  exactly zero sum gives NaN here. This is SurfClust's definition and is reproduced;
  weighting by `|v|` would be more useful and could be offered as a separate field.
- **2026-10-02 · Phase 6 (not covered).** SurfClust's central-node and
  weighted-central-node columns (run with `-no_cent`; slow in SUMA) are not
  computed. File-writing options (`-out_roidset`, `-out_clusterdset`, `-out_fulllist`,
  `-prepend_node_index`) are I/O concerns. `-thresh_col` (threshold from another
  column) is supported by the API (`ClusterInput::tail_values`, the caller applies the
  threshold) but is not exercised against the live program. SurfClust prints 2-3
  decimals, which limits how tightly values can be compared.
- **2026-10-02 · Phase 6 (sumaru migration).** sumaru's `cluster.rs` differs in ways
  that will change results or speed: its default is bisided (core's default is
  `Merged`, SurfClust's); it chooses ONE size metric (area or nodes) where SurfClust
  applies both limits; it starts seeds from the LOWEST node index; it allocates an
  `O(nodes)` distance array per search and uses `Vec::contains` in its ring search
  (quadratic overall); and `surfclust_command` builds a command string, which is
  consumer tooling and was not ported. Node areas there were `f32`.
- **2026-10-02 · Phase 6 (topology decisions).** Only an empty node set or an
  out-of-range index is an error; degenerate, duplicate, non-manifold,
  inconsistently wound and bow-tie defects are REPORTED so a viewer can still open a
  damaged mesh (`require_clean` is the strict form). `TopologyId` is an FNV-1a hash of
  the node count and triangle list IN ORDER: reordering triangles changes it. Decide
  whether an order-independent identity is wanted. Bow-tie detection is quadratic in
  node degree (fine for meshes, noted).
- **2026-10-02 · Phase 6 (benchmark, release build, 2026-10-02, Apple silicon).**
  One clustering with neighborhoods computed on demand, then a cache built once and 20
  re-clusterings (about a third of the nodes active):
  nodes 40,962 | 1 ring 2.0 ms | 3 rings 21 ms | 6 rings 46 ms | radius 2 edges 7 ms |
  radius 5 edges 42 ms. Building the full cache costs about 2x one on-demand run, then
  each cached run is about 1/10 of an on-demand run for wide settings, so the cache
  pays off after roughly three re-runs; for 1 ring it is never worth it. Smaller meshes
  scale linearly (2,562 nodes: 0.3 to 4 ms). The earlier concern that sumaru's search
  would be quadratic is removed by the reusable `NeighborhoodSearcher`.
- **2026-10-02 · Phase 6 (a Phase 5 diagnostic corrected).** The integration test
  tying clusters to overlays showed `OverlayDiagnostics::rejected_by_cluster` counted
  every non-surviving row, including rows already hidden for failing the threshold.
  It now counts only rows that PASSED the threshold and lost their cluster (documented
  on the field).
- **2026-10-02 · Phase 2 (lenient intent APIs).** `afni_io::gifti::DataArray::stat()` and
  `afni_core::stat::StatSpec::from_nifti_intent` return `Option`, so a malformed
  correlation looks the same as "not a statistic". The checked forms
  (`DataArray::stat_with_origin`, `Volume::stats_with_origin`, `StatSpec::from_intent`)
  report the problem. The lenient ones are kept only so callers that cannot surface an
  error keep compiling; decide whether to deprecate them.
