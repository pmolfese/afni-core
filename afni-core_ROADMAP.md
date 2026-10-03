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
| 3 | FDR/MDF curves and multiple-comparison helpers | ⬜ |
| 4 | Colors, continuous maps, and label tables | ⬜ |
| 5 | Thresholding, overlay evaluation, and compositing | ⬜ |
| 6 | Surface topology, geometry metrics, and clustering | ⬜ |
| 7 | Volume neighborhoods and clustering | ⬜ |
| 8 | File-neutral ROI model and operations | ⬜ |
| 9 | Time-series preprocessing and seed correlation | ⬜ |
| 10 | Graph, tract, and other derived semantic models | ⬜ |
| 11 | Consumer migration, performance, and API stabilization | ⬜ |

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

## Phase 3 — FDR/MDF curves and multiple-comparison helpers ⬜

- [ ] Move the semantic `ThresholdCurve` into `afni-core` or add a lossless
  core equivalent with validated `x0`, nonzero `dx`, and finite samples.
- [ ] Port AFNI's clamped four-point interpolation from `mri_floatvec.c`
  exactly. Sumaru currently uses a clamped Catmull–Rom polynomial, which is
  similar but not the same interpolation.
- [ ] Implement statistic-threshold to `z(q)` and `q`, plus inverse `q` to
  threshold, matching `thd_fdrcurve.c` edge behavior.
- [ ] Implement MDF lookup and document that its x-axis is `log10(p)`, unlike
  the statistic-threshold x-axis of an FDR curve.
- [ ] Port and verify sumaru's automatic FDR-curve construction against
  `mri_fdrize.c`, including ignored zeros/non-finite values, minimum sample
  counts, true-positive estimation, q floors, and the 101-sample curve.
- [ ] Provide a simple Benjamini-Hochberg helper for callers that have raw
  p-values but no stored curve; keep it distinct from AFNI's stored-curve
  algorithm.
- [ ] Keep p-values and q-values separate in all names and display-facing
  results.
- [ ] Compare interpolation and generated curves point-for-point with AFNI's
  `fdrval`, `3dFDR`, and committed `FDRCURVE_*` fixtures.

---

## Phase 4 — Colors, continuous maps, and label tables ⬜

- [ ] Port `Rgba`, `ColorStop`, `ContinuousColorMap`, and label-color lookup
  from `sumaru/src/color.rs` without `anyhow` or viewer dependencies.
- [ ] Reconcile `afni_io::labels::LabelTable` with the display label table:
  preserve file order and attributes in I/O, while core provides validated
  keys, names, optional colors, lookup, and fallback policies.
- [ ] Port AFNI's built-in color maps exactly, including byte endpoints and
  central gaps: `DC_spectrum_AJJ`, `DC_spectrum_ZSS`, and `display.c` bigmaps.
- [ ] Add explicit interpolation modes (continuous, stepped/direct, nearest)
  and define duplicate-stop behavior.
- [ ] Define whether interpolation occurs in encoded RGB or linear-light RGB;
  use AFNI/SUMA encoded-RGB behavior for parity and name any perceptual option
  separately.
- [ ] Port stable fallback label colors, key-zero/unlabeled policy, alpha
  handling, brightness scaling, and NaN/missing colors.
- [ ] Support building a color map from AFNI/SUMA label tables without losing
  integer keys or silently quantizing colors.
- [ ] Generate compact golden tables from AFNI for every built-in map and test
  all 256 entries where AFNI uses a 256-entry map.

---

## Phase 5 — Thresholding, overlay evaluation, and compositing ⬜

Port the pure parts of `sumaru/src/overlay.rs`; do not port viewer caches or
GPU resources as core state.

- [ ] Define a validated `Threshold` with Off, Above, Below, Between, and
  Outside modes, including exact inclusive-boundary behavior.
- [ ] Represent symmetric/absolute thresholding explicitly instead of relying
  on callers to manufacture `[-T, T]` ranges.
- [ ] Define range selection (recorded/computed auto range versus manual),
  symmetric intensity ranges, clipping, zero display, opacity, brightness
  columns, and missing/non-finite policies.
- [ ] Define `OverlaySpec` as immutable display intent and return an
  `OverlayEvaluation` containing colors, pass/fail mask, and diagnostics.
  Consumers may cache it; the core model should not hide invalidation inside a
  mutable viewer cache.
- [ ] Port continuous and discrete-label mapping, sparse row-to-sample
  expansion, brightness modulation, threshold mask modes, and cluster masking.
- [ ] Separate AFNI-compatible fade behavior from sumaru enhancements such as
  desaturation, darkening, and saturation boost.
- [ ] Add pure alpha compositing for anatomical underlay, ordered overlay
  planes, live RGBA overlays, and ROI annotations. Specify straight versus
  premultiplied alpha.
- [ ] Add threshold transfer by matched p-value as a core helper: preserve the
  source tail choice and convert through the destination `StatSpec`.
- [ ] Produce GPU-ready lookup tables/uniform values, but keep shader code and
  buffers in the viewer.
- [ ] Test CPU evaluation against SUMA for threshold modes, `shw_0`, intensity
  clipping, brightness, direct color mapping, and alpha. Then use the same
  cases as shader conformance tests in each viewer.

---

## Phase 6 — Surface topology, geometry metrics, and clustering ⬜

Surface clustering currently depends on substantial reusable work in
`sumaru/src/surface.rs`. Extract that foundation before moving the cluster
algorithm.

- [ ] Define a validated triangle topology: node/face counts, node neighbors,
  face neighbors, edges, boundary/non-manifold diagnostics, and stable topology
  identity.
- [ ] Compute face normals, vertex normals, face areas, per-node areas, total
  area, bounds, and edge lengths without a graphics dependency.
- [ ] Provide connected components, ring neighborhoods, and bounded geodesic
  graph searches reusable by clustering and ROI tools.
- [ ] Port `ClusterParams`, `ClusterInput`, `ClusterLabels`, and summaries from
  `sumaru/src/cluster.rs`.
- [ ] Support minimum node count and minimum surface area, edge-ring and
  millimeter-radius connectivity, largest-first labels, and deterministic tie
  breaking.
- [ ] Keep positive and negative tails separate when requested. Document that
  sumaru's bisided behavior intentionally differs from `SurfClust`, which can
  merge touching opposite signs after masking.
- [ ] Add peak node/value, min/max, area, node count, centroid/center of mass,
  and optional coordinate summaries useful to viewers and reports.
- [ ] Compare against `SurfClust` for every mode it can express and retain
  separate tests for corrected bisided semantics.
- [ ] Benchmark wide ring/radius searches; precompute reusable neighborhoods
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
