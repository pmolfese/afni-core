# afni-core

File-neutral data models and algorithms for AFNI/SUMA tools: what decoded data
*means*: statistics and p-values, FDR q-values, color maps, thresholds, clusters,
and regions of interest. It is shared by the
[`sumaru`](https://github.com/pmolfese/sumaru) surface viewer, a future 2D slice
viewer (`afniru`), command-line tools, and tests.

Reading and writing files is **not** done here: that is
[`afni-io`](https://github.com/pmolfese/afni-io). The dependency points one way,
`afni-io` -> `afni-core`, and a test enforces that `afni-core` contains no file
parsing, GUI, GPU, socket, or process-launching code.

AFNI's C source is the behavioral reference. Where this crate matches AFNI the
tests compare against AFNI's own output; where it deliberately differs, the
difference is written down (see the discovery log in
[`afni-core_ROADMAP.md`](afni-core_ROADMAP.md)).

## Status

| Phase | What | State |
|---|---|---|
| 0 | Error/numeric conventions, conformance harness, dependency enforcement | done |
| 1 | Domains (surface/volume), sparse and dense datasets, typed columns, adapters | done |
| 2 | Tail-aware p-values and critical values for every AFNI statistic | done |
| 3 | FDR/MDF curves, q-values, Benjamini-Hochberg | done |
| 4 | Colors, continuous maps, AFNI's built-in scales, label colors | done |
| 5 | Thresholds, transparent thresholding, overlay evaluation, compositing | done |
| 6 | Mesh topology and geometry, surface clustering (SurfClust-compatible) | done |
| 7 | Voxel connectivity and volume clustering (`3dClusterize`-compatible) | done |
| 8 | ROI model, node-set operations, undoable edits, ROI-to-dataset | done |
| 9 | Detrend/bandpass/orts (`THD_bandpass_vectors`-exact), seed correlation with its statistic | done |
| 10 | Graph (network) and tract models, with `afni-io` readers/writers | done |
| 11 | Consumer migration, performance, API stabilization | planned |

The full plan, with a dated log of every surprise found along the way, is in
[`afni-core_ROADMAP.md`](afni-core_ROADMAP.md). Every place where results differ from
AFNI/SUMA, and why, is listed in
[`docs/DIFFERENCES_FROM_AFNI.md`](docs/DIFFERENCES_FROM_AFNI.md).

## Layout

| Path | Contents |
|---|---|
| `src/lib.rs` | Crate docs, scope, module list |
| `src/error.rs`, `src/numeric.rs` | One `Error` type; `f32`/`f64` rules, checked indices, NaN/Inf policy |
| `src/domain.rs`, `src/mapping.rs`, `src/column.rs`, `src/dataset.rs` | Surface/volume domains, dense/sparse row maps, typed columns, the `Dataset` |
| `src/stat.rs`, `src/stats.rs`, `src/special.rs` | `StatSpec`; p-values and critical values with an explicit `Tail`; log-space special functions |
| `src/curve.rs`, `src/fdr.rs` | Validated FDR/MDF curves with AFNI's interpolation; q-values, curve construction |
| `src/threshold.rs`, `src/overlay.rs`, `src/composite.rs` | Thresholds (exact boundaries, AFNI/SUMA fades, matched-p transfer); data + display spec to colors and a pass mask; alpha compositing |
| `src/volume_cluster.rs` | NN1/2/3 voxel connectivity, thresholds, size limits, ranked clusters with peak, centroid and bounding box (world coordinates through the affine) |
| `src/graph.rs`, `src/tract.rs` | Validated networks (full / triangular / sparse edge layouts, measures, ranges, thresholds) and tracts (length, tangents, bounds, selection) |
| `src/signal.rs`, `src/instacorr.rs` | Detrend, Legendre regressors, orts, L2 normalize, FFT bandpass (no dependency); SUMA-style seed correlation returning `Correl(samples, 1, removed_dof)` |
| `src/roi.rs`, `src/roi_ops.rs`, `src/roi_edit.rs` | The ROI model (lossless codes), `NodeSet`, ROI-to-dataset; grow/shrink/boundary/components/shortest path/fill; undoable edit commands |
| `src/topology.rs`, `src/mesh.rs`, `src/cluster.rs` | Validated triangle-mesh connectivity and diagnostics; normals, areas, distance searches; connected-cluster labeling |
| `src/color.rs`, `src/afni_colors.rs`, `src/suma_colormaps.rs`, `src/labels.rs` | `Rgba` and continuous maps; AFNI's nine built-in scales, exactly; SUMA's nine standard maps (`bgyr19`, `byr64`, ...); label tables and label colors |
| `tests/common/` | Shared AFNI-comparison test harness |
| `tests/data/conformance/` | Committed AFNI reference values, with the AFNI version that made them |
| `tests/data/regenerate_*.sh` | Scripts that rebuild those references from AFNI |
| `docs/ARCHITECTURE.md` | Dependency direction and design decisions, phase by phase |

## A taste

```rust
use afni_core::stat::StatSpec;
use afni_core::stats::{critical_value, p_value, Tail};
use afni_core::afni_colors::AfniColorScale;
use afni_core::column::ColumnRange;
use afni_core::color::Rgba;
use afni_core::overlay::{evaluate_rows, OverlayColors, OverlayInputs, OverlaySpec, RangeSelection};
use afni_core::threshold::Threshold;

// A t statistic with 23 degrees of freedom, as AFNI writes `Ttest(23)`.
let spec = StatSpec::parse("Ttest(23)").unwrap();
// The tail is always your choice; there is no hidden default.
let p = p_value(&spec, 2.5, Tail::TwoSided)?;
let t = critical_value(&spec, 0.01, Tail::TwoSided)?; // |t| for p = 0.01

// AFNI's red-to-blue scale, colored by value over a display range.
let scale = AfniColorScale::SpectrumRedToBlue.to_color_map(256)?;
let range = ColumnRange::new(-4.0, 4.0)?;
let color = scale.color_for_value(t, &range, Rgba::TRANSPARENT);

// A whole overlay: color values, hide what fails |value| >= 2, in one call.
let mut overlay = OverlaySpec::new(OverlayColors::Continuous(scale));
overlay.intensity_range = RangeSelection::Manual(range);
overlay.threshold = Threshold::AbsoluteAbove(2.0);
let values = [-3.5, -1.0, 0.5, 2.0, 3.9];
let result = evaluate_rows(&overlay, &OverlayInputs { intensity: &values, ..Default::default() })?;
assert_eq!(result.passed, [true, false, false, true, true]);
# Ok::<(), afni_core::Error>(())
```

## Design in brief

- **Explicit, never guessed.** The tail, the correlation parameterization, a
  surface's node count, and what to do with NaN are always stated by the caller;
  malformed metadata is an error, not a plausible-looking guess.
- **Exact where it can be.** Probabilities keep their logarithm, so p = 1e-400 is
  representable. AFNI's interpolation, FDR construction, color tables, and
  SUMA's value-to-color mapping are ported to reproduce AFNI's numbers, and the
  tests say how closely they do.
- **Deterministic.** Neighbor lists, edge lists, cluster ranks and tie-breaking are
  all defined, so the same input always gives the same output in the same order.
- **Small.** One dependency (`thiserror`), `unsafe` forbidden, Rust 1.74.

## Testing

```sh
cargo test                   # committed fixtures only; AFNI is not needed
AFNI_CORE_LIVE=1 cargo test  # also replay cases against AFNI on PATH
```

Some conformance tests need data files that live in `afni-io`'s test fixtures
(`afni-io/tests/data/`); those tests run in the `afni-io` repository, which depends
on this crate. Regenerating a reference needs AFNI (and for some, its source
tree):

```sh
tests/data/regenerate_conformance.sh         # cdf conventions
tests/data/regenerate_nifti_stats.sh         # nifticdf at full precision
tests/data/regenerate_afni_colorscales.sh    # display.c color scales
tests/data/regenerate_suma_colormaps.sh      # SUMA standard maps (MakeColorMap -std)
tests/data/regenerate_volume_clusters.sh     # 3dClusterize maps and reports
tests/data/regenerate_roi_refs.sh            # SurfDist distances, ROIgrow growth
tests/data/regenerate_signal_refs.sh         # THD_bandpass_vectors (builds a C harness on libmri)
tests/data/regenerate_instacorr_refs.sh      # SUMA InstaCorr steps with AFNI functions
tests/data/regenerate_scaletomap.sh          # SUMA's value-to-color mapping
tests/data/regenerate_surface_refs.sh        # SurfaceMetrics / SurfMeasures / SurfClust
cargo run --release --example cluster_bench  # clustering timings
```

## Building next to afni-io

`afni-core` builds on its own. `afni-io` currently depends on it by relative
path (`afni-core = { path = "../afni-core" }`), so to build or test `afni-io`,
keep both repositories checked out side by side:

```text
afni_rust/
  afni-core/
  afni-io/
```

(A cloned `afni-io` alone will not build until that dependency is changed to a
git or registry dependency; see the roadmap discovery log.)

## License

Public domain. United States Government work (17 U.S.C. § 105); outside the US,
rights are waived under [CC0 1.0](LICENSE).
