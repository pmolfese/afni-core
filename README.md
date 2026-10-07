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
| 11 | Consumer migration, performance, API stabilization | in progress |

The full plan, with a dated log of every surprise found along the way, is in
[`afni-core_ROADMAP.md`](afni-core_ROADMAP.md). Every place where results differ from
AFNI/SUMA, and why, is listed in
[`docs/DIFFERENCES_FROM_AFNI.md`](docs/DIFFERENCES_FROM_AFNI.md).

## Layout

| Path | Contents |
|---|---|
| `src/lib.rs` | Crate docs, scope, module list |
| `src/error.rs`, `src/numeric.rs` | One `Error` type; `f32`/`f64` rules, checked indices, NaN/Inf policy |
| `src/domain.rs`, `src/mapping.rs`, `src/column.rs`, `src/dataset.rs` | Surface/volume domains, dense/sparse row maps, typed columns, the `Dataset`, and checked immutable transformations |
| `src/stat.rs`, `src/stats.rs`, `src/special.rs` | `StatSpec`; p-values and critical values with an explicit `Tail`; log-space special functions |
| `src/curve.rs`, `src/fdr.rs` | Validated FDR/MDF curves with AFNI's interpolation; q-values, curve construction |
| `src/threshold.rs`, `src/overlay.rs`, `src/composite.rs` | Thresholds (exact boundaries, AFNI/SUMA fades, matched-p transfer); data + display spec to colors and a pass mask; alpha compositing |
| `src/volume_cluster.rs` | NN1/2/3 voxel connectivity, thresholds, size limits, ranked clusters with peak, centroid and bounding box (world coordinates through the affine) |
| `src/graph.rs`, `src/tract.rs` | Validated networks (full / triangular / sparse edge layouts, measures, ranges, thresholds) and tracts (length, tangents, bounds, selection) |
| `src/signal.rs`, `src/instacorr.rs` | Detrend, Legendre regressors, orts, L2 normalize, explicit power spectra, FFT bandpass (no dependency); SUMA-style seed correlation returning `Correl(samples, 1, removed_dof)` |
| `src/timeseries.rs` | Allocation-free row/time-series views over column-major datasets, including reusable work buffers and sparse sample lookup |
| `src/mask.rs` | Domain-sized sample masks, numeric membership rules, sparse row mapping, and checked boolean composition |
| `src/processing.rs` | Mask-aware multi-column voxel/node mapping, time-series reduction, and same-length time-series transformation |
| `src/reduction.rs` | One-pass numeric column summaries: counts, extrema with sample locations, means, and population/sample variance |
| `src/roi.rs`, `src/roi_ops.rs`, `src/roi_edit.rs` | The ROI model (lossless codes), `NodeSet`, ROI-to-dataset; grow/shrink/boundary/components/shortest path/fill; undoable edit commands |
| `src/topology.rs`, `src/mesh.rs`, `src/cluster.rs` | Validated triangle-mesh connectivity and diagnostics; normals, areas, distance searches; connected-cluster labeling |
| `src/color.rs`, `src/afni_colors.rs`, `src/suma_colormaps.rs`, `src/labels.rs` | `Rgba` and continuous maps; AFNI's nine built-in scales, exactly, plus its default `Reds_and_Blues_Inv`; SUMA's nine standard maps (`bgyr19`, `byr64`, ...); label tables and label colors |
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

## Time-series access

A `Dataset` stores one column per time point so column types and metadata stay
intact. `TimeSeriesView` provides the row-oriented access expected by AFNI
programs without transposing the complete dataset:

```rust
use afni_core::{dataset::Dataset, Result};

fn process(dataset: &Dataset) -> Result<()> {
    let time = dataset.time_series()?;
    let mut work = vec![0.0; time.time_point_count()];

    for row in 0..time.series_count() {
        // Reuse one contiguous buffer for detrending, filtering, FFTs, etc.
        time.copy_series_into(row, &mut work)?;
        let domain_sample = time.sample_for_series(row)?;
        println!("sample {domain_sample}: {} time points", work.len());
    }
    Ok(())
}
```

For read-only calculations, `time.series(row)?` is an allocation-free iterator.
Sparse datasets retain their row order and expose the corresponding volume
voxel or surface-node index through `sample_for_series`.

## Sample masks

`SampleMask` replaces an ambiguous raw `&[bool]` with one boolean per sample
plus the exact surface or volume domain those booleans describe. A mask made
from a sparse dataset column expands into domain order; samples with no stored
row are unselected:

```rust
use afni_core::{
    dataset::Dataset,
    mask::{MaskRule, SampleMask},
    Result,
};

fn selected_rows(dataset: &Dataset) -> Result<Vec<usize>> {
    let mask = SampleMask::from_column(
        dataset,
        &dataset.columns()[0],
        MaskRule::NonZero,
    )?;

    // `rows` maps the domain-sized mask back through this dataset's SampleMap.
    Ok(mask
        .rows(dataset)?
        .enumerate()
        .filter_map(|(row, selected)| selected.then_some(row))
        .collect())
}
```

`and` and `or` reject structurally different domains even when their lengths
match. Preserve `DomainId` for surfaces when identity matters: two anonymous
surfaces with the same node count cannot be distinguished. `values()` provides
a temporary compatibility bridge for older APIs that still accept a
domain-ordered boolean slice.

## Volume grids and coordinates

`VolumeDomain` is the canonical semantic volume grid. It provides checked
i-fastest `linear_index`/`ijk` conversion and both affine directions:

```rust
use afni_core::domain::VolumeDomain;

fn round_trip(domain: &VolumeDomain, ijk: [f64; 3]) -> afni_core::Result<[f64; 3]> {
    let world = domain.ijk_to_world(ijk)?;
    domain.world_to_ijk(world)
}
```

`world_to_ijk` returns fractional, possibly out-of-grid coordinates and rejects
a missing or singular affine. `crop(start, dims)` and `pad(before, after)`
return grid descriptions with adjusted affine origins; they do not themselves
move or allocate voxel values. Since either operation changes the sample set,
the resulting domain intentionally clears its domain identifier.

## Spatial affine transforms

AFNI `aff12.1D` matrices are spatial mappings, not dataset voxel-to-world grid
affines. `AffineTransform` represents one validated, invertible mapping and
`AffineTransformSeries` represents a non-empty sequence such as one motion
transform per volume:

```rust
use afni_core::affine::{AffineTransform, AffineTransformSeries};

let shift = AffineTransform::from_aff12_row([
    1.0, 0.0, 0.0, 2.0,
    0.0, 1.0, 0.0, 3.0,
    0.0, 0.0, 1.0, 4.0,
])?;
assert_eq!(shift.apply_point([1.0, 1.0, 1.0])?, [3.0, 4.0, 5.0]);

let series = AffineTransformSeries::new(vec![shift.clone(), shift])?;
let inverse = series.inverse()?;
assert_eq!(inverse.len(), 2);
# Ok::<(), afni_core::Error>(())
```

`first.then(&second)` applies `first` and then `second`, matching AFNI
`cat_matvec` application order. Coordinate convention is part of the type;
RAI/RAS conversion is explicit. File parsing remains in `afni-io`.

## Safe dataset transformations

`Dataset` and `DataColumn` do not expose mutable fields that can break their
validated row/domain relationship. Transformations return a new dataset and
leave the source untouched if the closure fails or produces the wrong number
of rows:

```rust
use afni_core::{
    column::ValueMetadataPolicy,
    dataset::Dataset,
    Result,
};

fn double_first_column(dataset: &Dataset) -> Result<Dataset> {
    dataset.transform_column(0, |column| {
        column.map_numeric_to_f64(
            ValueMetadataPolicy::DiscardValueMetadata,
            |value| value * 2.0,
        )
    })
}
```

The metadata policy is deliberately required. `Preserve` is appropriate for a
representation-only change whose scientific meaning is identical;
`DiscardValueMetadata` clears units, statistic metadata, FDR/MDF curves, label
tables, and the recorded range after arithmetic. Column append, replacement,
selection/reordering, and removal are also checked, and `column_at` reports an
error rather than panicking on a bad index.

## Voxelwise and time-series processing

The processing helpers combine checked datasets, domain-aware masks, and
reusable work buffers. This example detrends every selected voxel or node while
leaving unselected series unchanged:

```rust
use afni_core::{
    column::ValueMetadataPolicy,
    dataset::Dataset,
    mask::SampleMask,
    processing::{transform_time_series, UnselectedSeries},
    signal::Detrend,
    Result,
};

fn detrend_selected(dataset: &Dataset, mask: &SampleMask) -> Result<Dataset> {
    transform_time_series(
        dataset,
        Some(mask),
        UnselectedSeries::Preserve,
        ValueMetadataPolicy::DiscardValueMetadata,
        |_context, series| {
            Detrend::Linear.apply(series);
            Ok(())
        },
    )
}
```

Common 3dTstat-style summaries need no callback. Every requested statistic is
computed during the same pass over each selected series, and the result is a
new scalar `Dataset` on the same volume or surface domain:

```rust
use afni_core::{
    dataset::Dataset,
    mask::SampleMask,
    numeric::NonFinitePolicy,
    processing::TimeSeriesStatistic,
    reduction::VarianceNormalization,
    Result,
};

fn temporal_stats(dataset: &Dataset, mask: &SampleMask) -> Result<Dataset> {
    dataset.summarize_time_series(
        Some(mask),
        [
            TimeSeriesStatistic::Mean.output("mean", 0.0)?,
            TimeSeriesStatistic::StandardDeviation(VarianceNormalization::Sample)
                .output("stdev", 0.0)?,
            TimeSeriesStatistic::Slope.output("slope_per_second", 0.0)?,
        ],
        NonFinitePolicy::Skip,
    )
}
```

When command-line flags determine the outputs at runtime, collect the requests
in a `Vec` and use the dynamic form. It performs the same single-pass
accumulation without requiring a match on every possible output count:

```rust
use afni_core::{
    dataset::Dataset,
    numeric::NonFinitePolicy,
    processing::{TimeSeriesStatistic, TimeSeriesSummary},
    Result,
};

fn requested_stats(
    dataset: &Dataset,
    want_mean: bool,
    want_maximum: bool,
) -> Result<Dataset> {
    let mut requests: Vec<TimeSeriesSummary> = Vec::new();
    if want_mean {
        requests.push(TimeSeriesStatistic::Mean.output("mean", 0.0)?);
    }
    if want_maximum {
        requests.push(TimeSeriesStatistic::Maximum.output("maximum", 0.0)?);
    }
    dataset.summarize_time_series_dynamic(None, &requests, NonFinitePolicy::Skip)
}
```

The built-ins are count, sum, sum of squares, L2 norm, mean, population or
sample variance/standard deviation, RMS, slope per second, minimum, and
maximum. Slope requires a dataset time step. Undefined selected-row results are
NaN; each output explicitly declares the value written outside the mask.

`combine_columns` is the equivalent of a C voxel loop reading several
sub-bricks. For calculations not covered by the built-ins,
`summarize_time_series_with` returns one or several columns from one callback,
so an expensive fit or spectrum is performed only once per voxel/node.
`summarize_time_series_dynamic_with` provides the runtime-sized equivalent; its
callback fills a reusable output slice, avoiding a new allocation per sample.
Callbacks receive both the stored row and the full-domain sample index, which
keeps sparse datasets unambiguous. Mask domains, output descriptions, and
column types are validated before processing begins.

`Dataset::derive(kind, columns)` packages newly calculated columns as a related
dataset. It preserves the spatial domain, sparse mapping, and domain/geometry
parent identifiers while clearing the source object's own ID and time axis.
Use `replace_columns` when the result is still the same scientific dataset and
its existing kind and timing remain valid.

`signal::power_spectrum` uses the same arbitrary-length FFT as the bandpass
implementation. Its options make zero-padding, full versus one-sided output,
raw versus FFT-length-normalized power, and non-finite handling explicit.

## Column summaries

`Dataset::summarize_column` replaces the repeated C loop that counts usable
values and computes a masked range, mean, and variance. Both choices that can
silently change a result—non-finite handling and the variance denominator—are
required explicitly:

```rust
use afni_core::{
    dataset::Dataset,
    mask::SampleMask,
    numeric::NonFinitePolicy,
    reduction::{ColumnSummary, ColumnSummaryOptions, VarianceNormalization},
    Result,
};

fn masked_stats(dataset: &Dataset, mask: &SampleMask) -> Result<ColumnSummary> {
    dataset.summarize_column(
        0,
        Some(mask),
        ColumnSummaryOptions::new(
            NonFinitePolicy::Skip,
            VarianceNormalization::Sample,
        ),
    )
}
```

The result distinguishes selected, finite, non-finite, included, zero, and
nonzero counts. Empty selections produce `None` rather than a made-up zero for
the numeric reductions. Extrema retain the stored row and complete-domain
voxel or surface-node index, so sparse datasets remain unambiguous. Ties use
the first stored row.

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
tests/data/regenerate_afni_default_scale.sh  # AFNI default overlay scale (afni.c + pbardefs.h)
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
