# How afni-core differs from AFNI and SUMA

This is the running list of every place where `afni-core` does something other than
the AFNI/SUMA C code, and why. It exists so that nobody has to guess whether a number
from this crate should equal a number from AFNI: look the feature up here first.

**Keep it current.** Any change that makes results differ from AFNI/SUMA (a decision
you make, a quirk we copy, a quirk we cannot copy, something not ported) gets a row
here **in the same commit**, with its reference ID, and a dated entry in the roadmap's
discovery log (`afni-core_ROADMAP.md`). Decisions made by the project owner also get an
`R` number in the roadmap's "Resolved" table; the **Owner** column below points to it.

## How to read the tables

| Kind | Meaning |
|------|---------|
| **DEVIATION** | afni-core deliberately does something different from AFNI/SUMA. Results can differ. |
| **COPIED** | AFNI/SUMA does something odd and afni-core reproduces it on purpose, so numbers match. |
| **UNREPRODUCIBLE** | AFNI/SUMA's result depends on platform, compiler or rounding noise, so no implementation can match it exactly. afni-core picks a deterministic rule. |
| **NOT PORTED** | An AFNI/SUMA feature that afni-core does not implement. |
| **ADDED** | Something afni-core offers that AFNI/SUMA does not (it never changes a parity result). |

**Effect** says how large the difference is where we measured it. **Checked** says what
the comparison was made against: a live program, AFNI's own library called from a C
harness, AFNI's source read line by line, or nothing yet. Entries marked *unverified*
were not compared against a running AFNI.

**What "matches AFNI" means in the tests.** Reference values come from the live
programs (`cdf`, `3dFDR`, `ScaleToMap`, `SurfClust`, `3dClusterize`, `ROI2dataset`,
`SurfDist`, `ROIgrow`, `ConvertDset`, `MakeColorMap`), or, where no program exposes the
function, from AFNI's own C functions compiled into small harnesses (`display.c` colour
scales, `THD_bandpass_vectors`, the SUMA InstaCorr call sequence, FATCAT's `TrackIO.c`).
A harness proves agreement with the **source tree**, not that a shipped binary was
built from it. Every comparison uses the tolerance of the reference's own printing or
arithmetic (named in the test). Set `AFNI_CORE_LIVE=1` to re-run the live ones.

---

## 1. Across the whole crate

| ID | Kind | AFNI / SUMA | afni-core | Why | Owner |
|----|------|-------------|-----------|-----|-------|
| G-1 | DEVIATION | Mostly 32-bit `float` arithmetic. | `f64` for math, `f32` only for stored display/mesh buffers, and `f32` **inside the parity ports** (FDR curves, colour scales, overlay panes, colormaps, band indices) where the exact float path decides the result. | Accuracy, except where parity needs the float path. | design |
| G-2 | DEVIATION | Sentinel values and silent substitution: `99.99`/`999.99` statistics, `-1` distances, a non-positive time step silently set to 1.0. | Typed errors (`Error`). | "Do not silently invent values." | design |
| G-3 | DEVIATION | NaN/infinity are mostly ignored by the C comparisons (`<`/`>` are false for NaN) or cast to `int` (undefined behaviour). | Explicit policies (`NonFinitePolicy`, `MissingThreshold`, `missing_color`); invalid rows reported as `NaN` or errors. | Predictable, documented behaviour. | design |
| G-4 | DEVIATION | A tail (one- or two-sided) is often implied by the program or the statistic code. | `Tail` is always explicit; never inferred from a code. | `cdf -t2p` alone is two-sided for some codes and upper-tail for others (see S-1). | design |
| G-5 | DEVIATION | A surface's node count is guessed from file contents in places. | Never guessed: sparse datasets need a caller-supplied `node_count`. | AFNI "dense" files still carry an index list. | R1 |
| G-6 | DEVIATION | Volume sub-bricks keep their integer datum type. | Sub-bricks become `Float32` columns (scaled values); label table and `Label` role are kept. | One column model. Revisit if exact integer storage is needed. | design |
| G-7 | DEVIATION | `COLMS_RANGE` is stored and recomputed by AFNI. | A recorded range is read-time provenance only (`RangeReport`); never trusted. | Ranges go stale. | design |
| G-8 | DEVIATION | Plain `==` on data with NaN. | IEEE `==` plus a separate `eq_nan_aware`. | Datasets with missing samples compared reasonably. | R2 |
| G-9 | ADDED | n/a | `Dataset` preserves an index list of `0..n-1` as indexed, never collapsing it to dense. | Round-trip fidelity. | R1 |

## 2. Statistics (`stats.rs`, `special.rs`, `stat.rs`)

| ID | Kind | AFNI / SUMA | afni-core | Effect / Checked | Owner |
|----|------|-------------|-----------|------------------|-------|
| S-1 | DEVIATION | `cdf -t2p` is **two-sided** for `fizt`, `fitt`, `fico` but **upper-tail** for `fift`/`fict`. | Explicit `Tail {Lower, Upper, TwoSided}`; `TwoSided` only for symmetric distributions. | Calls differ from `cdf` unless the same tail is asked for. Checked live (`cdf`, `nifticdf`). | G-4 |
| S-2 | DEVIATION (naming) | Gamma "scale" parameter. | Called `rate`; the math is identical (`x * scale` goes into the incomplete gamma). | None. `nifti_stats -q 2 GAMMA 3 2` equals `Q(3, 4)`. | — |
| S-3 | DEVIATION | Binomial/Poisson tails at **non-integer** statistics use a continuous extension (`P(X > 2.5)` differs from `P(X > 2)`). | True step function (`floor(x)`); integer quantile (smallest `k` with `P(X > k) <= p`). | Differs only at non-integer inputs. Conformance compares integer statistics. | — |
| S-4 | DEVIATION | `correl_p2t` returns 0 when `nort < 1` although `correl_t2p` accepts it. | Accepts `nort >= 0` both ways. | Removes an asymmetry. | — |
| S-5 | COPIED | `p2dsetstat` on a correlation with `nfit > 1` is the **upper tail of R** (a multiple correlation, `R >= 0`), `-1sided` doubles p first. | Same: `nfit > 1` is one-sided. | `Correl(30,2,1)`: p=0.01 gives 0.537614, p=0.02 gives 0.501569. Checked live. | — |
| S-6 | DEVIATION | NIfTI correlation intents (1 parameter, `dof`) and AFNI correlation (`samples, nfit, nort`) share one code and are told apart by nothing explicit. | Classified **by structure** (AFNI always has `nfit >= 1`; the standard leaves `p2 = p3 = 0`); standard `Correl(dof)` is stored as `Correl(dof+1, 1, 0)`; malformed headers are **errors**; an explicit writer (`IntentOrigin`) may be stated. The old `Option`-returning APIs are `#[deprecated]` but kept. | Verified against `3dAFNItoNIFTI -pure` and the GIfTI writer. | R5, R6 |
| S-7 | DEVIATION | CDFLIB noncentral chi-square and F (`cumchn`, `cumfnc`) are only accurate to about 1e-5 absolute (`F(3,20;2)` at 3: AFNI 0.8237480, exact 0.8237555). | Implemented exactly (two independent implementations agree). | AFNI differs from truth by up to ~1e-5; tests use 2e-5 for these two. | R4 |
| S-8 | DEVIATION | Noncentral t: AS 243's tail cancels in the far tail on the wrong side of the shift. | Detects digit loss and integrates the definition instead. | Agrees with AFNI to ~2e-10 forward, ~2e-9 inverse where AFNI is accurate. | R4 |
| S-9 | DEVIATION | CDFLIB inverses are good to ~1e-8 absolute (`-1 0.5 CORREL 18` gives -5.0e-9, exact 0). | Exact to double precision. | Inverse tests use 5e-8 absolute. | — |
| S-10 | DEVIATION | Uniform reports `q = 1 - u` outside its support (e.g. 1.25). | Clamped probability. | Out-of-support inputs only. | — |
| S-11 | COPIED | Inverse Gaussian upper tail is `1 - cdf` (accuracy lost below ~1e-12). | Same; every other distribution keeps both tails accurate in log space. | Below 1e-12 only. | — |
| S-12 | DEVIATION | `student_t2p` returns p = 1 for `dof < 1`. | Evaluates it. | `dof < 1` only. | — |
| S-13 | DEVIATION | NIfTI `LOGPVAL`/`LOG10PVAL` read in various ways. | Interpreted as `exp(-abs(v))` / `10^-abs(v)` (the reference library writes positive `-log p`). | Sign convention only. | — |

## 3. FDR / MDF curves (`fdr.rs`, `curve.rs`)

| ID | Kind | AFNI / SUMA | afni-core | Effect / Checked | Owner |
|----|------|-------------|-----------|------------------|-------|
| F-1 | COPIED | Stored-curve algorithms are single precision. | `f32` arrays; p floored at `1e-15` and rounded to `f32`; `qsmal`/`m1`/`qfac` logic; **`0.1666667` (not 1/6)** in the cubic; AFNI's clamped four-point cubic (not Catmull-Rom). | `3dFDR` z-scores match with a worst difference of **exactly 0**. Checked live. | — |
| F-2 | COPIED | Ties in the step-up procedure are ordered by AFNI's unstable quicksort (`cs_sort_ff.c`), which decides the curve's `x0`. | Ported line for line (`afni_qsort_float_float`) and the **default** (`deterministic_ties: false`). A deterministic order is an option. | A stable order gave `x0 = 0.000391` vs AFNI's `0.000473`. | R7 |
| F-3 | COPIED | `interp_floatvec` returns the first sample for every `x` when the curve has two samples. | Copied. | Two-sample curves only. | — |
| F-4 | COPIED | `estimate_m1` bins `(int)((p-0.15)*20)`, truncating toward zero, so p in (0.10, 0.15) lands in bin 0. | Copied. | Affects the m1 estimate slightly. | — |
| F-5 | COPIED | `qsmal` uses 0 as "unset", so a first small q at rank 0 is lost. | Copied. | Harmless in practice. | — |
| F-6 | DEVIATION | `mri_fdr_curve` divides by `t2 - t1` (NaN for a degenerate curve); `interp_inverse_floatvec` divides by `yp - ym` on a flat segment. | Guarded (fraction 0 / left end). | Degenerate inputs only. | — |
| F-7 | DEVIATION | `fdrval -qinput` clamps only `q <= 0` to `1e-9` and `q >= 1` to `0.99999`. | `threshold_for_q` takes `q` in `(0, 1]`; `q = 0` is an error, `q = 1` gives threshold 0. | Edge q only. | — |
| F-8 | DEVIATION | For a q below the curve's range AFNI uses the **stored** brick maximum (`DSET_BSTAT_MAXABS`). | Computes max abs(stat) **from the data**. | Differs only if a file's stored statistic is stale. | R8 |
| F-9 | DEVIATION | FDR curves for the classic statistics 2-10 only, with an implicit rule. | Same default; any other statistic needs an explicit `FdrOptions::tail`. | No inherited undefined rule. | — |
| F-10 | NOT PORTED | `3dFDR -old` (count p = 1 voxels, skip the m0 correction). | Not implemented. | Legacy mode. | R9 |
| F-11 | DEVIATION | `student_t2p` etc. return sentinel statistics (`99.99`, `999.99`) for tiny p. | Used only to build test data; the sentinels exercise the p floor. | — | — |

## 4. Colors (`color.rs`, `afni_colors.rs`, `suma_colormaps.rs`, `labels.rs`)

| ID | Kind | AFNI / SUMA | afni-core | Effect / Checked | Owner |
|----|------|-------------|-----------|------------------|-------|
| C-1 | UNREPRODUCIBLE / DEVIATION | `display.c` colour scales are floating point; with fused multiply-add (default on Apple silicon) 6 entries across 45 scale/size combinations differ by **one byte** at the 0/360 hue wrap (e.g. `Reds_and_Blues` at 64 entries, index 31: `[255,0,11]` vs `[255,11,0]`). Only at sizes other than 256 (and 128). | Strict IEEE (what Rust computes; `-ffp-contract=off`). | At most one byte in 255, invisible on screen. Golden file built with strict IEEE. | R15 |
| C-2 | COPIED | AFNI's `pow` is `exp(y*log(x))`, 0 for `x <= 0`; channels `(int)(255*pow + 0.5)`. | Ported exactly. | The math library's `pow` can differ by a byte. | — |
| C-3 | COPIED | `Reds_and_Blues` subtracts `NBIG_MTOP + 1` in its upper half; `Reds_and_Blues_w_Green` paints only two entries green; `DC_spectrum_AJJ` uses `s` where its mirror uses `sb`. | Copied (sumaru's wider green band is a viewer choice and is not in core). | — | — |
| C-4 | DEVIATION | n/a | Tables are kept in AFNI order (index 0 = top = highest value) for byte comparison; `to_color_map` flips so position 0 is the lowest. | Orientation only. | — |
| C-5 | NOT PORTED | n/a | Golden tables were produced by compiling the exact `display.c` functions, not by running the binary (needs X11/Motif). | Proves agreement with the source, not the shipped binary. | — |
| C-6 | COPIED | SUMA's nine standard maps are built by fiducial interpolation in `float`. | `suma_colormaps` does the same in `f32`. | Checked against `MakeColorMap -std` (two printed decimals, so within 0.005). `ngray20` equals `gray20` in SUMA's source. | R11 |
| C-7 | NOT PORTED | AFNI's other `pbardefs.h` scales (`amber_circle`, gray circles, ...); non-linear colormaps with per-colour fractions (`SUMA_Linearize_Color_Map`). | Not implemented. | Unsupported maps give an explicit error. | R20 |
| C-8 | DEVIATION | Not AFNI definitions: sumaru's `fire`, `afni_p2_spanned`, `blue_white_red`. | Not ported. `amber_monochrome` matches `pbardefs.h`. The label-colour fallback palette (10 colours, key 0 gray) is a sumaru design, kept and named as such; label keys are `i64`; unlabeled keys transparent; duplicate keys are an error. | — | R11 |
| C-9 | DEVIATION | Interpolation space. | Default encoded RGB like AFNI/SUMA; linear-light is an option. At a shared stop position the earlier stop wins. | — | — |

## 5. Thresholds and overlays (`threshold.rs`, `overlay.rs`, `composite.rs`)

| ID | Kind | AFNI / SUMA | afni-core | Effect / Checked | Owner |
|----|------|-------------|-----------|------------------|-------|
| T-1 | COPIED | SUMA's modes are inconsistent at their boundary: `ABS_LESS_THAN` hides `-t < v < t` strictly (so `|v| = t` **passes**), `OUTSIDE_RANGE` hides `[lo, hi]` inclusively (ends **fail**). `LESS_THAN` is really "show `v >= t`". | Modelled honestly: `AbsoluteAbove` (`|v| >= t`), `Outside` (ends fail), `Above`, `Between`. `Below` exists for sumaru. | Boundary tests; confirmed independently by `SurfClust` (`-athresh`, `-ex_range`). sumaru's `Outside(-T,T)` is really `AbsoluteAbove(T)`. | R16 |
| T-2 | DEVIATION | A sample whose THRESHOLD value is NaN is never hidden and is drawn at full strength (`<`/`>` are false for NaN). | **Hidden** by default; `MissingThreshold::Show` reproduces AFNI/SUMA. | NaN threshold values only. | R13 |
| T-3 | DEVIATION | A NaN intensity is cast to `int` in AFNI (undefined behaviour). | Configurable `missing_color` (transparent by default). | NaN intensities only. | — |
| T-4 | COPIED | SUMA: colour `i` of `N` at `i/N`, last colour held from `(N-1)/N` (pane rule); AFNI volumes index panes from the **top**. | `PaneRule::{SumaInterpolated, SumaBanded, AfniBanded, Direct}` in `f32`. A continuous map with `N` stops at `i/(N-1)` is a separate type. | Checked against `ScaleToMap` (54 cases, exact). AFNI lands one pane lower exactly on a pane boundary. | R17 |
| T-5 | COPIED | AFNI fade: `255*((1-floor)*|v|/t)^k + 255*floor` as a byte, `rintf`, **clamped to 222**, `v == 0` rejected. SUMA fade: `min(1, |v|/t)`, no cap, no floor, only `ThreshRange[0]`. | Both as `FadeModel::{Afni, Suma}`. sumaru's default (0.85/0.5/0.35) and cubic/quartic curves are a viewer style (`FadeStyle`), not AFNI. | `FadeStyle::max_alpha` caps any model, so with `Afni` the cap is `min(222/255, cap)`. | R18 |
| T-6 | DEVIATION (*unverified*) | For `Above(t)` AFNI neither fades nor rejects a NEGATIVE threshold value, so it can draw opaque. | Negative values are transparent. | Could not be checked without the AFNI GUI. Revisit in Phase 11. | R19 |
| T-7 | DEVIATION | Brightness above 1 or below 0 can leave 0..1; SUMA leaves it to the GL clamp. | Clamps immediately. | Display only. | — |
| T-8 | NOT PORTED | `ScaleToMap -anr`, `-perc_clp`, `top_frac`, "no color" rows, coordinate-bias and contour options; `-br` only `(0, 1]` on the CLI. | The library allows `(0, 2]`; the rest are not implemented. | Covered cases pass exactly. | R20 |

## 6. Surfaces: topology, geometry, clustering (`topology.rs`, `mesh.rs`, `cluster.rs`)

| ID | Kind | AFNI / SUMA | afni-core | Effect / Checked | Owner |
|----|------|-------------|-----------|------------------|-------|
| M-1 | DEVIATION | `SurfClust -rmm r` (and `ROIgrow -lim`) build layered neighbourhoods whose distances are not shortest paths (a squared-length bug in `SUMA_getoffsets2`'s precursor choice); SUMA reaches **fewer** nodes than true graph distance on irregular meshes. | True Dijkstra distance. | At 14 mm on the irregular 642-node sphere SurfClust finds **35** clusters, afni-core **32** (never more). `ROIgrow` missed 6, 11 and 1 nodes in three cases. Identical on the regular sphere and at small radii. Pinned by tests. | R21 |
| M-2 | COPIED | Edge-ring neighbourhoods (`-rmm -N`). | Same. | Every ring case matches (counts, areas, means, extremes, variance, SE, centroid, centre of mass). | — |
| M-3 | UNREPRODUCIBLE | On regular meshes many clusters have equal areas whose `float` sums differ in the last bit; that noise orders them. | Stable sort in discovery order (`f64` sums). | Regular-mesh cases compared as sets; the irregular mesh matches rank for rank. | — |
| M-4 | COPIED | A node whose value is exactly 0 is inactive; clusters start from the **highest** node index downward; negative `-amm2` without `-n` means "at least that many nodes". | `exclude_zero_values` (off by default, on in the conformance cases); same seed order (`seed_node` recorded); mapped. | — | R23 |
| M-5 | DEVIATION | SurfClust merges touching positive and negative regions. | Default `Tails::Merged` (SurfClust's). `Tails::Separate` is the sumaru / `3dClusterize -bisided` behaviour; callers choose it explicitly. | — | R22 |
| M-6 | COPIED / ADDED | Centre of mass `sum(v xyz)/sum(v)` with **signed** values (absurd or NaN when values cancel). | Reproduced; `center_of_mass_abs` (|v|-weighted) added alongside. | — | R14 |
| M-7 | COPIED | Node normal is the normalized sum of **unit** triangle normals; node area is one third of adjacent triangle areas. | Same (distinguished from area- and angle-weighted). | `SurfaceMetrics`/`SurfMeasures`, ~1e-6 relative (SUMA uses `float`). | — |
| M-8 | DEVIATION | n/a | Topology defects (degenerate, duplicate, non-manifold, bow-tie, inconsistent winding) are **reported**, not fatal; only an empty node set or bad index is an error (`require_clean` is strict). | Lets a viewer open a damaged mesh. | — |
| M-9 | DEVIATION | n/a | `TopologyId` is an FNV-1a hash in triangle order; reordering triangles changes it. | Cheap. | R25 |
| M-10 | NOT PORTED | SurfClust central-node and weighted-central-node columns; file-writing options; CLI string building. | Not implemented (consumer / I/O business). | `-thresh_col` is supported by the API but not exercised live. | R24 |
| M-11 | DEVIATION | Extremes ties go to the first-discovered node. | Lowest node number. | Test data avoid exact ties. | — |
| M-12 | DEVIATION | SurfClust applies both size limits. | Same (`min_area` and `min_nodes` both). | sumaru used only one metric. | — |

## 7. Volume clustering (`volume_cluster.rs`)

| ID | Kind | AFNI / SUMA | afni-core | Effect / Checked | Owner |
|----|------|-------------|-----------|------------------|-------|
| V-1 | DEVIATION | `3dClusterize -clust_vol V` is read as **at least `V` voxels** (the value is negated and truncated like `-clust_nvox`); the report's "Volume threshold" line multiplies by voxel volume afterwards, hiding it. | `min_volume` is a real volume (count x voxel determinant); use `min_voxels` for parity. | `-clust_vol 10` dropped a 3-voxel, 45 µl cluster. Read from the source and confirmed live. | — |
| V-2 | DEVIATION | The report column labelled "Volume" prints the **voxel count**. | `voxel_count` and `volume` are separate. | Labelling only. | — |
| V-3 | COPIED | Clusters ranked by a stable bubble sort on voxel count; for `-bisided` the right tail is clustered first, then the left, then the lot is sorted. | Same (stable, right tail before left). | Cluster maps identical in 18 cases. | — |
| V-4 | DEVIATION | With 3333 or more clusters the sort is skipped (`** TOO MANY CLUSTERS TO SORT **`). | Always ranked. | Rare. | — |
| V-5 | COPIED | A voxel whose DATA value is exactly zero is never in a cluster. | `exclude_zero_data` (default on). | — | — |
| V-6 | COPIED | Thresholds are inclusive: `v >= t`, `v <= t`, `v <= left || v >= right`, `lo <= v <= hi`. | Same (note this is *not* core's `Threshold::Outside`). | — | — |
| V-7 | COPIED | Neighbours by voxel index only, even for oblique data; the report uses the **cardinal** grid for coordinates, not the oblique matrix. | Same; pass the cardinal matrix for byte parity, the real one for true positions. | Oblique case equals its cardinal twin. | — |
| V-8 | DEVIATION | Thresholds are 32-bit floats. | `f64`; use `t as f32 as f64` for exact agreement on 32-bit data (the test does). | Boundary voxels only. | — |
| V-9 | COPIED / DEVIATION | Mean signed; SEM from a one-pass sum of squares; peak = largest |value| (first met on a tie); CM weights by |value|; coordinates printed with one decimal. | Same definitions; SEM by the two-pass formula; neighbours visited in AFNI's k,j,i order so ties agree. | Report columns within their printed precision; maps exact. | — |
| V-10 | NOT PORTED | `p=` thresholds, `-abs_table_data`, the global totals line, `-binary`, `-pref_dat`, `-orient`, `-mask_from_hdr`. | Not implemented (callers derive them from the labels, `stats`, a matrix). | — | — |

## 8. ROIs (`roi.rs`, `roi_ops.rs`, `roi_edit.rs`)

| ID | Kind | AFNI / SUMA | afni-core | Effect / Checked | Owner |
|----|------|-------------|-----------|------------------|-------|
| O-1 | COPIED | `ROI2dataset -nodelist` writes every node as drawn (junction repeats kept); `-nodelist.nodups` keeps first occurrences; SUMA's dataset path (`SUMA_NodesInROI`) drops a node repeating the last node of the previous stroke. | `drawn_nodes`, `drawn_nodes_unique`, `ordered_nodes` respectively (three names, three behaviours). | 16 real files reproduced exactly. | — |
| O-2 | UNREPRODUCIBLE | When ROIs with **different labels** share a node, `ROI2dataset` keeps one by sorting with the C library's `qsort`, so the winner is platform luck (about half of 90 contested nodes went each way). | `OverlapPolicy::FirstWins` by default (`LastWins`, `Error` available). | The test checks the node set, every uncontested label, and that a contested label is one of its claimants. | — |
| O-3 | COPIED | `-pad_to_node N` writes rows 0..=N and refuses `N` below the largest node; `-pad_label` defaults to 0. | Same. | Matches. | — |
| O-4 | DEVIATION | `ROIgrow` under-reaches (see M-1). | `dilate_by_distance` is the true distance. | See M-1. | R21 |
| O-5 | DEVIATION (*unverified*) | SUMA's fill (`SUMA_FillToMask`) floods from a seed up to the mask; behaviour with a boundary that leaks is not known. | `fill_enclosed` floods from a seed; the editor **refuses** a fill that reaches the rim of an open surface. | Not compared on real data. | — |
| O-6 | DEVIATION | SUMA edits through an action stack in the viewer. | Plain `RoiCommand` values with inverses and an `RoiEditor`. | Replay of strokes was not needed. | — |
| O-7 | NOT PORTED | n/a | `.niml.roi` triangle paths and per-stroke distances: present on the core type, empty after reading a file (the file does not carry them). | — | — |

## 9. Time series and seed correlation (`signal.rs`, `instacorr.rs`)

| ID | Kind | AFNI / SUMA | afni-core | Effect / Checked | Owner |
|----|------|-------------|-----------|------------------|-------|
| P-1 | COPIED | `THD_bandpass_vectors`: even FFT length, band bins in `float` with ties to even, 0.5 (or 0.05 with orts) edge taper, mean and Nyquist dropped, orts filtered then projected out, removed-dimension count scaled by `len/nfft`. | Ported step by step. | 21 cases: removed count **exact**, series within ~2e-6. AFNI called directly through libmri. | — |
| P-2 | COPIED | This AFNI build pads to the next EVEN length (`USE_FFTN`: 137 → 138, 211 → 212). | Same. A build without it pads to a 2-3-5 smooth size and would differ. | — | — |
| P-3 | UNREPRODUCIBLE / DEVIATION | With a linear detrend and Legendre orts (SUMA's setup) the filter reduces the constant and linear orts to rounding noise; AFNI's pseudo-inverse keeps those noise directions and projects out **two arbitrary extra dimensions** that depend on 32-bit rounding. | Drops negligible columns. The removed-dimension **count** (which includes those orts) is still AFNI's. | Correlations differ by 0.014 (SUMA default, 120 samples) to 0.035 (polort 1, 100 samples). Cases without such orts agree to 4e-7. | — |
| P-4 | DEVIATION | A non-positive `dt` is set to 1.0; `ftop <= fbot` assigns `fbot` (typo: becomes a low-pass). | Both are errors. | — | — |
| P-5 | COPIED / DEVIATION | `THD_bandpass_OK` and `remain_dim` treat a 1-2 bin band as "ignore the filter", `THD_bandpass_vectors` filters it. The count is not capped (102 removed for 100 samples). | The filter follows `THD_bandpass_vectors`; `bandpass_remaining_dimension` follows the other. The `Correl` statistic is attached only if `samples - 1 - removed > 0`. | sumaru capped at `samples - 1`. | — |
| P-6 | COPIED | Filtering needs at least 9 samples even with no band; `polort + 1 < samples - 3`. A mean detrend adds 0 to the count, linear 1, quadratic 2. | Same. | sumaru enforced 9 only with a band and added the polynomial count after scaling. | — |
| P-7 | COPIED | Dataset rows get a LINEAR detrend; an outside seed gets a MEAN detrend; an ROI seed is the mean of prepared rows scaled to unit length. | Same. | — | — |
| P-8 | DEVIATION | AFNI scales an outside seed to unit length even with `normalize_dset` off. | Raw mode stays raw (plain dot product, **no** statistic attached). | Raw mode only. | — |
| P-9 | DEVIATION | A row with no variance is left unscaled (`THD_normalize` skips squared length <= 1e-20), giving ~0. | Marked invalid, `NaN`. A non-finite row is zeroed before filtering and is `NaN` (two series share one FFT, so a NaN would otherwise leak). | Degenerate rows only. | — |
| P-10 | ADDED | AFNI attaches the statistic downstream. | Output carries `Correl(samples, 1, removed_dof)`. | — | — |
| P-11 | NOT PORTED | Volume-side `THD_instacorr` (blur, despike, `qdet = polort`, global orts); methods other than Pearson; Fisher-z. | Not implemented. | — | — |

## 10. Graphs and tracts (`graph.rs`, `tract.rs`, and `afni-io`)

| ID | Kind | AFNI / SUMA | afni-core | Effect / Checked | Owner |
|----|------|-------------|-----------|------------------|-------|
| N-1 | COPIED (sumaru DEVIATED) | Sparse `INDEX_LIST` names end nodes by node INDEX, not position. | Indices stored and resolved to positions; a missing node is an error. | sumaru treated them as positions (right only when indices are 0..n-1). Checked with indices 5..8. | — |
| N-2 | COPIED | Full matrices are column-major; triangular layouts pack `(1,0) (2,0) .. (n-1,0) (2,1) ..`. | Same; endpoints by binary search, cell lookup in closed form. | `tri_diag` has no real fixture (AFNI always reads that size as a bigger triangle first); it rests on the source and unit tests. | — |
| N-3 | COPIED | `-graph_XYZ_LPI` flips x and y to AFNI's frame. | Identical (`domain::flip_dicom_ras`). | Checked on all nodes. | — |
| N-4 | DEVIATION | Coordinates are in AFNI's DICOM frame. | Stored as in the file; `positions_ras` and `TractSet::flipped` convert. | sumaru converted at read time. | — |
| N-5 | DEVIATION | `matrix_size` of a sparse graph holds the EDGE count; layouts that imply their edges still write an empty `INDEX_LIST`; `COLMS_RANGE` goes stale. | `matrix_size` ignored on read, `n n` written; empty `INDEX_LIST` kept as an extra; ranges dropped on rewrite (AFNI recomputes). | `ConvertDset` accepts and re-saves the graphs we write. | — |
| N-6 | COPIED | `TAYLOR_TRACT_DATUM`'s second number counts VALUES (3 per point). | Same; a disagreeing count is an error. A wrong count that still divides by 3 cannot be detected (inherent to the format). | Real AFNI-written files, ASCII and binary. | — |
| N-7 | COPIED | `Tract_Length` is the polyline length. | Same (`f64`). | Equal to AFNI's `Tract_Length` within `f32` rounding. | — |
| N-8 | NOT PORTED | Per-point tract scalars (FATCAT files carry none), tract colours, resampling between surfaces. | Not implemented. | — | — |
| N-9 | *unverified* | No AFNI program reads a `.niml.tract` without diffusion data. | Tract **writing** is checked by round trips and by matching the structure AFNI wrote, not by AFNI reading our file. | — | — |

## 11. Project-level decisions that change processing

These are the choices **you** made, in order, with where they bite. (Full text in the
roadmap's "Resolved" table.)

| Ref | Decision | Changes results? |
|-----|----------|------------------|
| R1 | Keep a file's representation (index list stays indexed). | No. |
| R2 | NaN-aware equality alongside `==`. | No. |
| R3 | Malformed FDR/MDF curves: error by default, skip with warning on request. | Only in how bad files are handled. |
| R4 | Implement noncentral t, F, chi-square exactly. | Yes: differs from AFNI's CDFLIB for F and chi-square by up to ~1e-5 (S-7). |
| R5 | Correlation in files classified by structure; malformed is an error. | Only for malformed headers. |
| R6 | Deprecate (but keep) the lenient `Option` stat APIs. | No. |
| R7 | New FDR curves default to AFNI's exact tie order. | Parity: yes (F-2). |
| R8 | Max abs(stat) from the data, not stored statistics. | Only with stale stored statistics (F-8). |
| R9 | Skip `3dFDR -old`. | Feature absent (F-10). |
| R10 | Keep the `fdr` fixture dumps. | No. |
| R11 | Port SUMA's named colormaps; sumaru's three non-AFNI maps are sumaru's own. | No. |
| R12 | `afni-io` depends on `afni-core` through git (`main`). | No. |
| R13 | A NaN threshold value is hidden by default. | Yes vs AFNI/SUMA (T-2). |
| R14 | Add a |value|-weighted centre of mass. | Added field only. |
| R15 | Strict IEEE bytes are "AFNI's" colour scale. | At most one byte at non-default sizes (C-1). |
| R16 | Defer sumaru's `Outside` to `AbsoluteAbove` migration. | sumaru unchanged (T-1). |
| R17 | Defer sumaru's pane-rule adoption. | sumaru unchanged (T-4). |
| R18 | Defer sumaru's fade style change. | sumaru unchanged (T-5). |
| R19 | Leave negative one-sided fade values transparent until checked. | Possibly vs AFNI (T-6). |
| R20 | Defer non-linear colormaps. | Feature absent (C-7, T-8). |
| R21 | Keep true shortest-path radius; no SUMA-compatible mode. | Yes: fewer nodes in SUMA (M-1, O-4). |
| R22 | Two-sided clusters keep tails separate when asked; core default stays `Merged`. | Parity for SurfClust (M-5). |
| R23 | Keep SurfClust's seed and tie order. | Parity (M-4). |
| R24 | No central-node columns. | Feature absent (M-10). |
| R25 | `TopologyId` stays triangle-order sensitive. | No. |

Not yet decided (for Phase 11): whether to run the AFNI GUI check for T-6, and which
sumaru behaviours to migrate; both are listed in the roadmap.

## 12. Where the same differences also concern sumaru

`sumaru` is read-only here, so these are recorded, not changed. When it is migrated onto
`afni-core`, expect small, explainable differences in: FDR curves (`f64`, stable sort,
no MDF, Catmull-Rom), thresholds (`Outside` is `AbsoluteAbove`), colour panes
(`i/(N-1)` stops vs SUMA's `i/N`), fades (a viewer style), clustering (bisided default,
one size metric, lowest-index seeds), sparse graph edge ends (positions, not indices),
InstaCorr (min-sample rule, DOF scaling, cap), ROI node orders, and the `fire` /
`afni_p2_spanned` palettes (not AFNI definitions).
