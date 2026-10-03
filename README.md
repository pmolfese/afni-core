# afni-core

File-neutral data models and algorithms for AFNI/SUMA tools: what decoded data
*means* (statistics, colors, thresholds, clusters, ROIs), shared by the
`sumaru` and `afniru` viewers and command-line tools. File parsing lives in the
separate [`afni-io`](../afni-io) crate; this crate must never depend on it.

**Status:** Phase 0 (foundations: error/numeric conventions, conformance
harness, dependency enforcement). See `afni-core_ROADMAP.md`.

## Layout

| Path | Contents |
|---|---|
| `src/lib.rs` | Crate docs, scope, module list |
| `src/error.rs` | Crate-wide `Error` / `Result` |
| `src/numeric.rs` | `f32`/`f64` rules, checked indices, NaN/Inf policy |
| `tests/common/` | Shared AFNI-comparison test harness |
| `tests/data/conformance/` | Committed AFNI reference values |
| `docs/ARCHITECTURE.md` | Dependency direction and design decisions |

## Testing

```sh
cargo test                   # committed fixtures only; AFNI not required
AFNI_CORE_LIVE=1 cargo test  # also replay cases against AFNI on PATH
tests/data/regenerate_conformance.sh   # refresh fixtures from AFNI
```

## License

Public domain. United States Government work (17 U.S.C. § 105); outside the US,
rights are waived under [CC0 1.0](LICENSE).
