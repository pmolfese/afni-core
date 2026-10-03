#!/usr/bin/env bash
# Regenerate the committed AFNI conformance fixtures in tests/data/conformance/.
#
# Normal `cargo test` never runs this: the generated files are committed so the
# tests work on machines without AFNI. Run this script only when you add cases
# or move to a new AFNI release, then review and commit the diff.
#
#   tests/data/regenerate_conformance.sh
#
# Each fixture starts with `#` provenance lines (AFNI version, date, generator)
# and then one case per line:
#
#     <cdf arguments> => <expected value>
#
# e.g. `-t2p fitt 2.0 10 => 0.073388` means `cdf -t2p fitt 2.0 10` printed
# `p = 0.073388`. The same lines are replayed against live AFNI by the
# `AFNI_CORE_LIVE=1` test in tests/conformance_harness.rs.
#
# NOTE: `cdf` prints only 6 significant digits, so these values are accurate to
# roughly 1e-6 relative. Tests must use a tolerance no tighter than that.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
out="$here/conformance/cdf.ref"

command -v cdf >/dev/null || { echo "cdf not found on PATH (is AFNI installed?)" >&2; exit 1; }

# `afni -ver` prints e.g. "Precompiled binary macos_13_ARM: ... (Version AFNI_26.2.08 'Gordian III')".
version="$(afni -ver 2>&1 | head -n 1)"

# Cases: arguments to `cdf`, one per line.
#
# Tail convention (observed with AFNI_26.2.08, recorded in the roadmap's
# discovery log): `cdf -t2p` is TWO-SIDED for the symmetric statistics (fizt,
# fitt, fico: z=0 gives p=1, z=1 gives 0.3173) but UPPER-TAIL for fift/fict.
# Phase 2 must make the tail explicit instead of inheriting this.
cases=(
  "-t2p fizt 0.0"
  "-t2p fizt 1.0"
  "-t2p fizt 1.96"
  "-t2p fizt 3.0"
  "-t2p fitt 0.0 10"
  "-t2p fitt 1.0 10"
  "-t2p fitt 2.0 10"
  "-t2p fitt 3.0 10"
  "-t2p fift 4.0 3 20"
  "-t2p fict 3.84 1"
  "-t2p fict 10.0 5"
  "-p2t fizt 0.025"
  "-p2t fitt 0.05 10"
  "-p2t fift 0.05 3 20"
  "-p2t fict 0.05 1"
)

{
  echo "# AFNI conformance fixture for the 'cdf' program (probability <-> statistic)."
  echo "# afni_version: $version"
  echo "# generated: $(date +%Y-%m-%d)"
  echo "# generator: tests/data/regenerate_conformance.sh"
  echo "# precision: cdf prints 6 significant digits"
  for c in "${cases[@]}"; do
    # shellcheck disable=SC2086  # intentional word splitting of $c
    result="$(cdf $c | head -n 1)"      # e.g. "p = 0.073388"
    echo "$c => ${result##*= }"
  done
} > "$out"

echo "wrote $out"
