#!/usr/bin/env bash
# Regenerate tests/data/conformance/suma_colormaps.ref: SUMA's standard color
# maps as printed by AFNI's own `MakeColorMap -std <name>` (which builds them with
# SUMA_MakeStandardMap in SUMA_Color.c). Used by tests/suma_colormaps_conformance.rs.
#
#   tests/data/regenerate_suma_colormaps.sh        (needs AFNI on PATH)
#
# Line format: <name> | r,g,b r,g,b ...    (entry 0 = lowest value; two decimals,
# the precision MakeColorMap prints)
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
out="$here/conformance/suma_colormaps.ref"
command -v MakeColorMap >/dev/null || { echo "MakeColorMap not found (is AFNI installed?)" >&2; exit 1; }
export AFNI_DONT_LOGFILE=YES

{
  echo "# SUMA standard color maps from AFNI's MakeColorMap -std"
  echo "# afni_version: $(afni -ver 2>/dev/null | head -1 || echo unknown)"
  for name in rgybr20 bgyr19 gray02 gray_i02 gray20 ngray20 bw20 byr64 bgyr64; do
    rows="$(MakeColorMap -std "$name" 2>/dev/null | awk 'NF==3 {printf "%s%s,%s,%s", (n++?" ":""), $1, $2, $3}')"
    echo "$name | $rows"
  done
} > "$out"
echo "wrote $out"
