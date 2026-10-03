#!/usr/bin/env bash
# Regenerate tests/data/conformance/scaletomap.ref: SUMA's own value-to-color
# mapping, run through AFNI's `ScaleToMap` program (SUMA_ScaleToMap in
# SUMA_Color.c). Used to check afni-core's overlay evaluation against SUMA for
# interpolated/banded/direct mapping, intensity clipping, masking, and the
# brightness factor.
#
#   tests/data/regenerate_scaletomap.sh        (needs AFNI on PATH)
#
# Line format (one case per line, fields separated by ' | '):
#   <name> | <ScaleToMap options> | <colormap: r,g,b r,g,b ... lowest value first>
#          | <input values> => <node:r,g,b ...>
# Nodes that ScaleToMap masked are ABSENT from the output (the cases use
# -nomsk_col); `ScaleToMap` masks exact zeros by default. Colormaps are passed with
# -frf so the first row is the lowest value. Output colors have 6 decimals.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
out="$here/conformance/scaletomap.ref"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
command -v ScaleToMap >/dev/null || { echo "ScaleToMap not found (is AFNI installed?)" >&2; exit 1; }
export AFNI_DONT_LOGFILE=YES

# Colormaps (rows "r,g,b"), lowest value first.
M4="1,0,0 0,1,0 0,0,1 1,1,0"
M5="1,0,0 0,1,0 0,0,1 1,1,0 0,1,1"
# A 16-color ramp made by awk so test and fixture agree: (i/15, 1-i/15, 0.25 or 0.75).
RAMP16="$(awk 'BEGIN{for(i=0;i<16;i++) printf "%s%.6f,%.6f,%.2f", (i?" ":""), i/15, 1-i/15, (i%2)?0.75:0.25}')"
# Eleven grays, each a different brightness.
GRAY11="$(awk 'BEGIN{for(i=0;i<11;i++) printf "%s%.6f,%.6f,%.6f", (i?" ":""), i/10, i/10, i/10}')"

U="0.001 0.1 0.2 0.25 0.4 0.5 0.6 0.75 0.8 0.999 1"
S="-1 -0.75 -0.5 -0.25 0.25 0.5 0.75 1"
D="0 1 2 3 4.2 4.7 -3 100"
W="0.01 0.15 0.3 0.45 0.62 0.77 0.93 0.4999999 0.5000001"

run_case() { # name opts map values
  local name="$1" opts="$2" map="$3" values="$4"
  local mapf="$work/map.1D" valf="$work/val.1D" raw colors
  echo "$map" | tr ' ' '\n' | tr ',' ' ' > "$mapf"
  echo "$values" | tr ' ' '\n' > "$valf"
  # `|| true`: a rejected option set must give an empty case, not abort the script.
  # shellcheck disable=SC2086  # intentional word splitting of the option string
  raw="$(ScaleToMap -input "$valf" -1 0 -cmapfile "$mapf" -frf $opts -nomsk_col 2>/dev/null || true)"
  colors="$(printf '%s\n' "$raw" \
    | awk '/^[0-9]+ [0-9.]+ [0-9.]+ [0-9.]+$/ {printf "%s%d:%s,%s,%s", (n++?" ":""), $1, $2, $3, $4}')"
  echo "$name | $opts | $map | $values => $colors"
}

{
  echo "# SUMA value-to-color mapping, from AFNI's ScaleToMap program"
  echo "# afni_version: $(afni -ver 2>&1 | head -n 1)"
  echo "# generated: $(date +%Y-%m-%d)"
  echo "# generator: tests/data/regenerate_scaletomap.sh"
  for mapname in M4 M5 RAMP16 GRAY11; do
    map="${!mapname}"
    run_case "interp_unit_$mapname"   "-interp -clp 0 1"       "$map" "$U"
    run_case "nointerp_unit_$mapname" "-nointerp -clp 0 1"     "$map" "$U"
    run_case "interp_signed_$mapname" "-interp -clp -1 1"      "$map" "$S"
    run_case "nointerp_signed_$mapname" "-nointerp -clp -1 1"  "$map" "$S"
    run_case "clip_narrow_$mapname"   "-interp -clp 0.2 0.8"   "$map" "$U"
    run_case "apr_$mapname"           "-interp -apr 0.8"       "$map" "$U"
    run_case "autorange_$mapname"     "-interp"                "$map" "$W"
    run_case "autorange_banded_$mapname" "-nointerp"           "$map" "$W"
    run_case "mask_range_$mapname"    "-interp -clp 0 1 -msk 0.4 0.6" "$map" "$U"
    run_case "bright_half_$mapname"   "-interp -clp 0 1 -br 0.5"  "$map" "$U"
    run_case "bright_quarter_$mapname" "-nointerp -clp 0 1 -br 0.25" "$map" "$U"
    run_case "combo_$mapname"         "-nointerp -clp 0.1 0.9 -br 0.75 -msk 0.45 0.55" "$map" "$U"
    run_case "direct_$mapname"        "-direct"                "$map" "$D"
  done
  # A value range of zero width: SUMA uses the middle color.
  run_case "flat_range_M5" "-interp" "$M5" "0.5 0.5 0.5"
  run_case "flat_range_M4" "-interp" "$M4" "0.5 0.5 0.5"
} > "$out"
echo "wrote $out ($(grep -vc '^#' "$out") cases)"
