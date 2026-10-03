#!/usr/bin/env bash
# Regenerate tests/data/conformance/roi_ops.ref: SUMA's own surface distances and
# ROI growth, used by tests/roi_conformance.rs to check afni_core::roi_ops.
#
#   tests/data/regenerate_roi_refs.sh        (needs AFNI on PATH, and python3)
#
# Uses the meshes made by regenerate_surface_refs.sh (tests/data/surfaces/).
# Lines:
#   distance <mesh> <from> <to> <graph distance>     from `SurfDist` (two decimals)
#   grow <mesh> <lim> <seed,seed,..> => <node node ...>   from `ROIgrow -lim`
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
surf="$here/surfaces"
# REGEN_OUT lets a test write somewhere else (the live tests do, so they never
# overwrite the committed file that the other tests are reading).
out="${REGEN_OUT:-$here/conformance/roi_ops.ref}"
command -v SurfDist >/dev/null || { echo "SurfDist not found (is AFNI installed?)" >&2; exit 1; }
export AFNI_DONT_LOGFILE=YES
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

{
  echo "# SUMA graph distances and ROI growth (regenerate with regenerate_roi_refs.sh)"
  echo "# afni_version: $(afni -ver 2>/dev/null | head -1 || echo unknown)"
  python3 - "$work" <<'PY'
import random, sys
rng = random.Random(8)
work = sys.argv[1]
# node pairs per mesh (642 nodes for ico3 and irr3, 162 for irr2)
for mesh, n in (("ico3", 642), ("irr3", 642), ("irr2", 162)):
    with open(f"{work}/{mesh}.pairs", "w") as f:
        for _ in range(40):
            a, b = rng.randrange(n), rng.randrange(n)
            f.write(f"{a} {b}\n")
PY
  for mesh in ico3 irr3 irr2; do
    (cd "$work" && SurfDist -i "$surf/$mesh.asc" -input "$mesh.pairs" 2>/dev/null) \
      | awk -v m="$mesh" '!/^#/ && NF == 3 { print "distance", m, $1, $2, $3 }'
  done
  # ROI growth from a few seeds by several distances.
  for spec in "ico3|10|8" "ico3|10,300|15" "ico3|5,6,200|25" "irr3|10|8" "irr3|10,300|15" \
              "irr3|5,6,200|25" "irr3|77|40" "irr2|3|20" "irr2|3,100|30"; do
    IFS='|' read -r mesh seeds lim <<< "$spec"
    echo "$seeds" | tr ',' '\n' > "$work/seeds.1D"
    rm -f "$work"/g.1D
    (cd "$work" && ROIgrow -i "$surf/$mesh.asc" -roi_nodes seeds.1D -lim "$lim" -prefix g >/dev/null 2>&1)
    nodes="$(grep -v '^#' "$work/g.1D" | awk 'NF { printf "%s%s", (n++ ? " " : ""), $1 }')"
    echo "grow $mesh $lim $seeds => $nodes"
  done
} > "$out"
echo "wrote $out"
