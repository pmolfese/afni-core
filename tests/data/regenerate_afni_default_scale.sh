#!/usr/bin/env bash
# Regenerate tests/data/conformance/afni_default_scale.ref: AFNI's DEFAULT overlay
# color scale, read from the AFNI source tree.
#
#   AFNI_SRC=~/Documents/Programming/afni tests/data/regenerate_afni_default_scale.sh
#
# AFNI starts with the scale named by AFNI_COLORSCALE_DEFAULT, which afni.c sets to
# `Reds_and_Blues_Inv` (24 May 2019). That scale is not computed by display.c like
# the nine "big" maps: pbardefs.h lists its 256 colors as hex strings. This script
# records the name from afni.c and the 256 colors from pbardefs.h, verbatim.
#
# Lines:
#   default <name>                        from afni.c
#   scale   r,g,b r,g,b ...  (256 entries, AFNI index order: index 0 = top of the bar)
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
src="${AFNI_SRC:-$HOME/Documents/Programming/afni}/src"
out="${REGEN_OUT:-$here/conformance/afni_default_scale.ref}"
[ -f "$src/afni.c" ] && [ -f "$src/pbardefs.h" ] || { echo "AFNI source not found at $src (set AFNI_SRC)" >&2; exit 1; }

name="$(sed -n 's/.*PUTENV("AFNI_COLORSCALE_DEFAULT","\([^"]*\)").*/\1/p' "$src/afni.c" | head -1)"
[ -n "$name" ] || { echo "AFNI_COLORSCALE_DEFAULT not set in afni.c" >&2; exit 1; }

NAME="$name" SRC="$src" python3 - > "$out" <<'PY'
import os, re
name, src = os.environ["NAME"], os.environ["SRC"]
text = open(os.path.join(src, "pbardefs.h")).read()
# The colorscale definition is a C string array that starts with its name.
start = text.index('"%s "' % name)
end = text.index("};", start)
colors = re.findall(r"#([0-9a-fA-F]{6})", text[start:end])
print("# AFNI's default overlay colorscale (regenerate with regenerate_afni_default_scale.sh)")
print("default", name)
print("scale", " ".join("%d,%d,%d" % tuple(int(h[i:i+2], 16) for i in (0, 2, 4)) for h in colors))
PY
echo "wrote $out"
