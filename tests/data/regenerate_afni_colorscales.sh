#!/usr/bin/env bash
# Regenerate tests/data/conformance/afni_colorscales.ref: AFNI's built-in color
# scales computed by AFNI's OWN C code.
#
# The scales live in display.c (`NJ_bigmaps_init`, `DC_spectrum_AJJ`,
# `DC_spectrum_ZSS`), which needs X11/Motif to build as a whole and is not
# available as a command-line tool. So this script copies exactly those
# functions (plus the macros they use from display.h) out of the AFNI source
# tree, byte for byte, into a tiny stand-alone program, and prints the tables.
# Nothing is re-implemented: if AFNI changes these functions, rerunning this
# script shows it.
#
#   AFNI_SRC=~/Documents/Programming/afni tests/data/regenerate_afni_colorscales.sh
#
# Output lines:
#   scale  <npane_big> <name> : r,g,b r,g,b ...        (AFNI index order: index 0 = top)
#   spec   <AJJ|ZSS> <hue> <gamma> : r,g,b
# Normal `cargo test` does not run this; the output is committed.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
src="${AFNI_SRC:-$HOME/Documents/Programming/afni}/src"
out="$here/conformance/afni_colorscales.ref"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
[ -f "$src/display.c" ] || { echo "AFNI source not found at $src (set AFNI_SRC)" >&2; exit 1; }

# Print one C function from a file, from the line matching $2 until its braces balance.
extract_function() {
  awk -v start="$2" '
    !on && $0 ~ start { on = 1 }
    on { print
         n = gsub(/\{/, "{"); m = gsub(/\}/, "}")
         depth += n - m
         if (seen == 0 && n > 0) seen = 1
         if (seen && depth == 0) exit }
  ' "$1"
}

{
  cat <<'PRELUDE'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef unsigned char byte;
typedef struct { byte r, g, b; } rgbyte;
#ifndef MIN
# define MIN(a,b) (((a)<(b)) ? (a) : (b))
# define MAX(a,b) (((a)>(b)) ? (a) : (b))
#endif
int npane_big = 256;                 /* AFNI's variable of the same name */
#define NPANE_BIG npane_big
PRELUDE
  # The macro block of display.h: NBIGMAP_INIT, NBIG_GAP/MBOT/MTOP, AJJ_*, BIGMAP_NAMES.
  sed -n '/^#define NBIGMAP_INIT/,/^int NJ_bigmaps_init/p' "$src/display.h" | sed '$d'
  extract_function "$src/display.c" '^static double mypow'
  extract_function "$src/display.c" '^rgbyte DC_spectrum_ZSS'
  extract_function "$src/display.c" '^rgbyte DC_spectrum_AJJ'
  extract_function "$src/display.c" '^int NJ_bigmaps_init'
  cat <<'MAIN'
int main(void) {
  int sizes[] = { 256, 128, 64, 100, 333 };
  for (int s = 0; s < 5; s++) {
    char **names; rgbyte **maps;
    npane_big = sizes[s];
    if (NJ_bigmaps_init(NBIGMAP_INIT, &names, &maps)) return 1;
    for (int m = 0; m < NBIGMAP_INIT; m++) {
      printf("scale %d %s :", npane_big, names[m]);
      for (int i = 0; i < npane_big; i++)
        printf(" %d,%d,%d", maps[m][i].r, maps[m][i].g, maps[m][i].b);
      printf("\n");
    }
  }
  /* Direct samples of the two spectrum functions, including wrap-around and gamma <= 0. */
  double hues[] = { -400, -90, -0.5, 0, 1, 30, 59.9, 60, 90, 119.99, 120, 121, 180, 239.99, 240,
                    241, 300, 359.5, 360, 360.5, 400, 725 };
  double gammas[] = { 0.0, 0.7, 0.8, 1.0, 1.7 };
  for (int h = 0; h < 22; h++) for (int g = 0; g < 5; g++) {
    rgbyte a = DC_spectrum_AJJ(hues[h], gammas[g]);
    rgbyte z = DC_spectrum_ZSS(hues[h], gammas[g]);
    printf("spec AJJ %g %g : %d,%d,%d\n", hues[h], gammas[g], a.r, a.g, a.b);
    printf("spec ZSS %g %g : %d,%d,%d\n", hues[h], gammas[g], z.r, z.g, z.b);
  }
  return 0;
}
MAIN
} > "$work/scales.c"

# -ffp-contract=off: strict IEEE arithmetic, no fused multiply-add. Clang on Apple
# silicon fuses `a - i*b` by default, which moves a hue by 1 ulp at the 0/360-degree
# wrap and changes ONE byte in a handful of non-default-size scales (never at the
# default 256 or at 128). Strict arithmetic is reproducible on every compiler and is
# what Rust computes; see the roadmap discovery log.
cc -O0 -ffp-contract=off -std=gnu99 -o "$work/scales" "$work/scales.c" -lm

commit="$(git -C "$src/.." log -1 --format=%h 2>/dev/null || echo unknown)"
version="$(afni -ver 2>&1 | head -n 1 || true)"
{
  echo "# AFNI built-in color scales from display.c, computed by AFNI's own code"
  echo "# afni_version: $version"
  echo "# afni_source_commit: $commit"
  echo "# generated: $(date +%Y-%m-%d)"
  echo "# generator: tests/data/regenerate_afni_colorscales.sh"
  echo "# index 0 is the top of the color bar (highest value)"
  "$work/scales"
} > "$out"
echo "wrote $out ($(grep -c '^scale' "$out") scales, $(grep -c '^spec' "$out") spectrum samples)"
