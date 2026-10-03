#!/usr/bin/env bash
# Regenerate tests/data/conformance/signal_bandpass.ref: AFNI's OWN
# `THD_bandpass_vectors` (thd_bandpass.c) run on synthetic series, including the
# number of dimensions it reports removing. Used by tests/signal_conformance.rs.
#
#   AFNI_SRC=~/Documents/Programming/afni tests/data/regenerate_signal_refs.sh
#
# `THD_bandpass_vectors` is a library function with no command-line program that
# prints its return value (1dBandpass hides it and rounds to 6 decimals), so this
# script compiles a small C program against AFNI's own libmri (built in
# $AFNI_SRC/build/targets_built) and calls the function directly. Nothing is
# re-implemented. Needs a C compiler and python3.
#
# Fixture format, one block per case:
#   case <name>
#   params <n> <nvec> <dt> <fbot> <ftop> <qdet> <nort> <band: 1 filter / 0 none>
#   in  <n numbers>     (nvec lines, then nort ort lines)
#   ndof <return value>
#   out <n numbers>     (nvec lines)
#   end
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
afni="${AFNI_SRC:-$HOME/Documents/Programming/afni}"
# REGEN_OUT lets a test write somewhere else (the live tests do, so they never
# overwrite the committed file that the other tests are reading).
out="${REGEN_OUT:-$here/conformance/signal_bandpass.ref}"
lib="$afni/build/targets_built"
[ -f "$afni/src/mrilib.h" ] && [ -d "$lib" ] || { echo "AFNI source/build not found at $afni (set AFNI_SRC)" >&2; exit 1; }
command -v cc >/dev/null && command -v python3 >/dev/null || { echo "need cc and python3" >&2; exit 1; }
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

inc="-I$afni/src -I$afni/src/nifti/nifti2 -I$afni/src/nifti/nifti_clib -I$afni/src/niml -I$afni/src/rickr -I$afni/src/f2c -I$afni/src/f2cdir -I$afni/src/eispack"
for d in "$afni"/src/nifti/*; do [ -d "$d" ] && inc="$inc -I$d"; done
for d in /opt/homebrew/include /opt/X11/include /usr/local/include; do [ -d "$d" ] && inc="$inc -I$d"; done

cat > "$work/harness.c" <<'C'
/* Reads cases on stdin, runs THD_bandpass_vectors, prints the fixture block. */
#include "mrilib.h"
int main(void){
  char tag[64], name[128]; int n,nvec,qdet,nort,band; double dt,fbot,ftop;
  while( scanf("%63s",tag)==1 ){
    if( strcmp(tag,"case")!=0 ) return 1;
    scanf("%127s",name);
    scanf("%*s %d %d %lf %lf %lf %d %d %d",&n,&nvec,&dt,&fbot,&ftop,&qdet,&nort,&band);
    float **vec=malloc(sizeof(float*)*nvec), **ort=nort?malloc(sizeof(float*)*nort):NULL;
    printf("case %s\nparams %d %d %.9g %.9g %.9g %d %d %d\n",name,n,nvec,dt,fbot,ftop,qdet,nort,band);
    for(int v=0; v<nvec+nort; v++){
      float *x=malloc(sizeof(float)*n);
      scanf("%*s");
      for(int i=0;i<n;i++){ double d; scanf("%lf",&d); x[i]=(float)d; }
      if(v<nvec) vec[v]=x; else ort[v-nvec]=x;
      printf("in");
      for(int i=0;i<n;i++) printf(" %.9g",x[i]);
      printf("\n");
    }
    /* "no band" is how SUMA asks for it: fbot 0, ftop huge. */
    int nd = THD_bandpass_vectors(n,nvec,vec,(float)dt,band?(float)fbot:0.0f,band?(float)ftop:99999.9f,qdet,nort,ort);
    printf("ndof %d\n",nd);
    for(int v=0; v<nvec; v++){
      printf("out");
      for(int i=0;i<n;i++) printf(" %.9g",vec[v][i]);
      printf("\n");
    }
    printf("end\n\n");
  }
  return 0;
}
C
cc -O1 $inc "$work/harness.c" -L"$lib" -lmri -Wl,-rpath,"$lib" -o "$work/harness" 2>"$work/cc.log" \
  || { cat "$work/cc.log" >&2; exit 1; }

python3 - > "$work/cases.txt" <<'PY'
import math, random
rng = random.Random(9)
def series(n, k):
    # a drift, a few tones and noise, rounded to float32-friendly decimals
    return [round(0.02*k*i + 2*math.sin(0.11*i*(k+1)) + math.cos(0.4*i) + rng.gauss(0, 1), 6) for i in range(n)]
def poly(n, count):
    step = 2.0/(n-1)
    out = []
    for m in range(count):
        col = []
        for i in range(n):
            x = step*i - 1.0
            p0, p1 = 1.0, x
            val = 1.0 if m == 0 else x
            for kk in range(2, m+1):
                p0, p1 = p1, ((2*kk-1)*x*p1 - (kk-1)*p0)/kk
                val = p1
            col.append(round(val, 6))
        out.append(col)
    return out
cases = []
def case(name, n, nvec, dt, fbot, ftop, qdet, band=1, orts=None):
    orts = orts or []
    cases.append((name, n, nvec, dt, fbot, ftop, qdet, band, [series(n, k) for k in range(nvec)], orts))
case("bandpass_even_quadratic", 100, 3, 2.0, 0.01, 0.1, 2)
case("bandpass_odd_length", 127, 2, 2.0, 0.01, 0.1, 2)
case("bandpass_prime_plus_one", 137, 3, 2.0, 0.01, 0.1, 2)
case("bandpass_linear", 100, 2, 1.5, 0.02, 0.2, 1)
case("bandpass_mean_only", 90, 2, 1.0, 0.03, 0.2, 0)
case("bandpass_no_detrend", 80, 2, 1.0, 0.03, 0.2, -1)
case("lowpass", 120, 3, 2.0, 0.0, 0.08, 2)
case("highpass_open_top", 120, 3, 2.0, 0.02, 99999.9, 2)
case("top_above_nyquist", 100, 2, 2.0, 0.05, 1.0, 2)
case("half_bin_edges", 100, 2, 1.0, 0.025, 0.115, 2)
case("odd_dt", 111, 4, 0.72, 0.01, 0.4, 2)
case("collapsed_band", 100, 2, 1.0, 0.0301, 0.0302, 2)
case("short_nine", 9, 2, 1.0, 0.05, 0.4, 2)
case("no_filter_linear", 60, 2, 2.0, 0.0, 0.0, 1, band=0)
case("no_filter_quadratic", 60, 3, 2.0, 0.0, 0.0, 2, band=0)
n = 100
# Orts that survive the filter: independent in-band regressors. Results must agree.
def inband_orts(n):
    return [[round(math.sin(0.20*i + ph) + 0.5*math.cos(0.45*i*ph), 6) for i in range(n)] for ph in (0.3, 1.1, 2.0)]
inband = inband_orts(n)
case("orts_independent_bandpass", n, 3, 2.0, 0.01, 0.1, 1, orts=inband)
# Legendre orts with a linear detrend (SUMA's own setup): the constant and linear
# columns are reduced to ROUNDING NOISE by the filter, and AFNI's pseudo-inverse then
# removes those noise directions too, which differ from run to run of 32-bit
# arithmetic. These are checked loosely; see the test.
case("degenerate_orts_polynomial_bandpass", n, 3, 2.0, 0.01, 0.1, 1, orts=poly(n, 3))
# The same orts minus the two the filter annihilates: this one must agree exactly.
case("orts_legendre_quadratic_only", n, 3, 2.0, 0.01, 0.1, 1, orts=poly(n, 3)[2:])
n = 101
extra = [round(math.sin(0.23*i) + 0.3*math.cos(0.9*i), 6) for i in range(n)]
case("degenerate_orts_poly_plus_extra", n, 2, 2.0, 0.01, 0.1, 1, orts=poly(n, 3) + [extra])
n = 80
case("degenerate_orts_without_filter", n, 2, 2.0, 0.0, 0.0, 1, band=0, orts=poly(n, 3))
case("orts_extra_no_detrend", n, 2, 2.0, 0.02, 0.2, -1, orts=inband_orts(n)[:2])
for name, n, nvec, dt, fbot, ftop, qdet, band, vecs, orts in cases:
    assert all(len(v) == n for v in vecs + orts), name   # every line must have n numbers
    print("case", name)
    print("params", n, nvec, dt, fbot, ftop, qdet, len(orts), band)
    for v in vecs + orts:
        print("in", " ".join(repr(x) for x in v))
PY
{
  echo "# AFNI THD_bandpass_vectors on synthetic series (regenerate with regenerate_signal_refs.sh)"
  echo "# afni_version: $(afni -ver 2>/dev/null | head -1 || echo unknown)"
  "$work/harness" < "$work/cases.txt" 2>/dev/null
} > "$out"
echo "wrote $out with $(grep -c '^case' "$out") cases"
