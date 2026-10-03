#!/usr/bin/env bash
# Regenerate tests/data/conformance/instacorr.ref: seed correlation done with AFNI's
# OWN functions, in the order SUMA's SUMA_dot.c calls them.
#
#   AFNI_SRC=~/Documents/Programming/afni tests/data/regenerate_instacorr_refs.sh
#
# SUMA's InstaCorr lives in the SUMA program (it needs a running SUMA, a surface
# viewer and a NIML dataset), so it cannot be driven from the command line. Its
# numerical steps are all calls into libmri, though, and the glue is a few lines
# (SUMA_dot.c: SUMA_DotXform_MakeOrts, SUMA_DotDetrendDset, the single-vector path
# around line 415, and THD_normalize). This script compiles a C program against
# AFNI's libmri that makes exactly those calls:
#   whole dataset:  THD_build_polyref(polort+1) [+ extra orts]
#                   THD_bandpass_vectors(qdet = 1)   then THD_normalize per row
#   external seed:  THD_bandpass_vectors(qdet = 0) with the same orts, THD_normalize
#   ROI seed:       mean of the prepared rows, THD_normalize (as the seed blur does)
# Nothing numerical is re-implemented. Needs a C compiler and python3.
#
# Fixture format, one block per case:
#   case <name>
#   params <n> <nrows> <dt> <fbot> <ftop> <polort> <nextra_orts> <band 1/0>
#   row <n numbers>      (nrows lines)  input series
#   ort <n numbers>      (nextra lines) extra regressors
#   seed <row indexes>   the ROI seed (the first one is also the single-row seed)
#   ext <n numbers>      an external seed series
#   ndof <n>             removed dimensions reported by the whole-dataset call
#   prep <n numbers>     (nrows lines) cleaned unit-length rows
#   corr_row / corr_roi / corr_ext <nrows numbers>
#   end
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
afni="${AFNI_SRC:-$HOME/Documents/Programming/afni}"
# REGEN_OUT lets a test write somewhere else (the live tests do, so they never
# overwrite the committed file that the other tests are reading).
out="${REGEN_OUT:-$here/conformance/instacorr.ref}"
lib="$afni/build/targets_built"
[ -f "$afni/src/mrilib.h" ] && [ -d "$lib" ] || { echo "AFNI source/build not found at $afni (set AFNI_SRC)" >&2; exit 1; }
command -v cc >/dev/null && command -v python3 >/dev/null || { echo "need cc and python3" >&2; exit 1; }
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

inc="-I$afni/src -I$afni/src/nifti/nifti2 -I$afni/src/nifti/nifti_clib -I$afni/src/niml -I$afni/src/rickr -I$afni/src/f2c -I$afni/src/f2cdir -I$afni/src/eispack"
for d in "$afni"/src/nifti/*; do [ -d "$d" ] && inc="$inc -I$d"; done
for d in /opt/homebrew/include /opt/X11/include /usr/local/include; do [ -d "$d" ] && inc="$inc -I$d"; done

cat > "$work/harness.c" <<'C'
#include "mrilib.h"
static double dotp(const float *a, const float *b, int n){ double s=0; for(int i=0;i<n;i++) s+=(double)a[i]*(double)b[i]; return s; }
static void readvec(float *x, int n){ for(int i=0;i<n;i++){ double d; scanf("%lf",&d); x[i]=(float)d; } }
static void printvec(const char *tag, const float *x, int n){ printf("%s",tag); for(int i=0;i<n;i++) printf(" %.9g",x[i]); printf("\n"); }
int main(void){
  char tag[64], name[128]; int n,nrows,polort,nextra,band; double dt,fbot,ftop;
  while( scanf("%63s",tag)==1 ){
    scanf("%127s",name);
    scanf("%*s %d %d %lf %lf %lf %d %d %d",&n,&nrows,&dt,&fbot,&ftop,&polort,&nextra,&band);
    float **rows=malloc(sizeof(float*)*nrows), **extra=malloc(sizeof(float*)*(nextra+1));
    printf("case %s\nparams %d %d %.9g %.9g %.9g %d %d %d\n",name,n,nrows,dt,fbot,ftop,polort,nextra,band);
    for(int r=0;r<nrows;r++){ rows[r]=malloc(sizeof(float)*n); scanf("%*s"); readvec(rows[r],n); printvec("row",rows[r],n); }
    for(int e=0;e<nextra;e++){ extra[e]=malloc(sizeof(float)*n); scanf("%*s"); readvec(extra[e],n); printvec("ort",extra[e],n); }
    int nseed; scanf("%*s %d",&nseed); int seeds[16]; printf("seed");
    for(int s=0;s<nseed;s++){ scanf("%d",&seeds[s]); printf(" %d",seeds[s]); } printf("\n");
    float *ext=malloc(sizeof(float)*n); scanf("%*s"); readvec(ext,n); printvec("ext",ext,n);

    /* SUMA_DotXform_MakeOrts: Legendre baseline (polort+1) plus the extra orts. */
    int nref = polort+1; float **ref = nref>0 ? THD_build_polyref(nref,n) : NULL;
    int nort = nref+nextra; float **ort = nort>0 ? malloc(sizeof(float*)*nort) : NULL;
    for(int i=0;i<nref;i++) ort[i]=ref[i];
    for(int e=0;e<nextra;e++) ort[nref+e]=extra[e];
    float fb = band?(float)fbot:0.0f, ft = band?(float)ftop:99999.9f;

    /* The external seed is cleaned first (it needs un-touched orts: the dataset call
       does not modify them, but keep the order explicit). qdet = 0 as in SUMA_dot.c. */
    float *eseed=malloc(sizeof(float)*n); memcpy(eseed,ext,sizeof(float)*n);
    THD_bandpass_vectors(n,1,&eseed,(float)dt,fb,ft,0,nort,ort);
    THD_normalize(n,eseed);

    /* SUMA_DotDetrendDset: qdet = 1, then THD_normalize on every row. */
    int nd = THD_bandpass_vectors(n,nrows,rows,(float)dt,fb,ft,1,nort,ort);
    for(int r=0;r<nrows;r++) THD_normalize(n,rows[r]);
    printf("ndof %d\n",nd);
    for(int r=0;r<nrows;r++) printvec("prep",rows[r],n);

    float *roi=calloc(n,sizeof(float));
    for(int s=0;s<nseed;s++) for(int i=0;i<n;i++) roi[i]+=rows[seeds[s]][i];
    for(int i=0;i<n;i++) roi[i]/=nseed;
    THD_normalize(n,roi);
    printf("corr_row"); for(int r=0;r<nrows;r++) printf(" %.9g",dotp(rows[r],rows[seeds[0]],n)); printf("\n");
    printf("corr_roi"); for(int r=0;r<nrows;r++) printf(" %.9g",dotp(rows[r],roi,n)); printf("\n");
    printf("corr_ext"); for(int r=0;r<nrows;r++) printf(" %.9g",dotp(rows[r],eseed,n)); printf("\n");
    printf("end\n\n");
  }
  return 0;
}
C
cc -O1 $inc "$work/harness.c" -L"$lib" -lmri -Wl,-rpath,"$lib" -o "$work/harness" 2>"$work/cc.log" \
  || { cat "$work/cc.log" >&2; exit 1; }

python3 - > "$work/cases.txt" <<'PY'
import math, random
rng = random.Random(77)
def make_rows(n, nrows):
    # three latent in-band signals plus a slow drift and noise; rows mix them.
    lat = [[math.sin(0.30*i + p) + 0.5*math.cos(0.11*i*(p+1)) for i in range(n)] for p in (0.0, 1.3, 2.9)]
    rows = []
    for r in range(nrows):
        w = [rng.uniform(-1, 1) for _ in range(3)]
        w[r % 3] += 1.5
        rows.append([round(50 + 0.04*i*(1 + 0.1*r) + sum(w[k]*lat[k][i] for k in range(3)) + rng.gauss(0, 0.7), 6) for i in range(n)])
    return rows
def ext(n):
    return [round(math.sin(0.3*i) + 0.4*math.cos(0.8*i) + 0.01*i + rng.gauss(0, 0.2), 6) for i in range(n)]
def case(name, n, dt, fbot, ftop, polort, band=1, extra=0):
    rows = make_rows(n, 12)
    orts = [[round(math.sin(0.21*i + 0.7*k) + 0.3*math.cos(0.6*i), 6) for i in range(n)] for k in range(extra)]
    print("case", name)
    print("params", n, len(rows), dt, fbot, ftop, polort, extra, band)
    for r in rows: print("row", " ".join(repr(v) for v in r))
    for o in orts: print("ort", " ".join(repr(v) for v in o))
    print("seed 3 0 1 2")
    print("ext", " ".join(repr(v) for v in ext(n)))
case("linear_only", 60, 2.0, 0.0, 0.0, -1, band=0)
case("bandpass_no_polort", 120, 1.5, 0.02, 0.2, -1)
case("bandpass_odd_length_no_polort", 127, 2.0, 0.02, 0.2, -1)
case("bandpass_extra_ort", 100, 2.0, 0.02, 0.2, -1, extra=2)
case("suma_default_polort2", 120, 2.0, 0.01, 0.1, 2)
case("polort1_extra_ort", 100, 2.0, 0.02, 0.2, 1, extra=1)
PY
{
  echo "# SUMA InstaCorr steps with AFNI's own functions (regenerate with regenerate_instacorr_refs.sh)"
  echo "# afni_version: $(afni -ver 2>/dev/null | head -1 || echo unknown)"
  "$work/harness" < "$work/cases.txt" 2>/dev/null
} > "$out"
echo "wrote $out with $(grep -c '^case' "$out") cases"
