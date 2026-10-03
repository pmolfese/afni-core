#!/usr/bin/env bash
# Regenerate tests/data/conformance/nifti_stats.ref: AFNI's own NIfTI statistical
# library (src/nifti/nifticdf, the code behind 3dcalc/p2dsetstat/3dFDR) evaluated
# at full double precision.
#
# Normal `cargo test` never runs this; the output is committed. Run it when you
# add cases or move to a new AFNI release, then review and commit the diff.
#
#   AFNI_SRC=~/Documents/Programming/afni tests/data/regenerate_nifti_stats.sh
#
# Why build it ourselves: the shipped `nifti_stats` demo is not installed with
# AFNI binaries, and it prints only 9 significant digits. We compile the same
# source with `%.17g` so the fixture holds every digit AFNI computes.
#
# Line format (same as the other conformance fixtures): `<arguments> => <value>`
#   <val> CODE p1 p2 p3        cdf          = P(X <= val)
#   -q <val> CODE p1 p2 p3     upper tail   = P(X >  val)
#   -1 <p> CODE p1 p2 p3       inverse cdf  = x with P(X <= x) = p
# Parameters use the NIfTI convention (CORREL takes one parameter, the dof).
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
src="${AFNI_SRC:-$HOME/Documents/Programming/afni}/src/nifti"
out="$here/conformance/nifti_stats.ref"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -f "$src/nifticdf/nifti_stats.c" ] || { echo "AFNI source not found at $src (set AFNI_SRC)" >&2; exit 1; }

sed 's/%\.9g/%.17g/' "$src/nifticdf/nifti_stats.c" > "$work/nifti_stats17.c"
cp "$src/nifticdf/nifticdf.c" "$src/nifticdf/nifticdf.h" "$src/nifticdf/nifticdf_version.h" "$src/nifti2/nifti1.h" "$work/"
cc -O2 -I"$work" "$work/nifti_stats17.c" "$work/nifticdf.c" -lm -o "$work/nifti_stats17"

version="$(afni -ver 2>&1 | head -n 1 || true)"
commit="$(git -C "$src/../.." log -1 --format=%h 2>/dev/null || echo unknown)"

# "CODE p1 p2 p3 : v1 v2 v3 ..." -- the values to evaluate for that distribution.
cases=(
  "CORREL 18          : -0.9 -0.5 -0.1 0.3 0.6 0.95"
  "CORREL 2.5         : -0.8 0.2 0.7"
  "TTEST 10           : -10 -3 -1.0 0.5 2 5 10"
  "TTEST 2.5          : -4 -0.5 1 6"
  "TTEST 100          : -2 1.96 4"
  "FTEST 3 20         : 0.2 0.5 1 3 10"
  "FTEST 1 1          : 0.5 1 4"
  "FTEST 10.5 7.25    : 0.8 2 5"
  "ZSCORE             : -8 -5 -1.96 0 1 3 8"
  "CHISQ 1            : 0.5 3.84 12"
  "CHISQ 5            : 1 5 20 60"
  "CHISQ 0.5          : 0.1 2"
  "BETA 2 5           : 0.05 0.3 0.7 0.95"
  "BETA 0.5 0.5       : 0.01 0.5 0.99"
  "GAMMA 3 2          : 0.1 1 3 10"
  "GAMMA 0.5 1        : 0.05 1 8"
  "BINOM 10 0.3       : 0 1 2 3 4 5 6 7 8 9 10"
  "BINOM 100 0.5      : 30 40 50 60 70"
  "POISSON 3          : 0 1 2 3 5 10 20"
  "POISSON 50         : 30 40 50 60 80"
  "NORMAL 1 2         : -6 -3 0 1 4 9"
  "LOGISTIC 0.5 2     : -20 -3 0.5 4 30"
  "LAPLACE 0.5 2      : -15 -2 0.5 3 20"
  "UNIFORM -1 3       : -2 -1.0 0 2 3 4"
  "WEIBULL 0.5 2 1.5  : 0.5 1 2 5 12"
  "CHI 3              : 0.3 1 2 4"
  "INVGAUSS 2 3       : 0.5 1 2 5 10"
  "EXTVAL 0.5 2       : -4 -2 0 1 4 12"
  "CHISQ_NONC 3 2     : 0.5 2 5 10 25"
  "CHISQ_NONC 10 15   : 5 15 25 40 70"
  "CHISQ_NONC 1.5 0.5 : 0.1 1 4 12"
  "FTEST_NONC 3 20 2  : 0.3 1 3 6 15"
  "FTEST_NONC 5 12.5 8 : 0.5 2 5 10 20"
  "TTEST_NONC 10 1.5  : -4 -1.0 0 1 2.5 5 9"
  "TTEST_NONC 5 -2    : -9 -3 -0.5 0 1 3"
  "TTEST_NONC 30 4    : 0.5 2 4 6 9"
  "TTEST_NONC 3.5 0   : -3 0.5 3"
)
# Probabilities for the inverse. Discrete distributions are skipped: AFNI inverts
# a continuous extension, while afni-core returns an integer quantile (documented).
inverse_p=(0.001 0.05 0.3 0.5 0.95 0.999)
discrete="BINOM POISSON"

{
  echo "# AFNI/NIfTI statistical library (nifticdf.c) reference values, full double precision."
  echo "# afni_version: $version"
  echo "# afni_source_commit: $commit"
  echo "# generated: $(date +%Y-%m-%d)"
  echo "# generator: tests/data/regenerate_nifti_stats.sh"
  echo "# note: CDFLIB-derived values; trust to ~1e-9 relative, not 1e-16"
  for c in "${cases[@]}"; do
    head="${c%%:*}"; vals="${c#*:}"
    head="$(echo $head)"   # collapse the alignment padding
    # NOTE: a bare "-1" value would be read as nifti_stats' inverse flag, so
    # cases must write it as -1.0.
    for v in $vals; do
      echo "$v $head => $("$work/nifti_stats17" "$v" $head)"
      echo "-q $v $head => $("$work/nifti_stats17" -q "$v" $head)"
    done
    code="${head%% *}"
    case " $discrete " in *" $code "*) continue ;; esac
    for p in "${inverse_p[@]}"; do
      echo "-1 $p $head => $("$work/nifti_stats17" -1 "$p" $head)"
    done
  done
} > "$out"

echo "wrote $out ($(grep -vc '^#' "$out") cases)"
