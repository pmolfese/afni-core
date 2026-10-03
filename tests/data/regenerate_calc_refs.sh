#!/usr/bin/env bash
# Regenerate tests/data/conformance/calc.ref: values of 3dcalc-style expressions
# as AFNI's own evaluator (parser.f, through `1deval`) computes them.
#
#   tests/data/regenerate_calc_refs.sh
#
# Normal `cargo test` never runs this. Each case is one line
#
#     <expression without spaces> [letter=number ...] => <value>
#
# e.g. `step(a-3)*step(b-2) a=4 b=2.5 => 1` means
# `1deval -a=4 -b=2.5 -num 1 -expr 'step(a-3)*step(b-2)'` printed 1.
# `1deval` prints 6 significant digits, so tests compare to ~1e-5 relative.
#
# `absextreme` is deliberately not here: AFNI's scalar evaluator never matches
# it (the 8-character opcode 'ABSEXTRE' is compared with 'ABSEXTREME'), so
# `absextreme(1,-5,3,7)` prints 4, the argument count. afni-core implements the
# documented behavior (largest absolute value); see DIFFERENCES_FROM_AFNI.md.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
out="${REGEN_OUT:-$here/conformance/calc.ref}"
command -v 1deval >/dev/null || { echo "1deval not found on PATH (is AFNI installed?)" >&2; exit 1; }
version="$(afni -ver 2>&1 | head -n 1)"

# expression|assignments (comma separated)
cases=(
  # --- grammar: precedence, associativity, unary signs, spaces, case ---
  '2+3*4|' '(2+3)*4|' '10/4|' '5-3-1|' '8/4/2|' '2^3^2|' '2**3**2|' '-2^2|' '-2**2|'
  '2^-1|' '2*-3|' '+5|' '[2+3]*4|' '2-(-3)|' '-(2+3)^2|' '2*3^2|' '2^3*2|' '1+2+3*4/6-1|'
  '-a^2|a=3' 'a^b^c|a=2,b=3,c=2' '-a-b|a=2,b=3' 'a--b|a=2,b=3' 'a/b/c|a=24,b=4,c=3'
  '1e3|' '1.5e-2|' '.5|' '3.|' '1d2|' 'PI|' 'pi*2|' 'Sin(PI/2)|' 'A+B|a=1,b=2' 'STEP(a)|a=2'
  # --- "designed not to fail" ---
  '1/0|' '0/0|' 'a/b|a=5,b=0' 'sqrt(-4)|' 'sqrt(a)|a=-9' 'log(0)|' 'log(-1)|' 'log(a)|a=-10'
  'log10(-100)|' 'log10(0)|' 'acos(2)|' 'asin(-3)|' 'atanh(2)|' 'acosh(0.5)|' 'asinh(-2)|'
  'asinh(50)|' 'acosh(20)|' 'acosh(3)|' 'atanh(0.5)|' 'exp(1000)|' 'exp(a)|a=2' 'exp(-200)|'
  'sinh(200)|' 'cosh(200)|' 'sinh(1)|' 'cosh(1)|' 'tanh(0.5)|' 'mod(7,0)|' '(0-2)**0.5|'
  '(0-2)**2|' '0**2|' '0**0|' '2**0.5|' 'atan2(0,0)|' 'atan2(1,1)|' 'atan2(-1,-1)|'
  # --- math ---
  'int(-2.7)|' 'int(2.7)|' 'mod(7,3)|' 'mod(-7,3)|' 'mod(7.5,2)|' 'abs(-3)|' 'min(3,2)|'
  'max(3,2)|' 'max(a,b)|a=-1,b=-2' 'sind(30)|' 'cosd(60)|' 'tand(45)|' 'cbrt(-8)|' 'cbrt(a)|a=10'
  'sin(1)|' 'cos(1)|' 'tan(1)|' 'asin(0.5)|' 'acos(0.5)|' 'atan(1)|'
  # --- masks ---
  'step(0)|' 'step(-1)|' 'step(0.001)|' 'step(a-3)|a=3' 'step(a-3)|a=3.0001' 'astep(3,2)|'
  'astep(-3,2)|' 'astep(2,2)|' 'astep(a,b)|a=-1.5,b=1' 'within(2,1,3)|' 'within(3,1,3)|'
  'within(1,1,3)|' 'within(3.1,1,3)|' 'within(a,b,c)|a=5,b=1,c=4' 'rect(0.5)|' 'rect(0.6)|'
  'rect(-0.3)|' 'bool(0)|' 'bool(0.1)|' 'bool(-5)|' 'notzero(a)|a=0' 'iszero(a)|a=0' 'not(0)|' 'not(2)|'
  'equals(2,2)|' 'equals(2,3)|' 'equals(a,b)|a=1.5,b=1.5' 'ifelse(1,5,6)|' 'ifelse(0,5,6)|'
  'ifelse(a-2,10,20)|a=2' 'and(1,1,0)|' 'and(1,2,3)|' 'and(a,b,c)|a=1,b=0,c=1' 'or(0,0,0)|'
  'or(0,0,3)|' 'or(a,b)|a=0,b=0' 'mofn(2,1,0,1)|' 'mofn(3,1,0,1)|' 'mofn(1,0,0,0)|' 'posval(-2)|'
  'posval(2)|' 'ispositive(0)|' 'ispositive(2)|' 'isnegative(-1)|' 'isnegative(0)|' 'tent(0.25)|'
  'tent(2)|' 'bell2(0.25)|' 'bell2(1)|' 'bell2(2)|'
  'step(a-3)*step(b-2)|a=4,b=3' 'step(a-3)*step(b-2)|a=4,b=1' 'step(a-3)+2*step(b-2)|a=4,b=3'
  'step(a-3)*(1-step(b-2))|a=4,b=1' 'a*step(b)|a=7,b=-1' 'astep(a,2)*within(b,0,5)|a=-3,b=2'
  'and(step(a-3),step(b-2))|a=4,b=3' 'or(step(a-3),step(b-2))|a=1,b=3' 'step(9-(a-1)*(a-1)-(b-2)*(b-2))|a=2,b=3'
  'step(9-(a-1)*(a-1)-(b-2)*(b-2))|a=9,b=3'
  # --- functions of several arguments ---
  'median(3,1,2)|' 'median(4,1,3,2)|' 'median(5)|' 'median(a,b,c,d,e)|a=5,b=1,c=4,d=2,e=3' 'mean(1,2,6)|'
  'mean(7)|' 'mean(1,2)|' 'stdev(1,2,3)|' 'stdev(5)|' 'stdev(2,4,4,4,5,5,7,9)|' 'sem(1,2,3)|'
  'mad(1,2,3,4,100)|' 'mad(1,3)|' 'mad(4)|' 'argmax(1,3,2)|' 'argmax(0,0)|' 'argmax(a,b,c)|a=0,b=-1,c=0'
  'argnum(1,0,2)|' 'argnum(0,0)|' 'choose(2,10,20,30)|' 'choose(4,10,20,30)|' 'choose(0,10,20)|'
  'amongst(2,1,2,3)|' 'amongst(5,1,2,3)|' 'orstat(1,5,3,9)|' 'orstat(2,5,3,9)|' 'orstat(3,5,3,9)|'
  'orstat(9,5,3,9)|' 'pairmin(3,2,7,5,-1,-2,-3,-4)|' 'pairmax(1,5,2,10,50,20)|' 'pairmax(1,1,2,9)|'
  'lmode(1,2,2,3,3)|' 'hmode(1,2,2,3,3)|' 'lmode(1,2,3)|' 'hmode(1,2,3)|' 'lmode(5)|'
  'minabove(2,1,3,5)|' 'minabove(9,1,3,5)|' 'maxbelow(4,1,3,5)|' 'maxbelow(0,1,3,5)|' 'extreme(1,-5,3)|'
  'extreme(0,0)|' 'extreme(a,b)|a=2,b=-2'
)

{
  echo "# AFNI conformance fixture for 3dcalc-style expressions (parser.f, via 1deval)."
  echo "# afni_version: $version"
  echo "# generated: $(date +%F)"
  echo "# generator: tests/data/regenerate_calc_refs.sh"
  echo "# precision: 1deval prints 6 significant digits"
  for c in "${cases[@]}"; do
    expr="${c%%|*}"
    vars="${c#*|}"
    args=()
    shown=""
    if [ -n "$vars" ]; then
      IFS=',' read -ra pairs <<< "$vars"
      for p in "${pairs[@]}"; do args+=("-${p%%=*}=${p#*=}"); shown+=" ${p}"; done
    fi
    value="$(1deval "${args[@]+"${args[@]}"}" -num 1 -expr "$expr" 2>/dev/null | tail -n 1 | tr -d ' ')"
    [ -n "$value" ] || { echo "1deval failed for: $expr ($vars)" >&2; exit 1; }
    echo "${expr}${shown} => ${value}"
  done
} > "$out"
echo "wrote $out"
