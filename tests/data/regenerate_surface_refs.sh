#!/usr/bin/env bash
# Regenerate the surface fixtures and AFNI/SUMA reference outputs used by
# tests/mesh_conformance.rs and tests/cluster_conformance.rs.
#
#   tests/data/regenerate_surface_refs.sh        (needs AFNI on PATH, and python3)
#
# Meshes (all made here, deterministically, from AFNI's CreateIcosahedron):
#   ico3.asc    642-node regular icosahedron sphere (radius 50)
#   irr3.asc    the same with every coordinate jittered by up to +-3 (irregular)
#   irr2.asc    162-node jittered sphere, for the geometry checks
#   patch2.asc  irr2 with the faces on one side removed: an OPEN surface
# Data (node-wise, deterministic pseudo-random):
#   dA.1D       ~45% of nodes active, signs mixed, the rest exactly 0
#   dB.1D       every node nonzero, values spread over about [-5, 5]
# References, all produced by AFNI programs from those files:
#   refs/irr2.*, refs/patch2.*   SurfaceMetrics (-area -edges -node_normals
#                                -face_normals -boundary_nodes) and SurfMeasures
#                                n_area_A
#   conformance/surfclust.ref    SurfClust tables for a matrix of options
# Normal `cargo test` does not run this; the output is committed.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
out="$here/surfaces"
conf="$here/conformance"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
command -v CreateIcosahedron >/dev/null || { echo "AFNI not found on PATH" >&2; exit 1; }
command -v python3 >/dev/null || { echo "python3 not found" >&2; exit 1; }
export AFNI_DONT_LOGFILE=YES
mkdir -p "$out/refs" "$conf"

# --- meshes and data ---------------------------------------------------------
( cd "$work"
  CreateIcosahedron -rd 3 -rad 50 -prefix ico3 >/dev/null 2>&1
  CreateIcosahedron -rd 2 -rad 50 -prefix ico2 >/dev/null 2>&1
  python3 - <<'PY'
# A tiny linear congruential generator, so output never depends on the Python version.
class Lcg:
    def __init__(self, seed): self.s = seed
    def next(self):
        self.s = (self.s * 1103515245 + 12345) & 0x7fffffff
        return self.s / 0x7fffffff

def read_asc(name):
    lines = open(name).read().split("\n")
    nv, nf = map(int, lines[1].split())
    V = [list(map(float, l.split()[:3])) for l in lines[2:2 + nv]]
    F = [list(map(int, l.split()[:3])) for l in lines[2 + nv:2 + nv + nf]]
    return V, F

def write_asc(name, V, F):
    used = sorted(set(i for f in F for i in f))
    new = {old: k for k, old in enumerate(used)}
    with open(name, "w") as fh:
        fh.write("#!ascii version in FreeSurfer format (SUMA)\n%d %d\n" % (len(used), len(F)))
        for old in used: fh.write("%.6f  %.6f  %.6f  0\n" % tuple(V[old]))
        for f in F: fh.write("%d %d %d 0\n" % tuple(new[i] for i in f))

def jitter(V, seed):
    r = Lcg(seed)
    return [[c + (r.next() * 2 - 1) * 3 for c in v] for v in V]

V3, F3 = read_asc("ico3.asc"); write_asc("ico3_out.asc", V3, F3)
write_asc("irr3_out.asc", jitter(V3, 7), F3)
V2, F2 = read_asc("ico2.asc"); V2 = jitter(V2, 11)
write_asc("irr2_out.asc", V2, F2)
cen = lambda f: [sum(V2[i][k] for i in f) / 3 for k in range(3)]
write_asc("patch2_out.asc", V2, [f for f in F2 if cen(f)[0] < 12])

r = Lcg(12345)
with open("dA.1D", "w") as fh:
    for _ in range(642):
        u, v = r.next(), (r.next() - 0.5) * 8
        fh.write("%g\n" % (0.0 if u < 0.55 or abs(v) < 1e-3 else round(v, 4)))
r = Lcg(777)
with open("dB.1D", "w") as fh:
    for _ in range(642):
        v = sum(r.next() for _ in range(4)) - 2   # roughly normal, spread about [-2, 2]
        v = round(v * 2.4, 4)
        fh.write("%g\n" % (v if abs(v) >= 1e-3 else 0.01))
PY
)
for m in ico3 irr3 irr2 patch2; do cp "$work/${m}_out.asc" "$out/$m.asc"; done
cp "$work/dA.1D" "$work/dB.1D" "$out/"

# --- geometry references (SurfaceMetrics, SurfMeasures) -------------------------
spec_for() { sed "s/ico3.asc/$1.asc/" "$work/ico3.spec" > "$work/$1.spec"; }
for m in irr2 patch2; do
  ( cd "$out" && SurfaceMetrics -i "$m.asc" -area -edges -node_normals -face_normals \
      -boundary_nodes -prefix "$work/$m" >/dev/null 2>&1 )
  spec_for "$m"; cp "$out/$m.asc" "$work/"
  ( cd "$work" && SurfMeasures -spec "$m.spec" -surf_A "$m.asc" -func n_area_A -out_1D "$m.nodearea.1D" >/dev/null 2>&1 )
  cp "$work/$m.area" "$out/refs/$m.area"
  cp "$work/$m.edges" "$out/refs/$m.edges"
  cp "$work/$m.NodeNormSeg.1D" "$out/refs/$m.NodeNormSeg.1D"
  cp "$work/$m.TriNormSeg.1D" "$out/refs/$m.TriNormSeg.1D"
  cp "$work/$m.nodearea.1D" "$out/refs/$m.nodearea.1D"
done
cp "$work/patch2.boundarynodes.1D.dset" "$out/refs/patch2.boundarynodes"
# Strip the volatile history lines so regenerating does not churn the diff.
sed -i.bak '/^#History/d' "$out"/refs/* && rm -f "$out"/refs/*.bak

# --- SurfClust references -------------------------------------------------------
# name | surface | data | SurfClust options  => rows (columns 0-22, comma-separated, ';' between rows)
cases=(
  "rings1_A|irr3|dA|-rmm -1"
  "rings2_A|irr3|dA|-rmm -2"
  "rings3_A|irr3|dA|-rmm -3"
  "rings1_A_regular|ico3|dA|-rmm -1"
  "rings2_A_regular|ico3|dA|-rmm -2"
  "rings1_B_thresh|irr3|dB|-rmm -1 -thresh 1.6"
  "rings1_B_athresh|irr3|dB|-rmm -1 -athresh 1.8"
  "rings2_B_athresh|irr3|dB|-rmm -2 -athresh 2.2"
  "rings1_B_athresh_regular|ico3|dB|-rmm -1 -athresh 1.8"
  "rings1_B_inrange|irr3|dB|-rmm -1 -in_range 0.8 2.2"
  "rings1_B_exrange|irr3|dB|-rmm -1 -ex_range -1.5 1.5"
  "rings2_B_exrange|ico3|dB|-rmm -2 -ex_range -1.5 1.5"
  "minnodes_B|irr3|dB|-rmm -1 -athresh 1.4 -n 4"
  "minarea_B|irr3|dB|-rmm -1 -athresh 1.4 -amm2 300"
  "negarea_nodes_B|irr3|dB|-rmm -1 -athresh 1.4 -amm2 -5"
  "sort_nodes_B|irr3|dB|-rmm -1 -athresh 1.4 -sort_n_nodes"
  "sort_none_B|irr3|dB|-rmm -1 -athresh 1.4 -sort_none"
  "mm_small_B|irr3|dB|-rmm 6 -athresh 2.0"
  "mm_mid_B|irr3|dB|-rmm 9 -athresh 2.0"
  "mm_wide_B|irr3|dB|-rmm 14 -athresh 2.2"
  "mm_A|irr3|dA|-rmm 7"
  "mm_regular_B|ico3|dB|-rmm 8 -athresh 2.0"
)
{
  echo "# SUMA SurfClust tables for the committed meshes and data"
  echo "# afni_version: $(afni -ver 2>&1 | head -n 1)"
  echo "# generated: $(date +%Y-%m-%d)"
  echo "# generator: tests/data/regenerate_surface_refs.sh"
  echo "# fields: name | surface | data | options => row;row;... (23 columns per row, SurfClust's own precision)"
  for c in "${cases[@]}"; do
    IFS='|' read -r name surf data opts <<<"$c"
    # shellcheck disable=SC2086  # intentional word splitting of the option string
    table="$(cd "$out" && SurfClust -i "$surf.asc" -input "$data.1D" 0 $opts -no_cent 2>/dev/null \
      | awk '/^ *[0-9]+ +[0-9]+ +[0-9.]+ / { gsub(/^ +/, ""); gsub(/ +/, ","); printf "%s%s", (n++ ? ";" : ""), $0 }' || true)"
    echo "$name | $surf | $data | $opts => $table"
  done
} > "$conf/surfclust.ref"
echo "wrote surfaces in $out and $conf/surfclust.ref ($(grep -vc '^#' "$conf/surfclust.ref") cases)"
