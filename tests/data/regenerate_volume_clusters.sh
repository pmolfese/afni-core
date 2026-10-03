#!/usr/bin/env bash
# Regenerate tests/data/conformance/volume_clusters.ref: voxel clustering done by
# AFNI's OWN `3dClusterize`, on small synthetic volumes made here with `3dUndump`.
# Used by tests/volume_cluster_conformance.rs to check afni_core::volume_cluster.
#
#   tests/data/regenerate_volume_clusters.sh        (needs AFNI and python3 on PATH)
#
# Each case records the volume (threshold values, optional separate data values,
# optional mask), the 3dClusterize options, the map of cluster numbers written by
# `-pref_map`, and the rows of the printed report. The affine is read back from
# AFNI itself (`3dmaskdump -xyz`), so coordinates follow AFNI's convention exactly.
# One case is made oblique (an IJK_TO_DICOM_REAL attribute is added): 3dClusterize
# still reports CARDINAL grid coordinates, which the test checks.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
# REGEN_OUT lets a test write somewhere else (the live tests do, so they never
# overwrite the committed file that the other tests are reading).
out="${REGEN_OUT:-$here/conformance/volume_clusters.ref}"
command -v 3dClusterize >/dev/null || { echo "3dClusterize not found (is AFNI installed?)" >&2; exit 1; }
command -v python3 >/dev/null || { echo "python3 not found" >&2; exit 1; }
export AFNI_DONT_LOGFILE=YES
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

ver="$(afni -ver 2>/dev/null | head -1 || echo unknown)"
OUT="$out" WORK="$work" VER="$ver" python3 - <<'PY'
import os, random, struct, subprocess

work, out, ver = os.environ["WORK"], os.environ["OUT"], os.environ["VER"]

def f32(x):
    """Round to the nearest 32-bit float (what AFNI stores)."""
    return struct.unpack("f", struct.pack("f", x))[0]

def fmt(x):
    return "%.9g" % x

def sh(args, **kw):
    return subprocess.run(args, cwd=work, capture_output=True, text=True, check=True, **kw)

def undump(prefix, dims, values, orient, dxyz):
    """Make a float dataset from a flat list of values (nonzero ones are written)."""
    nx, ny, nz = dims
    with open(os.path.join(work, prefix + ".txt"), "w") as f:
        for n, v in enumerate(values):
            if v != 0:
                f.write("%d %d %d %s\n" % (n % nx, (n // nx) % ny, n // (nx * ny), fmt(v)))
    sh(["3dUndump", "-dimen", *map(str, dims), "-orient", orient, "-prefix", prefix,
        "-ijk", "-datum", "float", "-overwrite", prefix + ".txt"])
    sh(["3drefit", "-xdel", str(dxyz[0]), "-ydel", str(dxyz[1]), "-zdel", str(dxyz[2]),
        prefix + "+orig"])

def dump(prefix):
    """Every voxel's value, in index order (x fastest)."""
    res = sh(["3dmaskdump", "-noijk", prefix + "+orig"])
    return [float(l.split()[0]) for l in res.stdout.splitlines()
            if l and not l.startswith("++")]

def affine(prefix, dims):
    """The 3x4 voxel-to-DICOM matrix, learned from the coordinates AFNI reports."""
    res = sh(["3dmaskdump", "-noijk", "-xyz", prefix + "+orig"])
    pts = [list(map(float, l.split()[:3])) for l in res.stdout.splitlines()
           if l and not l.startswith("++")]
    nx, ny, _ = dims
    o, di, dj, dk = pts[0], pts[1], pts[nx], pts[nx * ny]
    cols = [[di[a] - o[a] for a in range(3)], [dj[a] - o[a] for a in range(3)],
            [dk[a] - o[a] for a in range(3)]]
    return [cols[0][r:r+1] + cols[1][r:r+1] + cols[2][r:r+1] + [o[r]] for r in range(3)]

def parse_report(text):
    """Rows of the cluster table: 15 numbers each (Nvox, CM x3, min/max x6, mean, sem,
    max int, MI x3)."""
    rows = []
    for line in text.splitlines():
        if line.startswith("#") or line.startswith("++") or not line.strip():
            continue
        parts = line.split()
        if len(parts) == 16:
            rows.append(parts)
    return rows

rng = random.Random(20261002)

def noise(n, scale):
    # Sum of four uniforms: close enough to bell-shaped; rounded to a float.
    return [f32(scale * (sum(rng.uniform(-1, 1) for _ in range(4)) / 2)) for _ in range(n)]

cases = []

def case(name, dims, thr, options, dat=None, mask=None, orient="LPI",
         dxyz=(2.0, 2.5, 3.0), oblique=False):
    cases.append(dict(name=name, dims=dims, thr=thr, dat=dat, mask=mask, options=options,
                      orient=orient, dxyz=dxyz, oblique=oblique))

D = (12, 11, 9)
N = D[0] * D[1] * D[2]
thr = noise(N, 2.2)
# Plant values exactly on the thresholds so the inclusive ends are tested.
for pos, v in ((5, 2.0), (6, -2.0), (50, 2.0), (51, 2.0), (200, -2.0), (201, -3.0)):
    thr[pos] = v
dat = noise(N, 3.0)
for pos in range(0, N, 7):          # zeros in the data where the threshold may pass
    dat[pos] = 0.0
mask = [1 if rng.random() < 0.75 else 0 for _ in range(N)]

for nn in (1, 2, 3):
    case("right_nn%d" % nn, D, thr, ["-NN", str(nn), "-1sided", "RIGHT_TAIL", "2.0"])
case("left_nn1", D, thr, ["-NN", "1", "-1sided", "LEFT_TAIL", "-2.0"])
case("twosided_nn2", D, thr, ["-NN", "2", "-2sided", "-2.0", "2.0"])
case("bisided_nn2", D, thr, ["-NN", "2", "-bisided", "-2.0", "2.0"])
case("bisided_nn3", D, thr, ["-NN", "3", "-bisided", "-1.5", "1.5"])
case("bisided_nvox", D, thr, ["-NN", "2", "-bisided", "-1.5", "1.5", "-clust_nvox", "4"])
case("asymmetric_bisided", D, thr, ["-NN", "1", "-bisided", "-1.0", "2.5"])
case("within_range", D, thr, ["-NN", "1", "-within_range", "-1.0", "1.0"])
case("idat", D, thr, ["-NN", "2", "-bisided", "-2.0", "2.0"], dat=dat)
case("mask", D, thr, ["-NN", "2", "-bisided", "-2.0", "2.0"], mask=mask)
case("idat_mask_nvox", D, thr, ["-NN", "1", "-2sided", "-1.5", "1.5", "-clust_nvox", "3"],
     dat=dat, mask=mask)
case("oblique", D, thr, ["-NN", "2", "-bisided", "-2.0", "2.0"], oblique=True)
case("rai_voxels", (10, 9, 8), noise(720, 2.2), ["-NN", "3", "-bisided", "-2.0", "2.0"],
     orient="RAI", dxyz=(1.5, 1.5, 4.0))
# Many one-voxel clusters of both signs, to pin down the order of equal sizes.
ties = [0.0] * N
for n in range(0, N, 5):
    ties[n] = 3.0 if (n // 5) % 3 else -3.0
case("ties_bisided", D, ties, ["-NN", "1", "-bisided", "-2.0", "2.0"])
case("ties_twosided", D, ties, ["-NN", "1", "-2sided", "-2.0", "2.0"])
case("no_clusters", D, thr, ["-NN", "1", "-1sided", "RIGHT_TAIL", "50.0"])

lines = ["# Voxel clusters from AFNI's 3dClusterize (regenerate with regenerate_volume_clusters.sh)",
         "# afni_version: " + ver, ""]
for c in cases:
    nx, ny, nz = c["dims"]
    thr_p, dat_p, mask_p, map_p = ("thr_" + c["name"], "dat_" + c["name"],
                                   "msk_" + c["name"], "map_" + c["name"])
    undump(thr_p, c["dims"], c["thr"], c["orient"], c["dxyz"])
    aff = affine(thr_p, c["dims"])
    args = ["3dClusterize", "-inset", thr_p + "+orig", "-ithr", "0"]
    if c["dat"] is not None:
        undump(dat_p, c["dims"], c["dat"], c["orient"], c["dxyz"])
        # One dataset with two sub-bricks: [0] threshold, [1] data.
        sh(["3dTcat", "-prefix", "both_" + c["name"], thr_p + "+orig", dat_p + "+orig"])
        args = ["3dClusterize", "-inset", "both_" + c["name"] + "+orig", "-ithr", "0",
                "-idat", "1"]
    if c["mask"] is not None:
        undump(mask_p, c["dims"], [float(m) for m in c["mask"]], c["orient"], c["dxyz"])
        args += ["-mask", mask_p + "+orig"]
    if c["oblique"]:
        target = ("both_" + c["name"]) if c["dat"] is not None else thr_p
        sh(["3drefit", "-atrfloat", "IJK_TO_DICOM_REAL",
            "-1.9 0.3 0.2 4.0  0.2 -2.4 0.5 -3.0  0.1 0.4 2.9 1.5", target + "+orig"])
    args += c["options"] + ["-pref_map", map_p, "-overwrite"]
    res = subprocess.run(args, cwd=work, capture_output=True, text=True)
    rows = parse_report(res.stdout)
    labels = [int(float(v)) for v in dump(map_p)] if os.path.exists(
        os.path.join(work, map_p + "+orig.HEAD")) else [0] * (nx * ny * nz)
    lines.append("case " + c["name"])
    lines.append("dims %d %d %d" % c["dims"])
    lines.append("affine " + " ".join(fmt(v) for row in aff for v in row))
    lines.append("options " + " ".join(c["options"]))
    lines.append("thr " + " ".join(fmt(v) for v in c["thr"]))
    lines.append("dat " + (" ".join(fmt(v) for v in c["dat"]) if c["dat"] is not None else "-"))
    lines.append("mask " + (" ".join(str(m) for m in c["mask"]) if c["mask"] is not None else "-"))
    lines.append("map " + " ".join(str(v) for v in labels))
    for r in rows:
        lines.append("row " + " ".join(r))
    lines.append("end")
    lines.append("")

with open(out, "w") as f:
    f.write("\n".join(lines))
print("wrote", out, "with", len(cases), "cases")
PY
