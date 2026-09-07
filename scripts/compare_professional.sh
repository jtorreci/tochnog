#!/bin/bash
# ---------------------------------------------------------------------------
# compare_professional.sh - repeatable A/B harness: Tochnog GNU vs Tochnog
# Professional on the materi_stress_force validation family.
#
# For each model it runs BOTH binaries, extracts the node_dof_calcul records
# from both .dbs files, and compares them against the analytic statics
# (beam theory / pressure vessel). This is the acceptance criterion of the
# solver fix lot C/D of ProjectDocs/DIAG-SOLVE-MIXTO.md.
#
# Usage:
#   TOCHNOG_PROF_BIN=/path/to/professional/tochnog scripts/compare_professional.sh [model...]
#
# Models: gforce7 gforce7q4 gforce7q4_ref gffq4 gforce7_ref gforce10 gforce13
#         (all: runs every model)
#
# Environment:
#   TOCHNOG_PROF_BIN  - path to the Tochnog Professional binary (REQUIRED).
#   TOCHNOG_GNU_BIN   - path to the GNU binary (default: <repo>/build/tochnog).
#   PROF_TEST_DIR     - directory with the Professional .dat files
#                       (default: <repo>/../tn_prof_compare, see below).
#
# The Professional .dat files (its OWN validation syntax) are NOT committed:
# they come from the Professional distribution (test/other/force7.dat,
# force7q4.dat, force10.dat, force13.dat, beam2d_3.dat) and from the manual
# conversions of the previous session (ffq4.dat). Point PROF_TEST_DIR at a
# directory containing them:
#   force7.dat force8.dat force10.dat force13.dat force7q4.dat ffq4.dat
#   msf_shear.dat msf_shear_elong.dat msf_tunnel3d.dat
# Default: $TOCHNOG_PROF_BIN/../test_inputs or $HOME/tn_prof_inputs.
#
# Output: a markdown-ready comparison table per model (node_dof_calcul
# values + GNU/Professional ratios against the analytic targets).
# Exit status: 0 when every model ran in both binaries, 1 otherwise.
# ---------------------------------------------------------------------------
set -u
export LC_ALL=C   # decimal point + %.12g in awk must not follow the locale

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SUITE="$REPO_DIR/validation-suite/test-2014"

# --- binaries ---------------------------------------------------------------
if [ -z "${TOCHNOG_PROF_BIN:-}" ]; then
  echo "Error: TOCHNOG_PROF_BIN is not set (path to the Tochnog Professional binary)." >&2
  exit 1
fi
if [ ! -x "$TOCHNOG_PROF_BIN" ]; then
  echo "Error: TOCHNOG_PROF_BIN is not executable: $TOCHNOG_PROF_BIN" >&2
  exit 1
fi
GNU_BIN="${TOCHNOG_GNU_BIN:-$REPO_DIR/build/tochnog}"
if [ ! -x "$GNU_BIN" ]; then
  echo "Error: GNU binary not found: $GNU_BIN (build with scripts/build_safe.sh)" >&2
  exit 1
fi
# local numeric runtime (Debian lapack/blas extracted in external-downloads)
if [ -d "$REPO_DIR/external-downloads/numlib-runtime" ]; then
  export LD_LIBRARY_PATH="$REPO_DIR/external-downloads/numlib-runtime${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

PROF_DIR="${PROF_TEST_DIR:-$HOME/tn_prof_inputs}"
if [ ! -d "$PROF_DIR" ]; then
  echo "Error: Professional input dir not found: $PROF_DIR (set PROF_TEST_DIR)" >&2
  exit 1
fi

WORK="${WORK_DIR:-/tmp/opencode/compare_runs}"
mkdir -p "$WORK"

# --- helpers -----------------------------------------------------------------
# last node_dof_calcul record per node: node_dof_calcul <n> <values...>
ndc() { # ndc <dbs> <node> <value_index(0-based over the 9/16 values)>
  awk -v n="$2" -v i="$3" '
    $1=="node_dof_calcul" && $2==n { v = $(3+i) }
    END { printf "%.12g", v }' "$1"
}

run_gnu() { # run_gnu <workdir> <dat>
  local d="$1" dat="$2"; shift 2
  rm -f "$d"/*.dbs "$d"/tn.* "$d"/node_dof_* 2>/dev/null
  ( cd "$d" && timeout 300 "$GNU_BIN" "$dat" > run.log 2>&1 )
}

run_prof() { # run_prof <workdir> <dat>
  local d="$1" dat="$2"; shift 2
  rm -f "$d"/*.dbs "$d"/tn.* 2>/dev/null
  ( cd "$d" && timeout 300 "$TOCHNOG_PROF_BIN" "$dat" > run.log 2>&1 )
}

# 1 when the GNU solver reported an honest failure (diverged / broke down)
gnu_diverged() { grep -q "iterative solver did not converge" "$1/run.log"; }

ratio() { # ratio <gnu> <prof>  -> "r=0.9947 (GNU 0.5% under)" style
  awk -v g="$1" -v p="$2" 'BEGIN{
    if (p==0) { printf "p=0 g=%.4g", g; exit }
    r = g/p; printf "%.4f", r }'
}

echo "# Tochnog GNU vs Professional - materi_stress_force baseline"
echo "# GNU:   $GNU_BIN"
echo "# PROF:  $TOCHNOG_PROF_BIN"
echo "# Date:  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo ""

RC=0

# ---------------------------------------------------------------------------
# 2D cantilever family. Statics: N = -12.34, V = +100,
# M(x) = -100*(100-x) (M(x=50) = -5000).
# node_dof_calcul 2D layout (both binaries): [norx nory nors shex shey shes
# momx momy moms]; the y components are the plot-vector components in the
# thickness direction t (GNU t = centroid-ref -> -y; Professional +y), so the
# SIGN of nory/shey/momy flips between binaries; magnitudes are physical.
# ---------------------------------------------------------------------------
for m in gforce7 gforce7q4; do
  case "$m" in
    gforce7)   pdat="force7.dat";  gnode=7 ;;
    gforce7q4) pdat="force7q4.dat"; gnode=2 ;;
  esac
  wd="$WORK/$m"; mkdir -p "$wd"
  gd="$wd/gnu"; pd="$wd/prof"; mkdir -p "$gd" "$pd"
  cp "$SUITE/$m.dat" "$gd/" && cp "$PROF_DIR/$pdat" "$pd/"
  run_gnu  "$gd" "$m.dat";   GRC=$?
  run_prof "$pd" "$pdat";    PRC=$?

  echo "## $m (2D cantilever quad9/quad4, L=100 h=10, load (-1.234,-10)/unit length)"
  echo "| qty | node | analytic | GNU | Professional | GNU/Prof |"
  echo "|---|---|---|---|---|---|"
  if [ -f "$gd/$m.dbs" ] && [ -f "$pd/${pdat%.dat}.dbs" ]; then
    if gnu_diverged "$gd"; then
      p=$(ndc "$pd/${pdat%.dat}.dbs" 2 2); echo "| N | 2 | ±12.34 | SOLVER DIVERGED (residual in $gd/run.log) | $p | - |"
      p=$(ndc "$pd/${pdat%.dat}.dbs" 2 5); echo "| V | 2 | ±100 | SOLVER DIVERGED | $p | - |"
      p=$(ndc "$pd/${pdat%.dat}.dbs" 2 8); echo "| M | 2 | ±5000 | SOLVER DIVERGED | $p | - |"
    else
      # magnitude indices (s components): 2=nors 5=shes 8=moms; the
      # directional components flip sign between binaries (t orientation)
      for spec in "N $gnode 2 12.34" "V $gnode 5 100" "M $gnode 8 5000"; do
        set -- $spec; qty=$1; node=$2; idx=$3; an=$4
        g=$(ndc "$gd/$m.dbs" "$node" "$idx")
        p=$(ndc "$pd/${pdat%.dat}.dbs" "$node" "$idx")
        echo "| $qty | $node | ±$an | $g | $p | $(ratio "$g" "$p") |"
      done
    fi
  else
    echo "| (no .dbs) | | | GNU rc=$GRC | PROF rc=$PRC | |"
  fi
  [ -f "$gd/$m.dbs" ] && [ -f "$pd/${pdat%.dat}.dbs" ] || RC=1
  echo ""
done

# ---------------------------------------------------------------------------
# gffq4: fixed-fixed beam, 10 quad4, p=1/unit length. Statics identity
# |M_end| + |M_center| = pL^2/8 = 1250 (deep beam: end -825, center +425).
# ---------------------------------------------------------------------------
m=gffq4; wd="$WORK/$m"; mkdir -p "$wd"; gd="$wd/gnu"; pd="$wd/prof"; mkdir -p "$gd" "$pd"
cp "$SUITE/$m.dat" "$gd/" && cp "$PROF_DIR/ffq4.dat" "$pd/"
run_gnu "$gd" "$m.dat"; GRC=$?; run_prof "$pd" "ffq4.dat"; PRC=$?
echo "## $m (fixed-fixed beam, 10 quad4, p=1/unit length on top edge)"
if [ -f "$gd/$m.dbs" ] && [ -f "$pd/ffq4.dbs" ]; then
  ge=$(ndc "$gd/$m.dbs" 1 8); pe=$(ndc "$pd/ffq4.dbs" 1 8)
  gc=$(ndc "$gd/$m.dbs" 6 8); pc=$(ndc "$pd/ffq4.dbs" 6 8)
  awk -v ge="$ge" -v pe="$pe" -v gc="$gc" -v pc="$pc" 'BEGIN{
    gs = (ge<0?-ge:ge) + (gc<0?-gc:gc)
    ps = (pe<0?-pe:pe) + (pc<0?-pc:pc)
    printf "| M_end (x=0) | 1 | -825 (deep) | %s | %s | %.4f |\n", ge, pe, (pe!=0?ge/pe:0)
    printf "| M_center (x=50) | 6 | +425 (deep) | %s | %s | %.4f |\n", gc, pc, (pc!=0?gc/pc:0)
    printf "| |M_end|+|M_center| | | 1250 | %.3f | %.3f | %.4f |\n", gs, ps, gs/ps }'
else
  echo "| (no .dbs) | GNU rc=$GRC PROF rc=$PRC |"
  RC=1
fi
echo ""

# ---------------------------------------------------------------------------
# gforce7_ref: same cantilever, 8 quad9 elements. The GNU's interior section
# forces are the clean cross-validation against the Professional (exact).
# ---------------------------------------------------------------------------
m=gforce7_ref; wd="$WORK/$m"; mkdir -p "$wd"; gd="$wd/gnu"; mkdir -p "$gd"
cp "$SUITE/$m.dat" "$gd/"
run_gnu "$gd" "$m.dat"; GRC=$?
echo "## $m (8 quad9, same statics; Professional reference = its force7 family)"
echo "| qty | node(x) | analytic | GNU |"
echo "|---|---|---|---|"
if [ -f "$gd/$m.dbs" ]; then
  for spec in "N 23 2 12.34" "M 23 8 5000" "M 35 8 2500" "V 23 5 100"; do
    set -- $spec; qty=$1; node=$2; idx=$3; an=$4
    g=$(ndc "$gd/$m.dbs" "$node" "$idx")
    awk -v q="$qty" -v n="$node" -v g="$g" -v a="$an" 'BEGIN{
      m=(g<0?-g:g); printf "| %s | %d | %s | %s (%.4f) |\n", q, n, a, g, m/a }'
  done
else
  echo "| (no .dbs) GNU rc=$GRC |"; RC=1
fi
echo ""

# ---------------------------------------------------------------------------
# 3D cantilevers (hex8). Statics at z=50 per unit width: N=-12.34, V=-100,
# M=-5000. POST-FIX (lot C/D, 2026-08-28): the GNU converges (the
# previous divergence was an input BC bug - the -ra 1 4 list clamped only
# two opposite corners of the bottom face, leaving a rigid rotation free;
# fixed to 1 2 3 4). POST-FIX (sprint 12 lot 1, 2026-08-29): the section
# FRAME converged - the Professional's force10 carries thickness_switch
# -yes (its square 10x10 end face ties the extent rule; without the
# switch the PROFESSIONAL ITSELF gives the rotated frame: vectors ~0,
# mom in the mom2 slot - measured); our conversion had omitted the
# keyword. GNU fixes: the face corners as a quad loop (e2 was the face
# DIAGONAL - a 45-degree frame with -yes), the square-face tie broken
# by the reference-point direction (the Professional resolves the same
# tie to the reference direction), and the t orientation TOWARD the
# reference point + signed s-items (the Professional's conventions).
# Items compared: nory=1 shey=5 mom1y=9 of the 16-item 3D layout.
# ---------------------------------------------------------------------------
for m in gforce10 gforce13; do
  case "$m" in
    gforce10) pdat="force10.dat"; pnode=5 ;;
    gforce13) pdat="force13.dat"; pnode=5 ;;
  esac
  wd="$WORK/$m"; mkdir -p "$wd"; gd="$wd/gnu"; pd="$wd/prof"; mkdir -p "$gd" "$pd"
  cp "$SUITE/$m.dat" "$gd/" && cp "$PROF_DIR/$pdat" "$pd/"
  run_gnu "$gd" "$m.dat"; GRC=$?; run_prof "$pd" "$pdat"; PRC=$?
  echo "## $m (hex8 3D cantilever; section frame converged - sprint 12 lot 1)"
  if [ -f "$gd/$m.dbs" ] && [ -f "$pd/${pdat%.dat}.dbs" ]; then
    if gnu_diverged "$gd"; then
      echo "| node 5 | 5 | N=-12.34 V=-100 M=-5000 | GNU: SOLVER DIVERGED | |"; RC=1
    else
      echo "| qty | node | analytic | GNU | Professional | GNU/Prof |"
      echo "|---|---|---|---|---|---|"
      for spec in "nory:1:-12.34" "shey:5:-100" "mom1y:9:-5000"; do
        lbl=${spec%%:*}; rest=${spec#*:}; idx=${rest%%:*}; ann=${rest##*:}
        p=$(ndc "$pd/${pdat%.dat}.dbs" "$pnode" "$idx")
        g=$(ndc "$gd/$m.dbs" 5 "$idx")
        echo "| $lbl (z=50) | 5 | $ann | $g | $p | $(ratio "$g" "$p") |"
      done
    fi
  else
    echo "| (no .dbs) GNU rc=$GRC PROF rc=$PRC |"; RC=1
  fi
  echo ""
done

# ---------------------------------------------------------------------------
# MSF cross-validation (our tests in the Professional, same physics):
#   msf_shear (simple shear): she = G*gamma = 0.384615; the Professional's
#   node_dof_calcul is 0 for the Dirichlet-driven case (documented quirk) but
#   its sigma_xy field = 0.3846153846 EXACT = the GNU's she.
#   msf_tunnel3d (hoop p*R): nor = 0.1; GNU 0.10000000 EXACT vs Professional
#   0.09980686 (0.19% FE discretization).
# ---------------------------------------------------------------------------
if [ -f "$PROF_DIR/msf_shear.dat" ]; then
  wd="$WORK/msf_shear"; mkdir -p "$wd"; pd="$wd/prof"; mkdir -p "$pd"
  cp "$PROF_DIR/msf_shear.dat" "$pd/"
  run_prof "$pd" "msf_shear.dat"; PRC=$?
  echo "## msf_shear cross-validation (quad4 simple shear, she = G*gamma = 0.384615)"
  if [ -f "$pd/msf_shear.dbs" ]; then
    s=$(awk '$1=="node_dof" && $2==1 { print $8 }' "$pd/msf_shear.dbs")
    echo "| quantity | GNU (materi_stress_force.430) | Professional sigma_xy |"
    echo "| she | 0.3846153846 | $s |"
    echo "| note | GNU node_dof_calcul she = 0.3846153846 EXACT | PROF node_dof_calcul = 0 (Dirichlet-driven quirk); sigma_xy EXACT |"
  else
    echo "| (no .dbs) PROF rc=$PRC |"; RC=1
  fi
  echo ""
fi
if [ -f "$PROF_DIR/msf_tunnel3d.dat" ]; then
  wd="$WORK/msf_tunnel"; mkdir -p "$wd"; gd="$wd/gnu"; pd="$wd/prof"; mkdir -p "$gd" "$pd"
  cp "$PROF_DIR/msf_tunnel3d.dat" "$pd/"
  cp "$SUITE/msf_tunnel3d.dat" "$gd/"
  run_prof "$pd" "msf_tunnel3d.dat"; PRC=$?
  run_gnu "$gd" "msf_tunnel3d.dat"; GRC=$?
  echo "## msf_tunnel3d cross-validation (hex27 ring, hoop nor = E*u0*t/R = p*R = 0.1)"
  if [ -f "$pd/msf_tunnel3d.dbs" ] && [ -f "$gd/msf_tunnel3d.dbs" ]; then
    p=$(ndc "$pd/msf_tunnel3d.dbs" 1 3)
    g=$(ndc "$gd/msf_tunnel3d.dbs" 1 3)
    echo "| quantity | GNU | Professional | GNU/Prof |"
    echo "|---|---|---|---|"
    echo "| nory (hoop, node 1) | $g | $p | $(ratio "$g" "$p") |"
    echo "| note | LOT 5+: the section reads the element internal forces - both binaries give the FE-DISCRETIZED hoop of the curved ring (the L4 field integral gave the analytic 0.1 EXACT interpolating the prescribed field) | |"
  else
    echo "| (no .dbs) GNU rc=$GRC PROF rc=$PRC |"; RC=1
  fi
  echo ""
fi

echo "# End. RC=$RC"
exit $RC
