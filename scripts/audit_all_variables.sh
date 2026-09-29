#!/bin/bash
# ---------------------------------------------------------------------------
# audit_all_variables.sh - run every test on the GNU build AND on Tochnog
# Professional, compare EVERY output variable of the final state, evaluate
# each test's own targets on both dumps, and triage the differences.
#
# Report: ProjectDocs/AUDIT-ALL-VARIABLES-2026.md (method, clusters, verdicts).
#
# Usage:
#   TOCHNOG_PROF_BIN=/path/to/professional/src/tochnog \
#     scripts/audit_all_variables.sh [OUT_DIR]
#
# Environment:
#   TOCHNOG_PROF_BIN  Professional binary (REQUIRED).
#   TOCHNOG_GNU_BIN   GNU binary (default <repo>/build/tochnog; build it clean
#                     first: scripts/build_safe.sh --clean).
#   PROF_TEST_DIR     Professional test tree (default <prof>/../../test).
#   AUDIT_TIMEOUT     per test per binary, seconds (default 45).
#   AUDIT_LONG        comma list of Pro 'other'/'validation' tests re-run at
#                     AUDIT_LONG_TIMEOUT (default: the 2026 list below; set
#                     to "" to skip).
#   AUDIT_JOBS        parallel runs per binary (default 7; both binaries run
#                     at the same time).
# Outputs (OUT_DIR, default /tmp/tochnog-audit): runs/<corpus>/{gnu,gnuconv,
# pro}/<test>/ (inputs, logs, .dbs), runs/summary.json, runs/clusters.json,
# runs_iter/ (5-iteration diagnostic), report_tables.md.
# Nothing in the repository is written. LC_ALL=C throughout.
# ---------------------------------------------------------------------------
set -u
export LC_ALL=C
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/.." && pwd)"
PY="$SCRIPT_DIR/audit_all_variables.py"

PRO="${TOCHNOG_PROF_BIN:?set TOCHNOG_PROF_BIN to the Professional binary}"
GNU_SRC="${TOCHNOG_GNU_BIN:-$REPO/build/tochnog}"
PT="${PROF_TEST_DIR:-$(cd "$(dirname "$PRO")/.." && pwd)/test}"
OUT="${1:-/tmp/tochnog-audit}"
TO="${AUDIT_TIMEOUT:-45}"
J="${AUDIT_JOBS:-7}"
LONG_TO="${AUDIT_LONG_TIMEOUT:-180}"
LONG_DEFAULT="foundation1,dam_building,slope_classical_numerical,vajont_slope_fs,hypo7,hypo8,hypo9,diprisc1,ground11_nonsaturated,ground11_phreatic_level,mohr_coul_direct4,mohr_coul_direct5,large2,force9,force11,tensile_instability,dynamic8,hypo12,post7,slide1,slide4,validation_9,validation_14,validation_12,validation_8"
LONG="${AUDIT_LONG-$LONG_DEFAULT}"
LD="$REPO/external-downloads/numlib-runtime"

[ -x "$PRO" ] || { echo "Professional binary not executable: $PRO" >&2; exit 1; }
[ -x "$GNU_SRC" ] || { echo "GNU binary not found: $GNU_SRC" >&2; exit 1; }
mkdir -p "$OUT/bin" "$OUT/runs" "$OUT/runs_iter"
# freeze the GNU binary: a rebuild during the audit must not mix binaries
GNU="$OUT/bin/tochnog_gnu"; cp "$GNU_SRC" "$GNU"
echo "# GNU  $(sha256sum "$GNU" | cut -c1-16)  $(git -C "$REPO" rev-parse --short HEAD)"
echo "# PRO  $(sha256sum "$PRO" | cut -c1-16)  $PRO"

run() { python3 "$PY" run "$@"; }
both() { # both <corpus> <src> <convert> [extra run args...]
  local c=$1 src=$2 cv=$3; shift 3
  if [ -n "$cv" ]; then
    run --bin "$PRO" --convert "$cv" --src "$src" --out "$OUT/runs/$c/pro" --jobs "$J" --timeout "$TO" "$@" &
    run --bin "$GNU" --ld "$LD" --src "$src" --out "$OUT/runs/$c/gnu" --jobs "$J" --timeout "$TO" "$@" &
    wait
    run --bin "$GNU" --ld "$LD" --convert "$cv" --src "$src" --out "$OUT/runs/$c/gnuconv" --jobs "$((2*J))" --timeout "$TO" "$@"
  else
    run --bin "$PRO" --src "$src" --out "$OUT/runs/$c/pro" --jobs "$J" --timeout "$TO" "$@" &
    run --bin "$GNU" --ld "$LD" --src "$src" --out "$OUT/runs/$c/gnu" --jobs "$J" --timeout "$TO" "$@" &
    wait
  fi
}

# 1. GNU-dialect corpora: GNU on the original, GNU and Pro on the converted text
both sfnet "$REPO/external-downloads/sfnet/extracted/test" pro
both vs    "$REPO/validation-suite/test-2014" pro --exclude print_mesh_dof
both t26   "$REPO/validation-suite/test-2026" pro
# 2. Professional distribution tests, as distributed, on both binaries
for sub in other validation incremental_driver tutorial/tutorial_1 tutorial/tutorial_2 \
           tutorial/tutorial_3 tutorial/tutorial_4 easy_tochnog/axisymmetric_pile \
           easy_tochnog/dam easy_tochnog/easy easy_tochnog/shallow_foundation \
           easy_tochnog/sheet_pile easy_tochnog/slope; do
  ls "$PT/$sub"/*.dat >/dev/null 2>&1 && both "pro_${sub//\//_}" "$PT/$sub" ""
done
# 3. a handful of long Professional tests at a longer timeout
if [ -n "$LONG" ]; then
  for c in other validation; do
    TO=$LONG_TO both "pro_$c" "$PT/$c" "" --only "$LONG"
  done
fi
# 4. compare, then the 5-iteration diagnostic on the mismatching tests only
python3 "$PY" aggregate "$OUT/runs"
declare -A SRC=( [sfnet]="$REPO/external-downloads/sfnet/extracted/test:pro+iter"
                 [vs]="$REPO/validation-suite/test-2014:pro+iter"
                 [t26]="$REPO/validation-suite/test-2026:pro+iter"
                 [pro_other]="$PT/other:iter" [pro_validation]="$PT/validation:iter" )
for c in "${!SRC[@]}"; do
  only=$(python3 "$PY" mism-list "$OUT/runs" "$c")
  [ -z "$only" ] && continue
  src=${SRC[$c]%:*}; cv=${SRC[$c]##*:}
  run --bin "$PRO" --convert "$cv" --src "$src" --out "$OUT/runs_iter/$c/pro" --jobs "$J" --timeout "$TO" --only "$only" &
  run --bin "$GNU" --ld "$LD" --convert "$cv" --src "$src" --out "$OUT/runs_iter/$c/gnu" --jobs "$J" --timeout "$TO" --only "$only" &
  wait
done
python3 "$PY" aggregate "$OUT/runs_iter"
python3 "$PY" iter-effect "$OUT/runs" "$OUT/runs_iter" "$OUT/iter_effect.json"
# 5. triage clusters and the report tables
python3 "$PY" clusters "$OUT/runs" "$OUT/iter_effect.json"
python3 "$PY" stats "$OUT/runs" > "$OUT/report_tables.md"
echo "tables: $OUT/report_tables.md   per-test data: $OUT/runs/summary.json"
