#!/bin/bash
# CI smoke runner for tochnog.
#
# The full validation suite (validation-suite/test-2014, ~16 recovered
# tests run by scripts/build_safe.sh) is NOT versioned in this repo
# (see .gitignore), so CI runs this small committed subset instead:
# parse + mesh + elastic materi + solve + database dump, which is the
# core pipeline of every regression.
#
# Usage:
#   .github/ci/run_smoke.sh <path-to-tochnog-binary>
#
# Optional environment:
#   ASAN_OPTIONS / UBSAN_OPTIONS   passed through to the runs (sanitize job).
#
# Exit code 0 only if every input runs rc=0, writes <name>.dbs and the
# .dbs contains solved nodal dofs plus the end_data marker.
set -uo pipefail

BIN="${1:?usage: run_smoke.sh <path-to-tochnog-binary>}"
REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
case "$BIN" in
  /*) BIN_ABS="$BIN" ;;
  *)  BIN_ABS="$REPO_DIR/$BIN" ;;
esac
SMOKE_DIR="$REPO_DIR/.github/ci/smoke"
WORK="$(mktemp -d /tmp/tn_smoke.XXXXXX)"
trap 'rm -rf "$WORK"' EXIT

export LC_ALL=C
FAILED=0
TOTAL=0

for dat in "$SMOKE_DIR"/*.dat; do
  name="$(basename "$dat" .dat)"
  TOTAL=$((TOTAL + 1))
  cp "$dat" "$WORK/$name.dat"
  ( cd "$WORK" && timeout 300 "$BIN_ABS" "$name.dat" > "$WORK/$name.run.out" 2>&1 )
  rc=$?

  ok=1
  if [ "$rc" -ne 0 ]; then
    echo "FAIL $name: exit rc=$rc"
    tail -20 "$WORK/$name.run.out"
    ok=0
  elif [ ! -f "$WORK/$name.dbs" ]; then
    echo "FAIL $name: no $name.dbs written"
    ok=0
  elif ! grep -q "end_data" "$WORK/$name.dbs"; then
    echo "FAIL $name: $name.dbs has no end_data marker"
    ok=0
  elif ! grep -q "^node_dof" "$WORK/$name.dbs"; then
    echo "FAIL $name: $name.dbs has no solved node_dof records"
    ok=0
  fi

  if [ "$ok" -eq 1 ]; then
    echo "PASS $name (rc=0, .dbs with solved dofs)"
  else
    FAILED=$((FAILED + 1))
  fi
done

echo "=== smoke: $((TOTAL - FAILED))/$TOTAL PASS ==="
[ "$FAILED" -eq 0 ]
