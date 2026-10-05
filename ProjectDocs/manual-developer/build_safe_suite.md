# build_safe_suite

The internal regression suite of `scripts/build_safe.sh` (repaired in F0
T4, 2026-10-05). Developer tooling only: no user-visible feature.

## Implementation

- `scripts/build_safe.sh` builds `build/$TN_BIN` (sequential `make -j1`,
  memory limits) and then runs the suite over
  `validation-suite/test-2014` (inputs versioned since F0 T4).
- Run list: `scripts/build_safe_tests.txt`, one input name per line
  (`#` comments, blank lines ignored). Restored from the full
  pre-incident loop (`git show 2488c4a^:scripts/build_safe.sh`) minus
  `reset1` (input lost) and `materi_direct_auto` (expected-error test,
  checked by `expected_error_ok`), plus `gforce7` and the eleven
  `tsup_*` tests: 226 runs.
- Each run: `timeout 120` under `ulimit -v $TN_MEMLIMIT_KB`, stdout in
  `/tmp/<test>_safe.out`. A run with rc != 0 is recorded in
  `RUN_FAILED`.
- File checks: `check_ok` / `check_fail` (counted in `N_CHECK_OK` /
  `N_CHECK_FAIL`). Expected-error tests run outside the loop
  (`msf_error_ok`, `expected_error_ok`, `run_ok`).
- Gate at the end: the build fails when any check failed OR any run of
  the loop returned non-zero (`runs fallidos: N/M: <names>`). Summary
  line: `==> Checks: <ok> OK / <fail> FALLO; runs: <ok>/<total> OK`.
- Environment variables: `TN_EXTRA_FLAGS`, `TN_BIN`, `TN_MEMLIMIT_KB`,
  `TN_SKIP_SUITE` (build only) and `TN_CHECK_ONLY=1` (no build: run the
  suite on the existing `build/$TN_BIN`, error if it does not exist).

## Changes in F0 T4

- Removed the post-incident "reduced mode" (switch on the existence of
  `msf_beam2d.dat`): the full set of file checks always runs.
- Kept ONE set of `tsup_*` + `gforce7` checks (the former reduced-mode
  set, as `tsup_check` calls). The former full-mode duplicate was
  deleted: its targets described pre-incident inputs and its
  `rg -a "<record>" 1 0 <file>` calls never matched. `tsup_gap` and
  `tsup_tcap` now really check (the reduced versions exited on the
  first record / lacked an `END` exit).
- Restored `dpline_move.dat` (line 40 coordinates `0.5 0.5 0.5 0.6`) and
  `dpline_method.dat` (top nodes at y = 1, `node_start_refined` records):
  the regenerated copies held the post-run state (moved line, deformed
  mesh) and failed their checks.

## Result (2026-10-05, binary with the `element_dof` fix)

`TN_CHECK_ONLY=1`: runs 226/226 OK, checks 92 OK / 9 FALLO. The 9
failures are the inputs lost in the rm incident (not regenerated yet):
`msf_errors_{nogroup,bothdir,3dnodir,3dnoref,refcount,posttype,3dbothdir,3drefcount}`
and `tslv_zip`.

## Pending

- Regenerate the 9 lost inputs above (F0 T5).
- `ProjectDocs/QUALITY-CI.md` still describes the 16-test reduced mode.
- The tsup targets are the regenerated-input values; their derivation
  comments refer to the original Sprint 13 models.
