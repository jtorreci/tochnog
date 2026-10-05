# F0 — test infrastructure hygiene

## Objective
Make the regression safety net trustworthy before the next bug-fix phases.

## Tasks
- [x] T1 Version own `validation-suite/test-2026/*.dat` inputs (outputs stay ignored; `.git/info/exclude` adjusted locally). Route: inline. Commit `b147ee4`.
- [x] T2 Explore `scripts/build_safe.sh` file-check failures (~75 pre-existing): classify stale / broken check / real regression; propose cleanup. Route: delegated read-only explorer.
- [x] T3 Uninitialized `element_dof` denormals (~1e-310) in gen1, genbeam1, iface_gen*, conspr3, genera1/2, ho_othr1, spring2/3/5: find the source. Route: delegated read-only explorer.
- [x] T4 Apply the agreed fixes for T2/T3. Route: delegated direct (one writer; 2+ non-trivial files).
  - [x] T4.1 `generate.cc`: `array_set( tmp_element_dof, 0, mnolnuknwn )` after the 4 `get_new_dbl(mnolnuknwn)` (generate_spring / _beam_truss / _interface / _interface_triangle), as `create.cc`. Commit `a32ea74`. Evidence: clean builds of `b147ee4` (old) and the fix (new) in scratch; 2 runs of gen1, genbeam1, iface_gen* (6), conspr3, genera1/2, ho_othr1, spring2/3/5 byte-identical, 0 `e-3xx` (old: 1-574 lines/file); A/B diff = `element_dof` lines only. Docs `manual-*/generated_element_dof_init.md`, SEGUIMIENTO row, BACKLOG C23 (unrelated `dofid_no` `dof.21` garbage found by the A/B).
  - [x] T4.2 Version test-2014 inputs (`.gitignore` + local `.git/info/exclude`; `print_mesh_dof.dat` stays ignored: it is the meshdoff output dump, and the corpus "1 fail" of 233/1). Commit `2b6938c` (233 `.dat`, state as on disk incl. C17).
  - [x] T4.3 `scripts/build_safe.sh` repair: run list in `scripts/build_safe_tests.txt` (226 runs = pre-incident loop - reset1 - materi_direct_auto + gforce7 + tsup_*), reduced-mode switch removed, one tsup set (`tsup_check`; gap/tcap made non-vacuous), duplicate full-mode tsup blocks with broken `rg -a` removed, run-rc gate, `TN_CHECK_ONLY=1`, check counters. Commit `365ec1e` (with T4.4). Negative proof: a bogus list entry -> `runs fallidos: 1/227`, rc=1; `TN_CHECK_ONLY=1 TN_BIN=nope` -> error.
  - [x] T4.4 G5: `dpline_move.dat` line 40 coords `0.5 0.5 0.5 0.6`; `dpline_method.dat` nodes 3/4 at y=1 + `node_start_refined` 1-4. Both checks OK (vely.38 2 lines y2=0.995 vely=-0.00995, vely.39 1 line; disy.40 line 3 y=0.499 vs disy.41 y=0.5). Commit `365ec1e`.
  - [x] T4.5 `TN_CHECK_ONLY=1 ./scripts/build_safe.sh` on the new `build/tochnog`: runs 226/226 OK, checks 92 OK / 9 FALLO = exactly the G2 set below. rc=1 (expected until T5).
  - [x] T4.6 Corpora old vs new (clean builds): validation-suite test-2014 233/1 = 233/1, test-2026 47/0 = 47/0, sfnet 187/10 = 187/10 (same per-test status); changed outputs: `element_dof`-only, plus wall-clock `time_calculation` and the pre-existing `dof.21` (C23).
- [x] T5 Regenerate the inputs lost in the rm incident (G2): `msf_errors_nogroup`, `msf_errors_bothdir`, `msf_errors_3dnodir`, `msf_errors_3dnoref`, `msf_errors_refcount`, `msf_errors_posttype`, `msf_errors_3dbothdir`, `msf_errors_3drefcount`, `tslv_zip`. Route: inline (9 small mechanical inputs; intent from build_safe comments, SEGUIMIENTO lote-1/Sprint-12 rows, `calcul_force.cc` validate/`post_calcul_materi_stress_force`, `miscel.cc` zip). Each msf input = msf_errors_2dwarn (2D quad4) / msf_parse_3d (3D hex8) base with ONE defect, plus an `EXPECTED_ERROR:` marker for corpus classifiers. Evidence: each rc=1 with the intended message (no parse error); positive controls (same input with the defect removed) all rc=0; `tslv_zip` isolated run -> `tslv_zip.dbs.gz`, "zip -yes: output files zipped with gzip"; `TN_CHECK_ONLY=1 ./scripts/build_safe.sh` -> checks 101 OK / 0 FALLO, runs 226/226, rc=0. Commit `071af0c`.
- [x] T6 BACKLOG C23 (`dofid_no` `dof.21` garbage). Route: inline fix (one file, `print_hi.cc`) + delegated read-only explorer for the 192-site `get_new_dbl`/`new double[` audit. Root cause: NOT an uninitialized buffer; `print_dof()` looped `ipuknwn<nuknwn_` (value+derivative entries), repeated vector/stress blocks and read `node_dof` past the record end (heap garbage). Fix: primary-unknown loop as `print_vt.cc`. `build_safe.sh`: dof.21 check asserts 48 lines + no `e-3xx`; `freq_timeint` 176 -> 40 lines/print. Evidence: clean scratch build (33 s, rc 0) copied to `build/tochnog`; `dof.21` 768 lines/20 garbage -> 48, two runs byte-identical, 0 `e-3xx`; `TN_CHECK_ONLY=1` -> runs 226/226, checks 101 OK / 0 FALLO, rc=0; corpus A/B old (`a32ea74` build) vs new: same per-test status (test-2014 242 run + tslv_zip isolated: 241/1 = `print_mesh_dof.dat`; test-2026 47/0; sfnet 187/10), changed outputs only `dof.*` and wall-clock `time_calculation`. Audit: no confirmed a32ea74-pattern site; follow-ups listed as BACKLOG C25; same loop in `print_fr.cc`/`print_gm.cc` = C24. Docs `manual-*/print_dof_component_loop.md`, SEGUIMIENTO row, C23 closed. Commit: recorded in the next work unit.

## Notes
- `validation-suite/test-2014` inputs are versioned since `2b6938c`
  (including the C17 `materi_direct_mc` constraint fix).
- `build/tochnog` = clean build of `a32ea74` sources (generate fix).
- Stale doc outside T4 scope: `ProjectDocs/QUALITY-CI.md` still describes
  the 16-test reduced mode.
- Most test-2014 inputs were regenerated from `.dbs` (sorted records,
  `node_dof` lines); only dpline_move/dpline_method showed post-run state
  in a scan (non-zero `node_dof` in 23 files look like initial states).
- Branch `dev/f0-version-tests` from `8f6fbe3`.
