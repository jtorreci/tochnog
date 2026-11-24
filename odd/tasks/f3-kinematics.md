# F3 — kinematics defaults and memories (C8, C18, C19, C20)

## Objective
Converge the GNU kinematic defaults and memory formulations with the
Professional with one coherent design, because they all change the same
corpus results.

## Items
- C8: default `mesh` `-follow_material` (GNU) vs `-fixed_in_space` (Pro).
- C18: default memory with `materi_displacement` behaves as `-total_linear` in Pro (C16 column: GNU node-3 disx 0.472 vs Pro 0.5; elem.cc:569 coordinate update; A2 rejected because it broke mohr_coul_direct7 and sfnet temp5).
- C19: explicit `-updated_linear` on `-follow_material`: Pro current-geometry gradient (0.69565) vs GNU reference gradient (1.0).
- C20: `-updated_jaumann` (Pro manual 6.685/6.686) not implemented.

## Tasks
- [x] T1 Read-only design investigation: measure Pro rules (sweep), map GNU code paths, propose one design with corpus impact. Route: delegated read-only explorer. Outcome: design A, owner decision 2026-10-07.
- [x] T2 Implementation per design A. Route: delegated direct (one writer).
  - [x] WU1 `248ce83` (D1 mesh default `-fixed_in_space`, D2 memory default
    `-total_linear` with `materi_displacement`, D3 reference gradient of
    `-updated_linear` removed; C8, C18, C19).
  - [x] WU2 `09f0948` (D4 `-updated_jaumann`; C20).
  - [x] `7abd713` docs(agents): cap memory of every tochnog run (two
    session crashes by the global OOM killer, uncapped corpus runs).

## Rules
As F1/F2 (manual/analytic decides; Pro oracle when the manual is silent;
divergences need an owner decision and go to DIVERGENCES-PROFESSIONAL.md).

## Progress
- Branch `dev/f3-kinematics` from `c523094`.
- 2026-10-08 T2 done. Evidence:
  - Tests (test-2026, RED old rc=1 / GREEN new rc=0 / Pro rc=0):
    `mesh_default_fixed_in_space1` bar sigxx 0.407136 -> 0.5 (Pro 0.5);
    `mesh_default_fixed_in_space2` rotation 0.109567 -> 0 (Pro 0);
    `materi_memory_default_displacement1` column disx 0.472222 -> 0.5 (Pro 0.5);
    `materi_memory_updated_linear_current1` bar 0.5 -> 0.4071364 (Pro 0.4071364);
    `materi_memory_updated_linear_current2` column 100 steps 0.472222 -> 0.520480
    (Pro 0.520931, tol 1e-3, in-step geometry divergence);
    `materi_memory_updated_jaumann_shear1/2` = Pro exactly
    (0.125/0.234375/-0.125; 0.0618335/0.2396333/-0.0618335).
  - Clean build (git archive of `09f0948`), `build/tochnog` updated;
    `TN_CHECK_ONLY=1 ./scripts/build_safe.sh`: 138 OK / 0 FALLO, 226/226 runs.
  - Corpus A/B exit code, capped (ulimit -v 4 GB, 60 s, 4 jobs), old -> new:
    Pro other 239 -> 243 (+axisym9, distri3, hypo12, relax2, update3;
    -mohr_coul_direct7 = BACKLOG C41); sfnet 179 -> 171 (-axisym3, axisym4,
    contac1, contac2, examp12, temp4, temp5, visc_pl2 = BACKLOG B13);
    validation-suite 233 = 233; test-2026 74 -> 81 (the 7 new tests).
  - No corpus input was edited (no test-2014/test-2026 target changed).
  - Memory-heavy inputs (hit the 4 GB cap with BOTH binaries, pre-existing):
    Pro other `force11`, `large1`, `ground19_water_under_dam`,
    `mesh_activate1`.
  - Divergences added: F3-INSTEP, F3-CONDIF, F3-ULDISP
    (DIVERGENCES-PROFESSIONAL.md); C8/C18/C19/C20 rows removed (closed).
- Next step: owner review; push/merge are owner decisions.
