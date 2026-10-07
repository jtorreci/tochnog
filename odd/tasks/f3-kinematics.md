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
- [ ] T1 Read-only design investigation: measure Pro rules (sweep), map GNU code paths, propose one design with corpus impact. Route: delegated read-only explorer.
- [ ] T2 Implementation per the agreed design (after user decision).

## Rules
As F1/F2 (manual/analytic decides; Pro oracle when the manual is silent;
divergences need an owner decision and go to DIVERGENCES-PROFESSIONAL.md).

## Progress
- Branch `dev/f3-kinematics` from `c523094`.
