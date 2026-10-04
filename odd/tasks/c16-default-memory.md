# C16 — default material memory `-updated_linear` (as Professional)

## Objective
An unspecified `group_materi_memory` behaves as `-updated_linear` (Professional
default), and `control_materi_updated_apply -yes` switches unspecified memories
to `-updated` (manual 6.152).

## Problem (evidence: scratchpad c16/, c16b/; BACKLOG C16)
GNU hard-codes `memory=-UPDATED` as default at materi.cc (materi, set_deften_etc),
stress.cc, elem.cc and tendon.cc. Every explicit memory already matches Pro; only
the default differs (single quad simple shear: Pro linear, GNU Jaumann-like).

## Decisions
- 2026-10-04 user: adopt option A (measured: no GNU-corpus regressions; Pro
  targets 220 → 222, final-state matches 155 → 156, none lost).
- A2 (materi_displacement → `-total_linear`) rejected: breaks
  `mohr_coul_direct7` and sfnet `temp5`.
- 2026-10-04 user: `control_materi_updated_apply -yes` follows manual 6.152
  (option 2); the Pro binary ignores -yes → deliberate divergence (BACKLOG C21).

## Tasks
- [x] T1 Helper `materi_memory_default()` (general.cc) at the 5 sites, tendon assert, reference gradient only for explicit `-updated_linear` (elem.cc). Route: delegated (multi-file; reuse scratch source c16b/src). Commit `d8117aa`.
- [x] T2 Tests (simple shear default = Pro; `control_materi_updated_apply -yes` → `-updated`); corpus A/B; `hypo12` tolerance decision. Route: inline (single worker). Commit `0c7d4d3`.
- [x] T3 Docs (user + developer), SEGUIMIENTO, BACKLOG (close C16; new: default with materi_displacement behaves as `-total_linear` in Pro; explicit `-updated_linear` gradient current vs reference; missing `-updated_jaumann`). Commit `40dabbd`.

## Checks
Clean build; RED/GREEN; validation-suite 233/1, test-2026, sfnet 187/10 unchanged
except explained changes.

## Progress
- Branch `dev/c16-default-memory` from `9ebc367`.
- 2026-10-04 T1 applied in the working tree (uncommitted): option A from
  scratch c16b/src confirmed to be exactly HEAD 9ebc367 + A; tendon.cc edited
  in place (CRLF kept); stress.cc comment updated (-yes lives in
  materi_memory_default, -no gate stays in stress.cc). Clean scratch build OK.
- T2 evidence (uncommitted): materi_memory_default_shear1 (RED old rc=1,
  GREEN new rc=0, Pro 0/0.5/0) and materi_memory_default_stretch1 (old/new/Pro
  sigxx 0.695653; regression guard) in validation-suite/test-2026, registered
  in scripts/build_safe.sh. Corpus A/B new vs exp: identical except wall-clock
  noise; base vs new: hypo12 PASS->RUNFAIL, mohr_coul_direct6/8 RUNFAIL->PASS;
  16 load-induced timeouts rerun on base at 300 s: same statuses as new.
- BLOCKER (needs decision): Pro 02-08-2026 binary IGNORES
  `control_materi_updated_apply -yes` for an unspecified memory (simple shear
  stays 0/0.5/0 for dis, -follow_material, -fixed_in_space), contradicting
  manual 6.152 and the -yes branch of materi_memory_default(). Explicit
  -updated + -yes = -updated (both); -no turns explicit -updated linear (both).
  No corpus .dat uses -yes. Test (b) not written; no commits made.
- 2026-10-04 decision applied: test (b) `materi_memory_default_updated_apply1`
  targets the -updated values 0.231098 / 0.420739 / -0.228598 (GNU new rc=0;
  Pro misses them: 0 / 0.5 / 0). Not RED on the old binary (old default was
  already -updated); negative control: new binary without the -yes line rc=1.
- Commits: `d8117aa` feat (T1), `0c7d4d3` test registration (T2; .dat files
  gitignored, on disk), `40dabbd` docs (T3). New clean binary copied to
  build/tochnog. BACKLOG: C16 → R10; new C18-C22; 28 open (9 MEDIUM, 19 LOW).
- Review: RDD not run by this worker; pending parent decision.
- Next: parent review; push/PR are user decisions.
