# Corpus triage: last 5 sfnet RUNFAILs (parse/init/consolidation)

Date: 2026-09-14
Branch: `dev/parse-init-examp2` (worktree `/tmp/opencode/final5-20260914-160749`)
Baseline: `0d0ad0c` (`dev/masin-kernel-surface`), clean `build_safe.sh`.
Corpus: `external-downloads/sfnet/extracted/test` (197 `.dat`, **untracked**).
GNU baseline binary: `/tmp/opencode/final5-evidence-20260914-160940/baseline-tochnog`
(clean build from `0d0ad0c`).
Pro binary: `/tmp/opencode/tn_prof/tochnog_version_02-08-2026/src/tochnog`
(build 2-8-2026). jan-2014 reference GNU: `/tmp/opencode/tn2014/build/tochnog`.
All runs with `LC_ALL=C`.

## Result

| corpus | PASS | RUNFAIL | PARSE |
|---|---|---|---|
| baseline `0d0ad0c` | 186 | 11 | 0 |
| after `so.cc` band fix | **187** | **10** | **0** |

Per-test diff: **1 RUNFAIL -> PASS (`examp2`), 0 regressions**.
Suite `build_safe.sh`: **16/16**.

## Triage

| test | root cause | class | evidence |
|---|---|---|---|
| `examp2` | Bi-CG direct-LU retry band width short by `nprinc-1` -> out-of-bounds band fill -> wrong retry solution `-0.7206` | **FIXABLE (source, `so.cc`)** | LAPACK/SuperLU direct give the exact `-0.634519741323` = Pro `-0.6345197413352`; bisect -> `6cafafe`; instrumentation showed `band=2`, dof distance `3` |
| `pruebavel` | `group_condif_flow 0 -velx -vely` uses dof names; manual Pro 6.600 is numeric `beta1 beta2 beta3`; the validation-suite copy is 0 bytes | **OUT-OF-SCOPE (dead/malformed test)** | GNU, jan-2014 GNU and Pro all reject the record; no dof-name variant documented |
| `ho_mech2` | p-refinement instability test; current GNU gives `vely=-1238`, target `-641`; the Pro's own copy (`test/other/tensile_instability.dat`) uses `-updated_linear` + a different timestep control and a **different target `-1975`** | **OUT-OF-SCOPE (stale target/input)** | jan-2014 aborts (`number_of_intagration_points at least 9`); direct solver also `-1238`; sfnet target `-641` matches no binary |
| `taylor3` | axisymmetric Taylor impact; current GNU stalls at the minimal timestep (`kap=0`); jan-2014 runs but gives `2.56391` vs target `2.58442` (fails); the Pro's own copy uses `inertia_apply`, `mesh -follow_material`, `group_materi_memory -updated`, a 3-value timestep control and `_stop -no` | **OUT-OF-SCOPE (stale input + model/solver regression)** | Pro rejects the GNU input (`options_stabilization`); no binary reaches the target |
| `ho_mech4` | initia declares `materi_strain_total` but not the `materi_strain_plasti` required by `group_materi_plasti_vonmises` (manual Pro 6.750); the Pro's own copy (`test/other/cylinder_buckle.dat`) **does** declare it and uses a different mesh + target (`node_dof 1010 -velx = 0.0118`) | **OUT-OF-SCOPE (stale input/target)** | GNU 2014 source has the same check; with `materi_strain_plasti` added the GNU runs but gives `kap=0.06274` vs target `0.0747218`; the sfnet target is not the Pro's target |

## Fix — `so.cc` `solve()`

```c
- if ( diff*nprinc > band ) band = diff*nprinc;
+ if ( (1+diff)*nprinc > band ) band = (1+diff)*nprinc;
```

Two nodes separated by `diff` node numbers span `diff*nprinc + (nprinc-1)`
consecutive local indices. The old bound was short by `nprinc-1`, so the
first-to-last dof coupling of adjacent nodes (the pressure<->velocity
consolidation term, `nprinc=2`, dof distance 3 > band 2) wrote outside
`mat`/`ldmat` and corrupted the `dgbsv` retry. The Cuthill-McKee ordering
branch already used `(1+diff)*nprinc`. See
`manual-developer/solver_bicg_band_width.md`.

## Corpus input alignment (untracked, for reproducibility)

The corpus directory is gitignored. Two `.dat` edits are needed to run the
A/B (originals preserved in the evidence dir):

1. `incnav3.dat` — one line after `end_initia`:
   `inertia_apply -yes` (documented by the solver-timeout work unit;
   without it the model has no mass term and the solver aborts).
2. `examp2.dat` — two lines after `end_initia`:
   ```
   inertia_apply -yes
   groundflow_consolidation_apply -yes
   ```
   The coupled consolidation transient only exists with the coupling
   switch on (Professional default `-no`, manual 6.556); the Pro's own
   reference run also adds both. Without them the raw input gives
   `pres ~ 0` (no coupling), a different failure.

Baseline A/B uses the same aligned corpus for both binaries, so the only
difference is the `so.cc` fix.

## Reproducibility

- Baseline binary: clean build from `0d0ad0c`; A/B runner
  `/tmp/opencode/final5-evidence-20260914-160940/corpus-ab.sh` over the
  aligned corpus -> `186/11/0`.
- Fixed binary: same worktree rebuilt after the `so.cc` fix -> `187/10/0`.
- Classification files: `corpus-ab-baseline.txt`, `corpus-ab-fixed.txt`.
- `examp2` per-run evidence: `examp2-bandfix/`, `examp2-lapack/`,
  `examp2-superlu/`, `examp2-debug/` in the evidence dir.
