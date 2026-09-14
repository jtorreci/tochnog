# DIAG: physics/material-value and state-dof corpus subgroup

Date: 2026-09-14
Baseline commit: `d2aec60` (branch `dev/masin-kernel-surface`)
Subgroup: 8 sfnet corpus tests (material/mesh-value and state-dof failures)
Corpus baseline at `d2aec60`: **173 PASS / 24 RUNFAIL / 0 PARSE** (197 tests).
Evidence: `/tmp/opencode/matstate-evidence-20260914-124615/` (GNU d2aec60 runs,
Professional 02-08-2026 runs, and the pre-alignment binary
`/tmp/opencode/oldrepo/build/tochnog`, dated 2026-09-02).

## Result

All 8 tests are classified **OUT-OF-SCOPE**. Six of them are regressions
introduced by two deliberate, documented Professional-alignment commits; the
other two are long-standing contact/adaptive-refinement failures that also
fail on the pre-alignment binary.

| test | target | d2aec60 | root cause | class |
|------|--------|---------|------------|-------|
| examp11 | fscal 1e-4 | -1.36e-6 | `options_inertia` default -no (0560dcf) | OUT-OF-SCOPE |
| wave1 | scal -1 | 0 | `options_inertia` default -no (0560dcf) | OUT-OF-SCOPE |
| dens1 | dens 0.999002 | 1 | `options_inertia` default -no (0560dcf) | OUT-OF-SCOPE |
| truss6 | velx 0.1 | 2e8 | `options_inertia` default -no (0560dcf) | OUT-OF-SCOPE |
| viscos1 | temp 2 | 2e8 | `options_inertia` default -no (0560dcf) | OUT-OF-SCOPE |
| matrix8 | eptxx -0.01 | 0 | atmospheric cap + no-level forcing (741a490) | OUT-OF-SCOPE |
| examp6 | node 259 x 0.51 | 0.609375 | adaptive remesh + plasticity + contact | OUT-OF-SCOPE |
| examp17 | sigyy -20506 | -21103.3 | contact + p-refinement | OUT-OF-SCOPE |

## 1. Five tests: `inertia_apply` default flipped to -no (commit `0560dcf`)

### Evidence

The five tests are legacy GNU models that rely on the GNU's historic
`inertia_apply` default `-yes`. Commit `0560dcf`
("align inertia_apply and consolidation defaults with the Professional")
changed `general.cc` (and `truss.cc`) to `options_inertia=-NO`, matching the
Professional manual 6.785 ("Default, if inertia_apply is not specified, then
each of switch_0, switch_1 etc. is set to -no").

Adding a single `options_inertia -yes` record makes all five hit their targets
exactly (GNU d2aec60, evidence dir `inertia-<test>/`):

| test | target | with `options_inertia -yes` | old binary (2026-09-02) |
|------|--------|------------------------------|--------------------------|
| examp11 | 1e-4 +-1e-5 | 1.01444e-4 | PASS |
| wave1 | -1 +-1e-2 | -0.990179 | PASS |
| dens1 | 0.999002 +-1e-5 | 0.999001998 | PASS |
| truss6 | 0.1 +-1e-8 | 0.1 | PASS |
| viscos1 | 2 +-1e-4 | 2 | PASS |

The pre-alignment binary `/tmp/opencode/oldrepo/build/tochnog` (2026-09-02)
passes all five without modification, confirming `0560dcf` regressed them.

### Why it is not fixable

Reverting the default to `-yes` is not viable: it regresses the 13 static
corpus tests that `0560dcf` fixed (all in the 197 corpus):
`condif6`, `distri1`, `incnav2`, `matrix1`, `matrix3`, `merge1`, `nodfor1`,
`refine2`, `skip1`, `temp1`, `temp3`, `volume1`, `volume2`
(spot-checked: `temp1`, `condif6`, `matrix1`, `incnav2`, `volume1` all PASS
with the d2aec60 default and FAIL with `options_inertia -yes`).

The Professional behaves the same way: its own test suite explicitly adds
`inertia_apply -yes` to `wave1.dat`, `truss6.dat` and `viscos1.dat`
(`test/other/`), and the GNU corpus files are the older GNU variants without
that record. Running the Professional 02-08-2026 on its own `wave1.dat` with
the `inertia_apply` line removed gives `scal = 0` — identical to the GNU
default. The corpus targets are therefore stale with respect to the
Professional-aligned default.

There is no source-grounded way to satisfy both the stale legacy targets and
the Professional default: the only lever is the global `options_inertia`
default, and the two groups of tests want opposite values.

## 2. matrix8: atmospheric cap applied to the raw solved pressure (commit `741a490`)

### Evidence

`matrix8` (solid + groundflow) prescribes `pres = +2e-2` at node 1 and `0` at
node 2 and expects `eptxx = -0.01`. With d2aec60 the GNU returns `eptxx = 0`.

Commit `741a490` ("pressure post-item dependency monitors and no-level
mechanics pressure") changed `materi.cc` so that, when no phreatic-level /
static-pressure-height record exists, the raw solved pressure is routed
through `groundflow_phreatic_coord()` by forcing `level_found = 1`
(`materi.cc` lines 588-597). `groundflow_phreatic_coord()` clamps the pressure
at the atmospheric value (`groundfl.cc` lines 499-500:
`if (total_pressure>=pressure_atmospheric) total_pressure=pressure_atmospheric`).
With the default atmospheric pressure 0, a positive pore pressure (suction)
is clamped to 0, so the solid sees no pressure gradient and `eptxx` stays 0.

Before `741a490` the raw `pres` dof was used whenever no phreatic level was
found, which is why the pre-alignment binary passes `matrix8`.

### Why it is not fixable

The clamp matches the Professional: running the Professional 02-08-2026 on its
own `matrix8.dat` with the sign flipped to `pres = +2e-2` gives `node_dof 2 =
all zeros` — exactly the GNU result. The Professional's own `matrix8.dat` uses
`pres = -2e-2` and expects `+0.01`; the GNU reproduces that
(`+0.01` with the Professional input, rc=0). The corpus `matrix8.dat` is the
older GNU variant with the opposite pressure sign and a stale target. Reverting
the `level_found` forcing would make the GNU diverge from the measured
Professional behavior and would regress the non-corpus `ground11_nonsaturated`
alignment that `741a490` introduced.

## 3. examp6 and examp17: contact + adaptive mesh refinement

Both fail on the pre-alignment binary (2026-09-02) as well, so they are not
regressions.

- **examp6**: von Mises plasticity with hardening, contact penalty
  (`contact_penalty_velocity 1e16`), two global h-refinements plus
  `control_mesh_split`, a local refinement, a restart and three repeats. The
  target is the x coordinate of remesh-created node 259 (`0.51 +-0.02`). The
  GNU produces `0.609375` with the default and `0.53125` with
  `options_inertia -yes` (the pre-alignment binary also produces `0.53125`):
  the target is missed by 0.00125 even with the legacy default. The value is
  determined by the adaptive remeshing/split position, a large subsystem.
- **examp17**: 2D elasticity with a contact line (`contact_penalty_velocity
  1e8`), three h-refinements and one p-refinement. Target `sigyy = -20506
  +-100`; GNU `-21103.3` (both with and without inertia); the pre-alignment
  binary also fails. The residual is a contact / p-refinement convergence
  effect, a large subsystem.

## Validation

- Full 197 corpus with the d2aec60 worktree binary: 173 PASS / 24 RUNFAIL / 0
  PARSE (unchanged, no code change was made).
- No source files were modified.

## Conclusion

No minimal, source-grounded fix exists for this subgroup that does not
regress other corpus tests. The five inertia tests and matrix8 are the cost of
the deliberate Professional-alignment commits `0560dcf` and `741a490`; their
corpus inputs/targets predate that alignment. examp6/examp17 are pre-existing
contact/adaptive-refinement failures. Fixing any of them requires either
updating the stale corpus inputs to the Professional test-suite variants
(add `inertia_apply -yes`; flip the matrix8 pressure sign) or a dedicated
contact / p-refinement work unit.
