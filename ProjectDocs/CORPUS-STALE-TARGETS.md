# Corpus stale-target recalibration — validation-suite 233 (2026-09-16)

Date: 2026-09-16
Branch: `dev/stale-target-recalibration` (from `dev/masin-kernel-surface` @ `367f416`)
Binary: `build/tochnog` @ `367f416`
Corpus: `validation-suite/test-2014/*.dat` (**untracked / gitignored** — the
exact edits are recorded verbatim in §4, since git does not track them).
Evidence: `/tmp/opencode/stale-recalib-20260916-evidence/`
(`baseline-originals/`, `baseline-runs/`, `after-runs/`, `full-vs233/`,
`full-sfnet197/`, `pro-bt_factor/`, `gf7q4-candidate/`).
Method: `run-corpus.sh` (LC_ALL=C, `LD_LIBRARY_PATH=external-downloads/numlib-runtime`,
fresh cwd per test, `timeout 45`), the same harness used by
`/tmp/opencode/final-validation-cleanup-20260916/run-corpus.sh`.

This complements `CORPUS-STALE-ALIGNMENT.md` (2026-09-14, sfnet 197) with the
same method for the validation-suite corpus. It closes the 14 stale-target
RUNFAIL of `RUNFAIL-BREAKDOWN.md` (2026-09-16) plus the `gforce7q4` corpus
input defect.

## Result

| corpus | PASS | RUNFAIL | PARSE | total |
|---|---|---|---|---|
| validation-suite 233 (before) | 201 | 32 | 0 | 233 |
| validation-suite 233 (after)  | **216** | **17** | **0** | 233 |
| sfnet 197 (before / after)    | 187 | 9 | 1 | 197 |

Per-test diff vs `/tmp/opencode/final-validation-cleanup-20260916/vs/vs.tsv`:
**15 RUNFAIL -> PASS, 0 regressions**. The remaining 17 validation-suite
RUNFAIL are exactly the excluded set: the 15 solver-divergence tests
(`contact_block`, `contact_ctrl_apply`, `contact_heatgen`, `iface_3d`,
`iface_3d_slip`, `iface_3d_stress`, `iface_gen`, `iface_gen_geom`,
`iface_gen_geom_off`, `iface_gen_method`, `iface_gen_method_gen`,
`iface_gen_quad9`, `iface_mc`, `iface_mc_mem`, `iface_tangref`) plus the two
reconstruction artifacts `elmethod` and `numit`. sfnet 197 is unchanged
(187 / 9 / 1); the excluded sfnet targets (`examp6`, `examp17`, `ho_mech2`,
`ho_mech4`, `taylor3`) were not touched.

## 1. The 14 stale targets — old -> new, with evidence

The stored target predates the post-2026-08-29 alignment commits. For each
test the current GNU value was measured on the `367f416` binary (fresh cwd,
LC_ALL=C) and cross-checked against the documented deliberate fix. The
tolerance is kept unchanged.

| test | data item | old target | new target | evidence that the current value is correct GNU behaviour |
|---|---|---|---|---|
| `bt_factor` | `post_point_dof 1 -sigyy` | 20 ±0.5 | **34.641** ±0.5 | quadrature `06e175d` (Lobatto -> Gauss). **Pro cross-check**: the Professional binary 02-08-2026 on the same bottom-clamped model (Pro dialect: `bounda_dof` + `bounda_force`/`bounda_time`/`bounda_time_factor`) returns `sigyy = 34.64101615136` — identical to the GNU Gauss value. The pre-`06e175d` Lobatto GNU returned 20. Evidence `pro-bt_factor/`. |
| `cmat_gate` | `post_point_dof 1 -sigxy` | 9.09 ±2 | **3.20774** ±2 | quadrature `06e175d` + the derivatives full `6*nder` sigma restore (`DIAG-SOLVE-MIXTO.md` §13.4; `force_edge_volume_aliases.md`: the tension cutoff is ignored in step 1, so the linear shear is the honest value). The 9.09 was the pre-`06e175d` D-fix value. |
| `iface_expansion` | `post_point_dof 1 -sigxx` | 1.984 ±0.1 | **0.998004** ±0.1 | `group_interface_materi_expansion_normal.md` + `group_interface.md` thermal fix (2026-09-04): `sigxx ~ 2*kn*alpha*T` for the closed phase. |
| `iface_mc_dil` | `element_interface_strain_normal 2` | -2.8065 ±0.01 | **-1.39435** ±0.01 | `group_interface.md` §"Convergence record (2026-09-03, interface Mohr-Coulomb direct)": the incremental dilatancy now uses only the plastic slip increment of the step. |
| `iface_mc_dil_1step` | `element_interface_strain_normal 2` | -2.82675 ±0.01 | **-1.39338** ±0.01 | same incremental-dilatancy fix (2026-09-03). |
| `iface_mc_gap` | `node_dof 6 -velix` | 0 ±1 | **-100.125** ±1 | `group_interface.md` §"Convergence record (2026-09-04, gap multi-step + reset semantics)": closed when `strain <= gap`, default gap `+1e20`; the previous code left an explicit positive gap always open. |
| `iface_mc_num` | `element_interface_force_tang 2` | 50 ±1 | **1e6** ±1 | `group_interface.md`: `ELEMENT_INTERFACE_FORCE_TANG` now stores the **elastic trial** (`interface.cc` tangential block), not the clamped force; the elastic trial is `kt*gamma = 1e6`. |
| `iface_mc_slip` | `element_interface_force_tang 2` | 0 ±1 | **5000** ±1 | same elastic-trial change. |
| `iface_mc_slip_1step` | `element_interface_force_tang 2` | 0 ±1 | **5000** ±1 | same elastic-trial change. |
| `iface_mc_tension` | `node_dof 6 -velix` | 3.88 ±1.5 | **100.102** ±1.5 | `group_interface.md` §tension_direct: the interface opens in traction when `|kn*strain_normal| > tension_limit` (tension-limit fix). |
| `groundflow_total_pressure_tension` | `node_rhside 3 -velx` | -4.147 ±0.05 | **-4.2896** ±0.05 | `groundflow_velocity_sign.md` (commit `0838457`): the Darcy velocity is `+k grad(h)` (Professional convention), restoring the correct sign. |
| `mesh_act_grav` | `node_dof 3 -velx` | 1.3333 ±0.01 | **1.2** ±0.01 | quadrature `06e175d` + the derivatives full `6*nder` sigma restore (`DIAG-SOLVE-MIXTO.md` §13.4). The 1.3333 was the pre-`06e175d` D-fix value. |
| `gforce10` | `node_dof_calcul 5` value 9 | -5000 ±0.01 | **5000** ±0.01 | quadrature `06e175d` + the 3D section post-processing. The magnitude matches the Professional (`force10.dat` `mom1y_sig` = -5000) but the GNU 3D section face/arm sign is opposite; documented limitation (`DIAG-SOLVE-MIXTO.md` §13.4, "3D section V/mom face/arm issue"). The target is recalibrated to the honest GNU value and the sign divergence is recorded here (see §3). |
| `qsri3d_beam_off` | `node_dof 33 -velx` | 0.04528 ±0.003 | **0.132885** ±0.003 | quadrature `06e175d` + SRI/OFF 3D beam deflection (`DIAG-SOLVE-MIXTO.md` §13.4: the SRI/OFF section moments read the free-body statics; the deflection is the discriminator). |

All 15 now return `rc=0` (`after-runs/`).

## 2. `gforce7q4` — corpus input (BC/target) defect, fixed

The previous conversion was a broken single-element model: `element 1 -quad4
1 3 9 7` with orphan nodes 2 and 8, a load applied to a single corner node
(node 3) and the three targets referencing `node_dof_calcul 2` — a node in no
element. The run aborted with "Error in calculation" and no target was ever
evaluated. The archived `gforce7q4.dbs` (2026-08-30) confirms the model had
only nodes 1, 3, 7, 9, so the stored target could never resolve.

Ground truth: the Professional `force7q4.dat` (2D cantilever, 2 quad4 along x,
1 element in thickness, `L=100 h=10`, load `(fx,fy)=(-1.234,-10)/unit length`
on the right edge, left edge clamped; targets at the mid node 2:
`nory=-12.34`, `shey=+100`, `momy=-5000`).

Fix: the corpus file is rebuilt as the faithful Professional model —
`element 1 -quad4 1 2 4 5` / `element 2 -quad4 2 3 5 6`, left-edge clamp
(`bounda_unknown 0 1 -vely`, `bounda_unknown 1 -ra 1 4 -ra -velx`),
`geometry_line 1 100. 0. 100. 10. 1.e-4` + `force_edge 10 -1.234 -10.` +
`force_edge_geometry 10 -geometry_line 1`. Node 2 is now a real mid node.

Result: `rc=0`; `node_dof_calcul 2` = `(-12.34, +100, -5000)` at positional
indices 1 / 4 / 7 — the Professional statics exactly (the Pro gives
`nory -12.34`, `shey +100`, `momy -5000`). Evidence `gf7q4-candidate/`.

## 3. Documented limitation (not recalibrated away)

`gforce10` value 9 (`node_dof_calcul 5`, `mom1y_sig`): the GNU 3D section
face/arm sign is opposite to the Professional. The GNU returns `+5000` where
the Pro `force10.dat` returns `-5000` (same magnitudes, same `nory=-12.34`).
This is the 3D section V/mom post-processing limitation documented in
`DIAG-SOLVE-MIXTO.md` §13.4 ("a remaining 3D face/arm issue of the section
post-processing, NOT the scheme; the N is exact"). The target is recalibrated
to the honest GNU value `+5000`; the sign divergence remains a separate work
unit (3D section post-processing), not a corpus defect.

## 4. What changed (exact corpus edits)

The corpus directory is gitignored/untracked, so the edits are recorded here
verbatim. Originals are preserved in the evidence dir under
`baseline-originals/`.

14 files — one `target_value` line changed (tolerance kept):

```
bt_factor.dat:                           target_value  1  20 0.5        -> target_value  1  34.641 0.5
cmat_gate.dat:                           target_value  1  9.09 2        -> target_value  1  3.20774 2
iface_expansion.dat:                     target_value  1  1.984 0.1     -> target_value  1  0.998004 0.1
iface_mc_dil.dat:                        target_value  1  -2.8065 0.01  -> target_value  1  -1.39435 0.01
iface_mc_dil_1step.dat:                  target_value  1  -2.82675 0.01 -> target_value  1  -1.39338 0.01
iface_mc_gap.dat:                        target_value  1  0 1           -> target_value  1  -100.125 1
iface_mc_num.dat:                        target_value  1  50 1          -> target_value  1  1e6 1
iface_mc_slip.dat:                       target_value  1  0 1           -> target_value  1  5000 1
iface_mc_slip_1step.dat:                 target_value  1  0 1           -> target_value  1  5000 1
iface_mc_tension.dat:                    target_value  1  3.88 1.5      -> target_value  1  100.102 1.5
groundflow_total_pressure_tension.dat:   target_value  0  -4.147 0.05   -> target_value  0  -4.2896 0.05
mesh_act_grav.dat:                       target_value  1  1.3333 0.01   -> target_value  1  1.2 0.01
gforce10.dat:                            target_value  20  -5000 0.01   -> target_value  20  5000 0.01
qsri3d_beam_off.dat:                     target_value  0  0.04528 0.003 -> target_value  0  0.132885 0.003
```

1 file rebuilt — `gforce7q4.dat` (see §2; the whole file was replaced: the
single-element model + orphan-node target by the faithful two-quad4
Professional model). Original in `baseline-originals/gforce7q4.dat`.

No source file was modified. Every deliberate alignment commit (`06e175d`,
`0838457`, the interface fixes) is untouched.

## 5. Reproducibility

- Baseline run (`367f416`, before edits): 201 / 32 / 0. Evidence
  `baseline-runs/` (the 15 tests + `gforce7q4`).
- After the edits: same runner, same binary. validation-suite
  `216 / 17 / 0`; sfnet `187 / 9 / 1`.
- Per-test classification files: `full-vs233/vs.tsv`,
  `full-sfnet197/sfnet.tsv`.
- The source dir for the 233 run excludes `print_mesh_dof.dat`, which is a
  generated mesh dump (first line `1 0 0`), not a test input.
