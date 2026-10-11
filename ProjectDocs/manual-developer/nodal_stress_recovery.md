# Nodal stress recovery (absolute weighted average)

2026-09-15. Closes the `delete2` `node_dof` target (`node 2 -sigxx = 0.8`) and
aligns the whole `node_dof`/`post_point_dof` stress field of the
`delete`/`mpc` family with the Professional binary (02-08-2026).

## Problem

The GNU nodal stress recovery (`dof.cc parallel_new_dof_diagonal`) accumulated
the element stress INCREMENTS:

```
node_dof_new[stress] += node_rhside / node_lhside
```

with `materi.cc` assembling `node_rhside += volume*weight*(new_sig-old_sig)/dt`
and `general.cc` adding the incremental lumped inertia
`-volume*weight*(new_dof-past_dof)/dt` plus `volume*weight/dt` to `lhside`.

At a node shared by two elements whose increments cancel (+0.2 / -0.2) the
reported stress stays at the step-start value. On `delete2` (two `bar2` in
series, right bar deleted 40 % then reactivated) the shared node 2 reported
`0.6` while the Professional reports `0.8`.

The Professional ASSIGNS the absolute weighted average of the per-element
nodal stresses (its `element_dof`):

```
node_dof[stress] = sum_e (element_dof_e) / N_e
```

## The Professional's weights (measured from its `.dbs`)

`delete2` final state (Pro `delete2.dbs`):

| record | node 1 | node 2 | node 3 |
|---|---|---|---|
| `element_dof 1 -sigxx` | 0.8 | — | — |
| `element_dof 2 -sigxx` | — | — | 0.8 |
| `node_rhside -sigxx` | 0.8 | **1.6** | 0.8 |
| `node_lhside -sigxx` | 1.0 | **2.0** | 1.0 |
| `node_dof -sigxx` | 0.8 | **0.8** | 0.8 |

`node_lhside` for a stress dof equals the NUMBER of adjacent elements (each
element contributes weight 1; the per-element fit has `sum_points F = 1`).
The GNU 2026-09-13 alignment already returns that fit divided by the point
volume (`miscel.cc sri_stress_recovery_weight`), so `volume*weight = F` and
each element contributes `sum_g F = 1`. The only missing piece was that GNU
ACCUMULATED instead of ASSIGNED.

`delete3` (quad9, top element deleted 40 %, no reactivation) confirms the
deleted-element rule: the Pro `node_lhside` for the stress dofs is `1` on the
shared nodes (element 1 only) and `0` on the nodes exclusive to the deleted
element; `node_dof` is `0.6` on the shared nodes and keeps its reset value
`1.0` on the exclusive nodes.

## Implementation

- `polynom.cc` — new global `stress_recovery_absolute` (declared in
  `tochnog.h` / `tochnog-mod.h`).
- `general.cc` `general()` — sets the global from the same condition as the
  local `lagrange` flag:

  ```c
  stress_recovery_absolute = ( !materi_strain_total &&
    ( options_convection==-NO || materi_displacement ||
      ( materi_velocity && options_mesh[0]==-FOLLOW_MATERIAL ) ) );
  ```

  and SKIPS the incremental lumped block for `-MATERI_STRESS` dofs when the
  flag is set (the whole recovery is assembled in `materi.cc`).
- `materi.cc` `materi()` — with the flag set, the stress block assembles the
  ABSOLUTE numerator and the denominator:

  ```c
  tmp = volume * rec_w * new_sig[idim*MDIM+jdim];
  element_rhside[indx] += tmp;
  element_lhside[indx] += volume * rec_w;      // rec_w = F/volume
  ```

  When the flag is clear (convected stress dof) the historic increment is
  kept and the denominator/inertia come from `general()`.
- `dof.cc` `parallel_new_dof_diagonal()` — with the flag set, the stress dofs
  are ASSIGNED (`node_dof_new = node_rhside/node_lhside`) instead of
  accumulated. Every other non-principal dof keeps `+=`.
- `elem.cc` `elem()` — a currently deleted element is EXCLUDED from the
  recovery: after the frozen released load is applied, its stress slots in
  `element_rhside` are zeroed (its `element_lhside` was already zeroed). The
  frozen `element_rhside_delete` is a MOMENTUM right-hand side and must not
  enter the stress numerator. Without this, the reverted 2026-09-15 naive
  attempt polluted the average and regressed `delete3` (0.6 -> 0.12).

The numerator and denominator use the SAME 2026-09-13 least-squares weight,
so the node is the per-element arithmetic mean, exactly as the Professional.

## Gate rationale

The absolute assignment is only valid when the stress dofs are a PURE
recovered field. Two cases keep the historic incremental lumped update:

- the stress dof is CONVECTED (the `general()` `lagrange` flag is 0, e.g.
  `strconv1`): the assignment would discard the transport terms;
- an INDEPENDENT strain unknown exists (`materi_strain_total` and the
  plasti/elasti strain variants, e.g. `examp9`/`examp10`/`examp23`): the
  stress dof evolves with that unknown. The Professional cannot validate
  those models (unsupported `bounda_unknown` / `options_mesh` /
  `contact_geometry`), and the absolute form changes the adaptive-remesh
  path of `examp9`.

Measured: with the gate, the sfnet 197 and validation-suite 233
classifications are byte-identical to the 73fb304 baseline and the
`node_dof`/`post_point_dof` fields of all 191 sfnet tests that dump a `.dbs`
are unchanged (max abs diff 0). Only the `delete`/`mpc` family changes, and
there it matches the Professional exactly.

## Validation (binary built from 73fb304, clean `build_safe.sh`, `LC_ALL=C`)

`node_dof` + `post_point_dof` vs the Professional 02-08-2026 `.dbs`:

| test | GNU | Pro | max abs diff over the full field |
|---|---|---|---|
| delete1 | rc=0 | rc=0 | 0 |
| delete2 | rc=0, node 2 sigxx **0.8** | rc=0 | 0 |
| delete3 | rc=0, sigyy **0.6** | rc=0 | 1.1e-14 |
| delete4 | rc=0 | rc=0 | 0 |
| mpc1..mpc6 | rc=0 | rc=0 | <= 2.1e-10 |
| mpc7 | rc=0 | (no `.dbs`: `PRINT_APPLY -NO`) | — |

Corpus: sfnet 197 = 187 PASS / 9 RUNFAIL / 1 PARSE, 0 per-test class changes
vs the 73fb304 baseline; validation-suite 233 = 201 PASS / 32 RUNFAIL /
0 PARSE, 0 changes; internal suite 16/16 (the file-check FALLO list is
identical to the baseline). `strconv1` (gated) is byte-identical.

## Residual (separate work unit)

The one-iteration reporting lag of `mpc5` at `k=1` is NOT closed by this
change: the stress dofs are still non-principal and updated after the solve
from the PRE-solve constitutive stress, so `mpc5` reports `1.0` at `k=1` and
`0.5` at `k=2` while the Professional reports `0.5` at `k=1`. The test passes
at its own iteration count (2). Closing the lag needs a post-solve
constitutive re-evaluation (candidate (b) of
`mpc5_staggered_stress_lag.md`), not a different recovery weight.

## Simplex elements and the 1-point bar2 (F4 T6f WU1, 2026-10-11)

Until 2026-10-11 `sri_stress_recovery_weight()` returned `h_inol` for
every element outside the Gauss table (quad4/quad9/hex8/hex27/tet10). With
the lumped numerator/denominator multiplying the weight by the point
volume, the tria3/tria6/tet4/bar2 nodes were the AREA (volume, length)
weighted average of the adjacent elements. The Professional 02-08-2026
assembles them by the ARITHMETIC MEAN, like every other element (the
function's own design comment already said so).

Evidence (two-element probes of different size and stiffness, Professional
`element_dof` = per-element nodal value, `node_dof` compared with the mean
and the area-weighted mean of those values):

| element / rule | Professional `node_dof` | GNU before |
|---|---|---|
| tria3 1 point (default) | arithmetic mean, exact (1e-13) | area-weighted, exact |
| tria3 `-maximal` (3-point nodal rule) | arithmetic mean, exact | area-weighted, exact |
| tria6 (6-point nodal rule) | arithmetic mean, exact | area-weighted (element values ~1e-5 off, see BACKLOG C66) |
| tet4 1 point (default) | arithmetic mean, exact | volume-weighted, exact |
| bar2 1 point (default) | arithmetic mean, exact | length-weighted, exact |

Fix (`miscel.cc` `sri_stress_recovery_weight()`, before the Gauss-table
dispatch): for `tria3` with 1 or 3 points, `tria6` with 6 points, `tet4`
and `bar2` with 1 point the weight is `F/volume` with `F = 1` for the
1-point rules and `F = h_inol` (the Kronecker delta) for the nodal rules,
so each element contributes weight 1. The strain tensors follow through
`strain_recovery_weight()` (same function); the convected-stress path of
`general.cc` uses the same weight. Unmeasured rules (tet4 `-maximal`,
bar2 with 2 points, bar3/bar4, the quad Lobatto rules) keep `h_inol`
(BACKLOG C67).

Effect on `validation_14` (DIVERGENCES K15): the pass-1 hole-boundary
nodal stress / plastic strain (up to ~20 % at 8 hole nodes) was exactly
this weighting; with the fix the whole refined mesh (586 nodes, 1118
elements) equals the Professional and every `node_dof` is within 1e-4.

Tests (`validation-suite/test-2026`, targets = Professional, RED on the
previous binary, GREEN after): `tria3_nodal_recovery_mean1`,
`tria3_nodal_recovery_mean2`, `tria6_nodal_recovery_mean1`,
`bar2_nodal_recovery_mean1`, `tet4_nodal_recovery_mean1`; and
`validation-suite/test-2014/validation_14_pass1` (hole nodes 4 and 120).

## Files

- `dof.cc` — `parallel_new_dof_diagonal()`
- `general.cc` — `general()`
- `materi.cc` — `materi()`, the stress-dof block
- `elem.cc` — `elem()`, the `ELEMENT_DELETE_FACTOR` block
- `polynom.cc` — `stress_recovery_absolute`
- `tochnog.h` / `tochnog-mod.h` — the extern declaration
- `miscel.cc` — `sri_stress_recovery_weight()` (2026-09-13; simplex and
  1-point bar2 arithmetic mean F4 T6f, 2026-10-11)
