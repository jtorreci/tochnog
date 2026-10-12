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
- `general.cc` `general()` — sets the global (since F4 T6f C72,
  2026-10-13, `c1a1ae9`; before, a convected stress dof also cleared it):

  ```c
  stress_recovery_absolute = !materi_strain_total;
  ```

  and SKIPS the incremental lumped block for `-MATERI_STRESS` dofs when the
  flag is set (the whole recovery is assembled in `materi.cc`). With the
  flag set it also skips the convective and artificial-diffusion terms
  of the `-MATERI_STRESS` dofs (`continue` at the top of the `jdim`
  loop): `materi.cc` writes the absolute numerator into the same
  `element_rhside` slot, so a transport term there would pollute it.
- `materi.cc` `materi()` — with the flag set, the stress block assembles the
  ABSOLUTE numerator and the denominator:

  ```c
  tmp = volume * rec_w * new_sig[idim*MDIM+jdim];
  element_rhside[indx] += tmp;
  element_lhside[indx] += volume * rec_w;      // rec_w = F/volume
  ```

  When the flag is clear (independent strain unknown) the historic
  increment is kept and the denominator/inertia come from `general()`.
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
recovered field. One case keeps the historic incremental lumped update
(until 2026-10-13 a second one did, see below):

- an INDEPENDENT strain unknown exists (`materi_strain_total` and the
  plasti/elasti strain variants, e.g. `examp9`/`examp10`/`examp23`): the
  stress dof evolves with that unknown. The Professional cannot validate
  those models (unsupported `bounda_unknown` / `options_mesh` /
  `contact_geometry`), and the absolute form changes the adaptive-remesh
  path of `examp9`.

The CONVECTED stress dof (`-fixed_in_space` + `convection_apply -yes`,
`general()` `lagrange` 0) was the second case until F4 T6f C72
(2026-10-13). It was wrong twice: (1) the Professional 02-08-2026 does
not transport the nodal stress, its node stress is the recovery of the
integration-point stress, which it does not convect either (a stress
pulse stays bit-frozen; Engram 2775); (2) the incremental update was
inconsistent: its lumped mass used `sri_stress_recovery_weight()`
(F/volume, per-element arithmetic mean) while the convective term used
`h`, so the nodal stress moved at about v*volume_ip (measured 0.055 and
0.011 for volume_ip 0.0625 and 0.0125 at v = 1). Making the convective
term use the same weight was the alternative; it was rejected because a
correctly transported nodal stress would still differ from the
Professional and from the integration-point state that the GNU itself
keeps frozen. Evidence: `test-2026/convected_stress_recovery1` (40
quad4, uniform v = 1, elastic, sigzz pulse): post_point 0.5 / 1.0 / 0.5
at t = 2 = Professional (bit-frozen), before 0.348732 / 0.904484 /
0.6481; a parabolic-shear fluid variant also gives the Professional
0.5 / 1.0 / 0.5 (before 0.4827 / 0.9937 / 0.5173). A/B on every input
with `convection_apply -yes` (test-2014 `cd_method`, `cmat_gate`,
`materi_direct_wall`, test-2026 `convection_apply_yes1`, `stokes_*`,
Professional `validation_1`, `validation_8`, `condif2`, `condif3`):
velocities bit-identical (the nodal stress never feeds back), nodal
stresses move up to 1.3 % of their maximum on `validation_1` /
`stokes_*` only; sfnet `strconv1` (`materi_strain_total`) byte-identical.

Measured (2026-09-15): with the gate, the sfnet 197 and validation-suite 233
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

## Geometry-dependent quad4 / quad9 recovery (F4 T6f K21 WU1, 2026-10-12)

The K21 investigation (Engram 2751) fitted the Professional per-element
operator from its `.dbs` (`element_intpnt_dof` and `element_intpnt_coord`
-> `element_dof`, exact to 1e-10..1e-13) on the sfnet `examp1`/`incnav4`,
the Pro `tube1` and dedicated distorted probes:

| element / rule | Professional per-element nodal value | GNU before |
|---|---|---|
| quad4 2x2 Gauss | linear least squares `1,x,y` in PHYSICAL coordinates, `F_g(n) = [1,x_n,y_n] (A^T A)^-1 [1,x_g,y_g]^T` | iso bilinear `1,xi,eta,xi*eta` (table type 0) |
| quad9 2x2 Gauss (`-minimal`) | same physical linear fit | `h_inol` (volume weighted) |
| quad4 1 point (`-minimal`) | the point value | `h_inol` (volume weighted) |
| axisymmetric quad4 2x2 with a node EXACTLY at r=0 | `sum_g h_n(g) sig_g / sum_g h_n(g)` | table type 0 |

Assembly is the arithmetic mean over the adjacent elements in every case
(unchanged `F/volume` trick). The physical and iso linear fits are equal
on parallelograms (the fit is affine invariant), and the bilinear
interpolation also agrees when the integration-point field has no
`xi*eta` mode, which is why rectangular plane meshes never discriminated
them (C64, C42). The axis rule is triggered only by r == 0.: elements
starting at r = 1e-6 or 1e-3 use the linear fit. It was confirmed for a
viscous fluid (`tube1`) and an elastic axisymmetric probe (displacement
prescribed at every node; the two on-axis elements, one of them
distorted, exact to 1.5e-13 with the `h` rule, 0.31 of the scale off
with the linear fit).

Implementation (`miscel.cc`):

- `recovery_fit_kind()` selects the rule (2D only; quad4 1 point; quad4
  and quad9 with 4 points when `stress_recovery_gauss`).
- `recovery_fit_build()` evaluates the Lagrange shape functions at the
  2x2 Gauss points, maps the points to physical coordinates, and builds
  `F[inol][ipoint]`: the `h` rule on the axis, otherwise the 3x3 normal
  equations on centred and scaled coordinates (Gauss-Jordan, pivot
  < 1e-12 -> not applied, `h_inol` fallback).
- `recovery_fit_factor()` reads the start-of-step node coordinates
  (`NODE`, `VERSION_NORMAL`) and `GROUP_AXISYMMETRIC`, and keeps a
  `thread_local` cache of the last element (type, group, node numbers,
  coordinates), so the many calls per element (stress, strain tensors, IP
  states, `general.cc` denominators) build the operator once. Calls with
  the same node numbers skip the `db()` reads; the coordinates are re-read
  and compared on the `(inol 0, ipoint 0)` call, the first call of every
  element visit, so a moved or refined mesh is always caught (a first
  version that re-read the coordinates on every call cost +12 % wall on
  Pro `ground11_nonsaturated`). The element loop can run in threads
  (`options_processors`), hence `thread_local`.
- `sri_stress_recovery_weight()` and `strain_recovery_weight()` take two
  new arguments, `element_group` and `nodes[]`; every caller in
  `materi.cc` (stress, strain tensors, kappa, C68 IP states) and
  `general.cc` (strain/IP-state lumped denominators, convected stress)
  passes them, so numerator and denominator stay consistent.
- No gate (K21 WU1b, 2026-10-12): the fit applies to every stress dof,
  CONVECTED ones included. WU1 (`2305659`) kept the old operator when the
  stress dofs are convected (`stress_recovery_lagrange`, removed): the
  original test-2014 `cd_method` moved from sigxy 0.744 to 3.64 (target
  1 +- 0.3) with the fit. Both values were meaningless: `-ra 1 4 -ra` is a
  node LIST, so the element rotated rigidly (physical stress 0, the
  Professional gives 0 from t = 0.6), and neither value converged in dt
  or mesh. `cd_method`, `cd_geom` and `materi_direct_mc` were rewritten
  as well-posed simple shears (sigxy = c exactly with both operators).
  The convected (Eulerian, `-fixed_in_space` + convection) stress path:
  since F4 T6f C72 (2026-10-13) the nodal stress is recovered absolutely
  there too (see "Gate rationale"); neither code convects the
  integration-point state (BACKLOG C72 / C73).

Not changed (probes in K21 WU1, Professional = physical linear fit too,
GNU = iso linear table; BACKLOG C70): quad9 with 9 points on curved /
distorted elements (iso fit off by ~1.0 of a scale 8 field), hex8 2x2x2 on
distorted bricks (off by 1.1-1.5 of a scale 10 field). Unmeasured: hex27,
tet10 geometry dependence, the axis rule for quad9 / 1-point quad4 (the
point value is used), the Lobatto rules (C67).

Tests (`validation-suite/test-2026`, RED on 7eed0db, GREEN after; every
stress component of the whole node field = Professional to <= 1e-7 of
its scale): `recovery_examp1_quad4_point_mean1` (node 19 sigxx -0.2350738,
before -0.256815), `recovery_incnav4_quad9_minimal1` (nodes 1/3 sigxy
0.05, before 0.02/0.08), `recovery_tube1_axis1` (node 3 sigxx -159.5805,
before -66.40; nodes 2, 33, 40), `recovery_quad4_distorted_plane1`
(node 21 sigxx 2.0014257, before 1.1829).

Corpus A/B (omen, `f4t6g_rec_runs`, 7eed0db vs `2305659`): no rc change
except the 4 new tests (1 -> 0); 9 tests move closer to the Professional
(`tube1`, `pressed_rubber_cylinder`, `axisym5`, `examp1`, `ho_mech1`,
`incnav4`, `incnav5`, `axisym_quad4_hourglass1`,
`materi_memory_updated_linear_current2`), none farther; wall +0.9 %.

## Files

- `dof.cc` — `parallel_new_dof_diagonal()`
- `general.cc` — `general()`
- `materi.cc` — `materi()`, the stress-dof block
- `elem.cc` — `elem()`, the `ELEMENT_DELETE_FACTOR` block
- `polynom.cc` — `stress_recovery_absolute`
- `tochnog.h` / `tochnog-mod.h` — the extern declaration
- `miscel.cc` — `sri_stress_recovery_weight()` (2026-09-13; simplex and
  1-point bar2 arithmetic mean F4 T6f, 2026-10-11), `recovery_fit_kind()`,
  `recovery_fit_build()`, `recovery_fit_factor()` (K21, 2026-10-12)

## Differences with the Professional version

- Convected stress path (`-fixed_in_space` + `convection_apply -yes`):
  nodal output = Professional since `c1a1ae9` (absolute recovery of a
  non-convected integration-point state in both codes). Professional-side
  residuals on `validation_1` (drift, class B C72-V1; Stokes gap) and the
  original rigid-rotation `cd_method` stay open (BACKLOG C72).
- 3x3 quad9 / hex8 on curved or distorted elements: BACKLOG C70.
- Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
