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
`general.cc` uses the same weight. The quad4/quad9/hex8/hex27 Lobatto
rules are aligned since C79 (2026-10-13), the bar2 2-point, bar3, bar4
and quad16 rules since C67 (2026-10-16, next section), the tet4
`-maximal`, prism15 and prism6 `-minimal` rules since C67/C76
(2026-10-16, section "The last rules").

## Bars and quad16 (F4 T6f C67, 2026-10-16)

Measurement: 65 Professional 02-08-2026 probes (`gen67.py`: bar2 x4,
bar3 x3, bar4 x3 elements in 1D, quad16 2x2 elements; unequal element
sizes, non-separable distortion, interior bar nodes off the GLL
positions in physical space, nonlinear velocity prescribed at every node
of a viscous materi; methods default, `-gauss`, `-lobatto`, `-reduced`
factor 2 and 4, each with default/`-minimal`/`-maximal`). The
Professional `element_dof` was compared with candidate per-element
operators and its `node_dof` with the mean / length-weighted mean of
`element_dof`:

| rule | per-element value (residual) | node |
|---|---|---|
| bar2 1 point | point value (0) | mean (1e-13; length-weighted 3e-2) |
| bar2 2 points (Gauss, Lobatto, `-reduced -maximal`) | linear fit = Lagrange extrapolation = point value (1e-15) | mean |
| bar3/bar4/quad16 Gauss (default, `-minimal`) | physical linear fit `1,x(,y)` (5e-13; natural fit 1e-2, quadratic fit 5e-2) | mean |
| bar3/bar4/quad16 Lobatto `-minimal` | physical linear fit (4e-13) | mean |
| bar3/bar4/quad16 `-reduced` (Lobatto points, also with the points at the nodes) | physical linear fit (4e-13; point value 5e-2) | mean |
| bar2/bar3/bar4/quad16 Lobatto `-maximal` | point value at the node (0) | mean |
| quad16 axisymmetric, node at r=0 | 4x4 Gauss (C71 switch), `sum h sig / sum h` (4e-13) | mean |

The bar4/quad16 interior nodes sit at the Lobatto positions
`iso +-1/sqrt(5)` in both codes (`interpolation_polynomial()` builds the
npol>=3 Lagrange basis on the Lobatto points; the Professional
`element_intpnt_h` of the Lobatto rule is the identity), so the
Lobatto `-maximal` points are the nodes. `hex64` is rejected by both
codes (Professional "I do not know : -hex64"; GNU element record length
28 < 65). Truss/beam elements have no recovered nodal stress
(`element_truss_force`, `element_beam_force_moment`).

GNU before: `h_inol` for every rule above except the 1-point bar2, i.e.
a length/area-weighted mean of shape-function weighted values; the bar3
Lobatto `-minimal` mid nodes got no contribution (0). On the probes the
node_dof differed 1e-2..8e-1 of the field scale; after the fix every
node of the 56 supported probes (52 + 4 axisymmetric) equals the
Professional within 4e-12.

Implementation (`miscel.cc`): `recovery_tensor_npol()` (new) gives the
nodes per direction of the tensor elements per dimension (bars in 1D,
quad4/9/16 in 2D, hex8/27 in 3D); `recovery_fit_kind()` returns kind 3
(point value) for the Lobatto rule with `npoint==nnol` and kind 1 (fit)
for `npoint` = npol^ndim or (npol-1)^ndim (> 1) when
`stress_recovery_gauss`; `recovery_fit_shape()` is now dimension-generic
and places the nodes at the Lobatto positions of npol points (the old
`-1 + 2k/(npol-1)` spacing was the same for npol <= 3). `polynom.cc`
`pol()`: every `-reduced` rule that resolves to Lobatto points sets
`stress_recovery_gauss` (before: quad9/hex27 only). No change to the
point rules, so `msf_element_rule()` / `mc_element_rule()` need none.
Tests: `recovery_bar2_gauss2_mean1`, `recovery_bar3_gauss_fit1`,
`recovery_bar3_lobatto_minimal1`, `recovery_bar3_reduced1`,
`recovery_bar4_gauss_fit1`, `recovery_bar4_lobatto1`,
`recovery_quad16_gauss_fit1`, `recovery_quad16_lobatto1`,
`recovery_quad16_lobatto_minimal1`, `recovery_quad16_axisym_lobatto1`
(Professional targets 4e-6; RED rc 1 on 75bc0b5, GREEN).

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

The other Gauss rules followed in C70 (2026-10-13, next section).

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

## The other Gauss rules and the 3D elements (F4 T6f C70, 2026-10-13)

Probes on the Professional 02-08-2026 (prescribed velocity at every node,
viscous fluid, `element_intpnt_dof` + `element_intpnt_coord` ->
`element_dof`; distorted bricks, curved quad9/hex27/tet10 with displaced
mid-side nodes). Residual of the PHYSICAL linear least squares
`1,x,y(,z)` relative to the field scale:

| element / rule | Professional | residual | GNU before |
|---|---|---|---|
| quad9 3x3 | physical linear fit | 1e-12 | natural `1,xi,eta` table |
| hex8 2x2x2 | physical linear fit | 8e-13 | natural `1,xi,eta,zeta` table |
| hex8 1 point (`-minimal`) | the point value | exact | `h_inol` (volume weighted) |
| hex27 3x3x3 | physical linear fit | 4e-13 | natural table |
| hex27 2x2x2 (`-minimal`) | physical linear fit | 6e-13 | `h_inol` |
| tet10 5-point Keast | physical linear fit | 4e-13 | barycentric `1,L1,L2,L3` (= physical only on straight edges) |
| prism6 6 points | physical linear fit (of its own points) | 5e-13 | `h_inol` |
| axisymmetric quad9 3x3 / 2x2, node at r=0 | `sum_g h_n(g) sig_g / sum_g h_n(g)` | 9e-12 | natural table / physical fit |

The node is the arithmetic mean of the element values in every probe
(two elements sharing a face/edge, 1e-13). The natural-coordinate fits
differ by 2-15 % of the field scale on these elements.

Implementation (`miscel.cc`):

- `recovery_fit_kind()` covers quad4/quad9/hex8/hex27 Gauss rules (when
  `stress_recovery_gauss`), the 1-point quad4/hex8 (kind 2, point value),
  tet10 5 points and prism6 6 points (their `pol()` branches have fixed
  rules; prism6 leaves `stress_recovery_gauss` at 0).
- `recovery_fit_shape()` repeats the point layout and the shape functions
  of `pol()` (tensor index xi fastest, `m` Gauss points per direction from
  `npoint`; the tet10 Keast points; the prism6 3 x 2 points).
- `recovery_fit_build()` maps the points to physical coordinates and
  solves the `(ndim+1)` normal equations on centred, scaled coordinates;
  the axis rule applies to every 2D quad (quad4 and quad9).
- The static natural-coordinate table (`sri_recovery_factor()`,
  `sri_recovery_table_build()`, F4 T9) was removed: every case it served
  is now a geometry-dependent fit, and a singular fit falls back to
  `h_inol`.

Prism6 integration points (not changed, outside the recovery): the
Professional `element_intpnt_coord` of a prism6 are the triangle CORNERS
at the two Gauss levels in zeta, the GNU `pol()` uses the interior points
(2/3,1/6,1/6). The integration-point stresses therefore differ on a
non-linear field (distorted-prism probe: nodal sigzz 1.92 vs 1.87, 1.6 %
of the scale), although both recover with the same physical fit. The
test `recovery_prism6_linear1` uses a velocity in the prism6 span
`1,x,y,z,xz,yz` on right prisms (linear exact stress), so both point
sets give the exact nodal values (BACKLOG C75: kept as a deliberate
divergence, project owner 2026-10-13; the corner rule doubles the
bending error of a prism6 cantilever).

Tests (`validation-suite/test-2026`, `39d0138`; RED on 421a3f4 rc 1, GREEN rc 0,
the whole node stress field = Professional to <= 2e-12 of its scale):
`recovery_hex8_distorted1` (node 11 sigxx 1.3071691, before 1.1216),
`recovery_hex8_point_mean1` (node 2 sigxx 2.3571531, before 2.9277),
`recovery_hex27_distorted1` (node 41 sigxx 3.6670886, before 3.4134),
`recovery_tet10_curved1` (node 14 sigxx 2.6187305, before 2.5316),
`recovery_prism6_linear1` (node 1 sigxx 0.55 exact, before 1.0575),
`recovery_quad9_curved1` (node 11 sigxx 2.4332610, before 2.2832),
`recovery_quad9_axis1` (node 10 sigxx 0.9518217, before -0.1123).

Corpus A/B (omen `f4t6g_c70_runs`, 421a3f4 vs `39d0138`; Pro other,
sfnet, test-2026, test-2014): no rc / classification change; nodal
stresses move >= 1e-3 of their scale only in `force15` (0.25 -> 4.8e-11
vs the Professional), `force9` (3.4e-2 -> 1.1e-7), `force11` (3.4e-2 ->
7e-7, C57), `interface_patch`, `update2` and sfnet `silo1`; none moves
away from the Professional. Wall `large1` (hex8 30^3) +1.5 % locally
(the per-element fit replaces the F4 T9 table), corpus totals within
+-2 %.

## The last rules: tet4 `-maximal`, prism15, prism6 `-minimal` (F4 T6f C67 remainder + C76, 2026-10-16)

Enumeration. Every rule `pol()` can produce, against the dispatch of
`sri_stress_recovery_weight()` / `recovery_fit_kind()` at `9a0ca3d`: the
tensor rules (C67/C70/C79/C80), tria3 1 and 3 points, tria6 6 points,
tet4 1 point, bar2 1 point and tet10 5 points were covered; three fell to
the `h_inol` fallback: tet4 `-maximal` (4 nodal points), prism15 (21
points, every `group_integration_points` value) and prism6 `-minimal`
(record ignored, 6 points). tria3/tet4 axisymmetric: `pol()` forces the
1-point rule (covered); tet4 in an axisymmetric group does not arise (3D).

Professional probes (`scratchpad/c76/gen76.py`, `an76.py`; 81 runs:
tria3, tria6, tria3 axisymmetric off the axis and with a column of nodes
on it, tria6 axisymmetric, tet4, tet10, prism6, prism15; default /
`-gauss` / `-lobatto` x default / `-minimal` / `-maximal`; 12-24 distorted
elements of unequal size; a nonlinear velocity prescribed at every node of
a viscous materi). Per element the Professional `element_dof` is compared
with candidate operators applied to its `element_intpnt_dof`, and the
`node_dof` with the mean / volume-weighted mean of `element_dof`:

| rule | Professional points | per-element operator (residual) | assembly |
|---|---|---|---|
| tet4 `-maximal` | 5, the tet10 Keast points (manual 6.623) | physical linear fit (8e-14; constant stress) | mean (4e-13; vol-weighted 1e-1) |
| prism15 default/`-normal`/`-maximal` | 21, Dunavant 7 x 3 Gauss | physical linear fit 1,x,y,z (5e-13) | mean (3e-13; vol-weighted 9e-3) |
| prism15 `-minimal` | 8, 4-point triangle rule (centroid -27/48) x 2 Gauss | physical linear fit (7e-13) | mean |
| prism6 `-minimal` | 2, centroid x zeta 1/2 +- 1/(2 sqrt 3) | `h_n(g)/sum_g h_n(g)` (5e-13; linear fit 6e-2) | mean (3e-13) |

Implementation:

- `polynom.cc` `pol()`: tet4 `-maximal` = Keast 5 points (weights -4/5
  and 9/20, `h` of the linear tet4 at the points, `p` unchanged; point
  order as tet10: xi = L1 -> node 0, eta = L2 -> node 2, zeta = L3 ->
  node 1); prism6 `-minimal` = 2 points, weight 1/4 each (the default is
  now set to `-maximal` before the record is read: the function-local
  default `-MINIMAL` would otherwise select the 2-point rule); prism15:
  the 7-point triangle table corrected (`{1/3,1/3}, {b,a}, {a,b}, {a,a},
  {d,c}, {c,d}, {c,c}`, the Professional order; before `{b,b}` and
  `{d,d}`, a point outside the triangle and a rule not exact for linear
  integrands, monomial error 0.025). `stress_recovery_gauss = 1` for the
  three (dense operator: the materi() stress block visits every pair).
- `miscel.cc` `recovery_fit_kind()`: kind 1 for tet4 nnol 4 npoint 5 and
  prism15 nnol 15 npoint 21; new kind 4 for prism6 nnol 6 npoint 2.
  `recovery_fit_shape()`: tet4 Keast, prism15 21-point serendipity shape
  functions, prism6 2-point. `recovery_fit_build()`: kind 4 reuses the
  axis branch (`F = h_n(g)/sum_g h_n(g)`).
- The replicas `msf_element_rule()` (`calcul_force.cc`) and
  `mc_element_rule()` (`delete.cc`) only cover the tensor elements and are
  unchanged.
- Lumped masses: the tet4 Keast row sum is V/4 per node, as before; the
  prism6 2-point row sum is V/6 per node (positive).

Soundness of the two rules the GNU adopts (`scratchpad/c76/analytic/`):
the prism6 2-point stiffness has 7 zero eigenvalues (1 spurious mode per
element) but no global zero mode in any cantilever mesh tried (10x1x1 to
40x2x4, triangles in xy and xz); constant strain is integrated exactly
(patch test); tip-load cantilever error vs Timoshenko -18.0 / -7.1 / -2.8
% (6-point rule -34.7 / -13.1 / -4.85 %). Its weights are confirmed by a
20x2x2 elastic cantilever whose Professional nodal displacements equal a
2-point FE model to 1.5e-13. The Professional prism15 `-minimal` rule has
3 negative stiffness eigenvalues on a distorted prism15 (the negative
centroid weight) and is NOT adopted (BACKLOG C82).

Verification: every probe node = Professional within 6e-12 for the
rules changed (prism15 before: 24 x the stress scale; tet4 `-maximal`
1e-1; prism6 `-minimal` 3.6e-2); the prism15 elastic block
(`prism15_elastic_shear1`) displacements 3.7 % -> 3e-11. Tests
(`validation-suite/test-2026`): `recovery_tet4_maximal1`,
`recovery_prism15_fit1`, `recovery_prism6_minimal1`,
`prism6_minimal_cantilever1`, `prism15_elastic_shear1` (RED rc 1 at
`9a0ca3d`, GREEN).

Hardcoded data: the Keast, Dunavant and Gauss point tables are repeated
in `pol()` and `recovery_fit_shape()` (pending refactoring: one shared
rule table per element).

## Files

- `dof.cc` — `parallel_new_dof_diagonal()`
- `general.cc` — `general()`
- `materi.cc` — `materi()`, the stress-dof block
- `elem.cc` — `elem()`, the `ELEMENT_DELETE_FACTOR` block
- `polynom.cc` — `stress_recovery_absolute`
- `tochnog.h` / `tochnog-mod.h` — the extern declaration
- `miscel.cc` — `sri_stress_recovery_weight()` (2026-09-13; simplex and
  1-point bar2 arithmetic mean F4 T6f, 2026-10-11), `recovery_fit_kind()`,
  `recovery_fit_shape()`, `recovery_fit_build()`, `recovery_fit_factor()`
  (K21, 2026-10-12; 3D and 3x3 rules C70, 2026-10-13),
  `recovery_tensor_npol()` (bars and quad16, C67, 2026-10-16); kind 4 and
  the tet4/prism15/prism6 shapes (C67 remainder + C76, 2026-10-16)
- `polynom.cc` — `pol()`: tet4 `-maximal` Keast rule, prism6 `-minimal`
  2-point rule, prism15 triangle table (2026-10-16)

## Differences with the Professional version

- Convected stress path (`-fixed_in_space` + `convection_apply -yes`):
  nodal output = Professional since `c1a1ae9` (absolute recovery of a
  non-convected integration-point state in both codes). Professional-side
  residuals on `validation_1` (drift, class B C72-V1; Stokes gap) and the
  original rigid-rotation `cd_method` stay open (BACKLOG C72).
- 3x3 quad9, hex8, hex27, tet10, prism6: = Professional since C70
  (2026-10-13). The prism6 integration points differ (Professional
  triangle corners, GNU interior points; deliberate divergence, owner
  2026-10-13, BACKLOG C75), so prism6 nodal stresses agree only for a
  linear stress field. prism6 `-minimal` = Professional since C76
  (2026-10-16); every rule now has a measured recovery (quad/hex Lobatto
  C79, bars and quad16 C67, tet4 `-maximal` and prism15 C67 remainder).
- Integration rules still different (the recovery of each GNU rule is
  the Professional one): tria6 (GNU 6 nodal points, manual and
  Professional 7 / 4 points; axisymmetric refused by the GNU; BACKLOG
  C66), prism15 `-minimal` (Professional 8 points, indefinite stiffness;
  GNU 21; BACKLOG C82), tria3 `-maximal` (manual 4 points, both 3 nodal;
  axisymmetric GNU 1 point vs Professional 3 nodal / 1 on the axis) and
  tet10 `-maximal` (manual 10 points, both 5) (BACKLOG C83).
- `hex64`: rejected by both codes (the GNU source has the npol=4 branch
  in `pol()`, but the element record is too short); not a divergence.
- Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
