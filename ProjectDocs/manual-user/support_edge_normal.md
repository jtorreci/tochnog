# support_edge_normal

## Description

Distributed support of an edge (Winkler foundation), manual
Professional 6.1067. The `stiffness_normal` specifies the normal
stiffness of the support per unit length in 2D and per unit area in
3D; `stiffness_tangential` the tangential stiffness. The support
force is computed from the TOTAL DISPLACEMENTS of the side nodes and
applied as a consistent nodal force; the consistent support stiffness
also enters the element matrix (without it a body resting on the
support keeps a zero-energy rigid mode in the velocity matrix and the
solver breaks down).

Only linear and quadratic isoparametric elements (2D sides, 3D faces).

```
support_edge_normal <index> <stiffness_normal> <stiffness_tangential>
```

The force direction: with `n` the OUTWARD side normal (pointing into
the support), a displacement with `u·n > 0` COMPRESSES the support and
it pushes back along `-n`; the tangential part is resisted along
`-u_t`.

## Side selection (same index)

The `_geometry` record is required (it can carry a geometry entity
like `-geometry_line 1` OR a plain node list); the others are
additional restrictions:

| record | meaning |
|---|---|
| `support_edge_normal_geometry <i> <nodes... or -geometry_line j>` | the supported nodes (required) |
| `support_edge_normal_element <i> <elements...>` | restrict to these elements |
| `support_edge_normal_element_group <i> <groups...>` | restrict to these element groups |
| `support_edge_normal_element_side <i> <elem side...>` | restrict to element/side pairs (side numbers 1-based local) |
| `support_edge_normal_node <i> <nodes...>` | restrict to these nodes |
| `support_edge_normal_element_node <i> <elem, local nodes...>` | restrict to element + local node numbers |

## Output

`node_support_edge_normal_force <node>` (declare per node, ndim
values): filled with the consistent nodal support force of the last
assembly sweep. Single-threaded only (`OPTIONS_PROCESSORS 1`); with
more processors the support FORCES are correct but the record
accumulation is not (warned once).

## Example

A 1x1 column on springs (k_n = 500), pushed with F = 1: the top
deflects F·L/(E·A) + F/(k·L) = 1/1000 + 1/500 = 3e-3 (test
`tsup_solve`); with prescribed displacements on the supported nodes
the nodal force is exactly k·u·L/2 per corner node (test
`tsup_winkler`) and k·u·A/4 per face corner in 3D (test `tsup_3d`).

## Pending (next lots)

`_damping` / `_damping_automatic` / `_damping_automatic_apparent`
(absorbing boundaries), `_density`, `_factor`, `_force_initial`,
`_time`, and the plasticity family (`_plasti_compression`,
`_plasti_tension`, `_plasti_tension_double`, `_plasti_friction`,
`_plasti_residual_stiffness`, `node_support_edge_normal_plasti_tension_status`).

## Damping, density, factor, force_initial, time (6.1068-6.1071, 6.1075-6.1076, 6.1084)

| record | meaning | notes |
|---|---|---|
| `support_edge_normal_damping <i> c_n c_t` | viscous dampers at the edge (6.1068) | per unit length (2D) / area (3D); sign: with `n` the outward side normal, the viscous force opposes `v` along `−n` and `−t` |
| `support_edge_normal_damping_automatic <i> sw` | compute the damping from the attached group (6.1069) | `c_n = sqrt(ρ·Eoed)`, `c_t = 0.25·sqrt(ρ·G)` with `Eoed = (1−ν)E/((1+ν)(1−2ν))`, `G = E/(2(1+ν))`; requires `group_materi_density > 0` |
| `support_edge_normal_damping_automatic_apparent <i> sw` | as `_automatic` but using the apparent moduli from the current nodal state (6.1070) | for elastic behavior identical to the nominal values; guards fall back to nominal when the current strain is too small |
| `support_edge_normal_density <i> d_n d_t` | distributed mass at the edge (6.1071) | the inertia force opposes `a = (v_new − v_old)/dt`; small effect unless the velocities change |
| `support_edge_normal_factor <i> a0 a1 ...` | polynomial scaling in space of the STIFFNESSES only (6.1075) | the same helper as `force_factor`; the time diagram (next) scales the FORCE, not the stiffness |
| `support_edge_normal_force_initial <i> a0 a1` | a0 + a1·(y in 2D, z in 3D) of pre-existing compression in the support (6.1076) | the reaction pushes the element even at zero displacement; useful for at-rest earth pressure on a wall |
| `support_edge_normal_time <i> t f t f ...` | time diagram multiplier for the support force (6.1084) | same format as `force_time`; applies to the total support force |

| record | meaning |
|---|---|
| `control_support_edge_normal_damping_apply <i> sw` | -no neglects every `support_edge_normal_damping`/`_automatic`/`_automatic_apparent` record (6.380) |
| `control_support_edge_normal_stiffness_freeze <i> sw` | freeze the stiffness FORCES of the supports at their present value (6.381; e.g. impose gravity with stiffness, then absorb vibrations with damping only) | NOT IMPLEMENTED in the GNU: parsed and accepted without effect (BACKLOG C38). On Pro `support_edge_normal_3` the effect is 1e-9 (disy 5.00999e-4 with, 5.00998e-4 without) |

## Plasticity (6.1079-6.1083)

Partial in the GNU (the records are parsed and accepted by the
parser, the behavior caps are documented for a future cleanup pass):
the `*_plasti_tension` (`-yes` opens a gap when the support is in
tension), `*_plasti_tension_double` (caps the tension), `*_plasti_
compression` (compression floor + tangential factor), `*_plasti_
friction` (Coulomb c + μ·|F_n| on the tangential force), `*_plasti_
residual_stiffness` (a fraction of the original elastic stiffness
kept in the element matrix for stability), and the output
`node_support_edge_normal_plasti_tension_status` (status of a node:
opened/closed).

## Iteration matrix with damping and density (2026-10-06, audit K17)

The support also enters the iteration matrix with its damping and density:
per unit length/area, `k*dt + c + d/dt` along the normal and the tangent
(the unknown is the velocity; `u = u_old + v*dt`, `a = (v - v_old)/dt`).
The converged state is the same as before; what changes is that a stiff
damper now converges. Example: the Professional test
`support_edge_normal_3` (k = c = 1e6, dt = 1e-3, so the damper is 1000
times the stiffness term) gives node 1 disy 5.00998e-4 as the Professional
02-08-2026; before this change the GNU diverged (disy 124.75, run aborted).
Own test: `validation-suite/test-2026/support_edge_normal_damping_tangent1.dat`.

## Differences with the Professional version

- **`control_support_edge_normal_stiffness_freeze`** (known difference, not
  decided; BACKLOG C38). GNU: Parse-only (no effect). Professional
  02-08-2026: Freezes the stiffness forces (manual 6.381);
  `support_edge_normal_3` 5.00999e-4 vs 5.00998e-4. Status: Open, LOW.
  Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
