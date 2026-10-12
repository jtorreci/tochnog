# Nodal stress recovery

The stress reported at a node (`node_dof`, `post_point_dof` and the stress
fields written to `.flavia.res`) is the ABSOLUTE weighted average of the
per-element nodal stresses of the adjacent elements:

```
sigma_node = sum_e (element stress of e at the node) / (number of adjacent elements)
```

This is the Tochnog Professional convention (each element contributes weight
1, the per-element value is a least-squares fit of its integration-point
stresses). Before 2026-09-15 the GNU version accumulated the element stress
INCREMENTS instead, so at a node shared by two elements whose increments
cancel the reported stress kept its previous value.

The same arithmetic mean applies to every element type. Since 2026-10-11
(F4 T6f) this includes the `tria3` (default 1 point and `-maximal`), the
`tria6`, the `tet4` (1 point) and the 1-point `bar2`, which before were
AREA (volume, length) weighted: at a node shared by a small stiff element
and a large soft one the reported stress is now the plain mean of the two
element values, as in the Professional. Example:
`validation-suite/test-2026/tria3_nodal_recovery_mean1.dat` (two tria3 of
area 0.5 and 2 sharing an edge): node 2 `sigxx` 6.4643 = (1.0938 +
11.8348)/2 (before: 9.6866). The strain tensors follow the same rule.
Not yet measured, and therefore unchanged: tet4 `-maximal`, bar2 with 2
points, bar3/bar4.

## quad4 and quad9 per-element values (K21, 2026-10-12)

Since F4 T6f K21 the per-element nodal value of the 2D Gauss rules is the
one of the Professional 02-08-2026:

- `quad4` with 2x2 Gauss points (default) and `quad9` with 2x2 points
  (`group_integration_points -minimal`): a least-squares plane
  `a + b x + c y` through the integration-point values, fitted in the
  PHYSICAL coordinates and evaluated at the node. On rectangles and
  parallelograms this is the same as before; on distorted elements the
  nodal stresses change (up to tens of per cent of the local value).
- `quad4` with 1 point (`-minimal`): the point value, so the node is the
  plain mean of the adjacent element values (before: area weighted).
- axisymmetric `quad4` element with a node exactly on the axis (r = 0):
  the shape-function weighted average of its point values. Elements that
  only come close to the axis (r = 1e-6) use the plane fit.

Example: `validation-suite/test-2026/recovery_tube1_axis1.dat` (laminar
tube flow, axisymmetric quad4): node 3 sigxx -159.58 (before -66.40),
equal to the Professional. Other examples:
`recovery_quad4_distorted_plane1.dat`, `recovery_examp1_quad4_point_mean1.dat`,
`recovery_incnav4_quad9_minimal1.dat`. The strain tensors and the other
recovered integration-point states follow the same rule, also when the
stress is CONVECTED (`convection_apply -yes` on a fixed mesh with
`materi_velocity`; since 2026-10-12, K21 WU1b, no exception).

Difference with the Professional version: on the convected (Eulerian)
stress path the GNU does not convect the integration-point stress state,
and the nodal stress is an incremental transported output that is not
consistent with it; the Professional sweeps the stress out of the fixed
mesh (BACKLOG C72, class B). Use a Lagrangian mesh
(`options_mesh -follow_material`) when the stress history matters.

Not yet aligned (the Professional also fits in physical coordinates
there): `quad9` with 3x3 points and `hex8` on curved or distorted
elements; see BACKLOG C70.

## When it matters

The difference is visible in models where a node is shared by elements with
opposite stress increments, in particular partially deleted / reactivated
elements (`control_mesh_delete_geometry_factor` +
`control_mesh_delete_geometry_stop`). Example: `delete2.dat` (two `bar2` in
series, right bar deleted 40 % then reactivated, `control_reset_dof -sigxx 1`)
now reports the shared node `2` with `sigxx = 0.8`, matching the Professional;
before it reported `0.6`.

A currently deleted element is EXCLUDED from the recovery: it contributes
neither the numerator nor the denominator, so a node whose only adjacent
element is deleted keeps its previous stress and a shared node is recovered
from the active elements alone (`delete3.dat`: `sigyy = 0.6` on the shared
nodes, `1.0` on the nodes exclusive to the deleted element).

## Scope

The absolute average is applied when the stress is a purely recovered field.
Two cases keep the historic incremental update because the stress there is
not a pure post-processing recovery:

- the stress is CONVECTED (e.g. `strconv1`, a `materi_strain_total` transport
  model with `options_mesh -fixed_in_space`);
- an independent strain unknown exists (`materi_strain_total` and the
  plasti/elasti strain variants).

These cases have no Professional reference (their keyword sets are not
supported by the Professional binary) and are unchanged.

## Residual

`mpc5` still reports its stress one iteration behind the velocity (1.0 at
`k=1`, 0.5 at `k=2`; the Professional reports 0.5 at `k=1`). The test passes
at its own iteration count; closing the one-iteration lag needs a post-solve
constitutive re-evaluation (see the developer manual).
