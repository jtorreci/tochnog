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
