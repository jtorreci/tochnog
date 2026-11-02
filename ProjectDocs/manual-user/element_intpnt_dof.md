# element_intpnt_dof

## Description

Output record (manual Professional 6.437): "Unknowns as saved per element in
the element integration points". One record per element:

```
element_intpnt_dof <element> u_1(p1) ... u_n(p1) u_1(p2) ... u_n(p2) ...
```

`n` is the number of unknowns (`dof_label` order) and the points follow the
integration rule of the element (`element_intpnt_npoint`, see
[element_intpnt_records](element_intpnt_records.md)): `npoint*n` values,
point by point. It is the integration-point state the GNU remembers between
time steps when `global_element_dof_apply` (`options_element_dof`) is `-yes`
(the default): strains, stresses, kappa, damage, history variables, plus the
interpolated principal unknowns.

It is filled after every converged time step; it is not meant to be given
in the input file.

## Why it is useful

It shows the stresses and strains where the material law is evaluated,
without the nodal averaging of `node_dof`.

## Example

`element_intpnt_dof1.dat` (validation-suite/test-2026), one unit quad4 with
the default 2x2 Gauss rule and 10 unknowns: 40 values, equal to the
Professional 02-08-2026 except the displacement slots (see below).

```
element_intpnt_dof  1  0.000446581987384 -0.000186075828077 0.000446581987384 -0.000186075828077 2.33676533861 0.474126300588 0 0.0338661643282 0 0.711189450882 ...
```

## Differences with the Professional version

- The displacement slots (`disx`, `disy`, `disz`) hold the interpolated
  displacement in the GNU; the Professional writes 0 there (DIVERGENCES
  ELEMENT-DOF-DISPLACEMENT, class B).
- Springs, contact springs, trusses and beams get no record (as the
  Professional). Interface and empty elements keep the zero record of the
  initialisation with the maximal length (`npointmax` points, default 27).
- Format change of the GNU .dbs (F5 T2 V5, 2026-11-02): until then the GNU
  wrote this storage under the name `element_dof`, padded to 27 points.
  A .dbs of an older build fed back as input therefore gives the old
  integration-point values to the new (output only) `element_dof` record,
  and the new `element_intpnt_dof` starts from zero.
