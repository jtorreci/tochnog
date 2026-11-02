# dof_element_dof

## Description

Print-only record (manual Professional 6.409):

```
dof_element_dof dof_per_element_0 dof_per_element_1 ...
```

One value per dof, in `dof_label` order: `-yes` when the dof is by default
calculated per element (non-continuous field, remembered in the element
integration points, `element_intpnt_dof`), `-no` when it is a continuous
(nodal) field. It is written to the database after the input; it is not an
input record.

The `-yes` dofs are the integration-point states of the material update:
`materi_stress`, `materi_strain_elasti`, `materi_strain_plasti`,
`materi_strain_total`, `materi_strain_intergranular`, the history variables
(`materi_history_variables`, cam-clay and hypoplastic history),
`materi_plasti_kappa`, `materi_plasti_kappa_shear`, the cap history,
`materi_damage`, `materi_plasti_rho`, `materi_maxwell_stress` and the
softening variables. Every other dof (velocity, displacement, temperature,
pressure, work, ...) is `-no`.

The record shows the default classification: `global_element_dof_apply -no`
(elements use the averaged nodal values) does not change it.

## Why it is useful

It tells which `node_dof` values are averages of element values and which
ones are solved or integrated in the nodes, which matters when reading the
`element_intpnt_dof` / `element_dof` records.

## Example

`validation-suite/test-2026/dof_element_dof1.dat` (quad4 with velocity,
displacement, total strain, stress, kappa and temperature):

```
dof_label   -velx -vely -disx -disy -eptxx -eptxy -eptxz -eptyy -eptyz -eptzz -sigxx -sigxy -sigxz -sigyy -sigyz -sigzz -kap -temp
dof_element_dof   -no -no -no -no -yes -yes -yes -yes -yes -yes -yes -yes -yes -yes -yes -yes -yes -no
```

Same list as the Professional 02-08-2026. `dof_element_dof2.dat` gives the
same record with `global_element_dof_apply -no`.

## Differences with the Professional version

- The Professional marks every non-principal dof except the displacement
  `-yes` (also `materi_work`, `materi_plasti_f`, `materi_void_fraction`); the
  GNU marks only the dofs it really remembers in the integration points
  (the list above), so `work`, `f` and `void` are `-no` in the GNU
  (DIVERGENCES DOF-ELEMENT-DOF-SET, class B).
