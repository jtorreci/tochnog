# mesh_activate_gravity_time

## Description

`mesh_activate_gravity_time` gradually activates the **gravity** for
elements between `time_start` and `time_end`, from the bottom to the top of
the mesh (typical for dam/dumping construction). The element activation
time is interpolated from the global window and the lowest coordinate of
the element. Between the element start and end activation times the gravity
ramps from 0 to its full value.

The elements to which the record applies are selected by:
- `mesh_activate_gravity_element` (element range),
- `mesh_activate_gravity_element_group` (element groups),
- `mesh_activate_gravity_geometry` (elements inside a geometry).
Without a selection record, all elements apply.

`mesh_activate_gravity_time_initial` sets a time of birth before which the
element is completely inactive. `mesh_activate_gravity_time_strain_settlement
-yes` marks the record as combined with `strain_settlement_parameters`.
`mesh_activate_gravity_method` selects `-method1` (element inactive before
activation) or `-method2` (element active with reduced stiffness via
`mesh_activate_gravity_stiffness_factor`). `control_mesh_activate_gravity_apply`
selects which records are applied.

## Uso

Place it in the data part:

```
mesh_activate_gravity_time 0  0.0 1.0
mesh_activate_gravity_element 0  -ra 1 5 -ra
```

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `mesh_activate_gravity_time` | `time_start time_end` | Global activation window. |
| `mesh_activate_gravity_element` | element range | Elements to activate. |
| `mesh_activate_gravity_element_group` | groups | Element groups to activate. |
| `mesh_activate_gravity_geometry` | geometry | Elements inside the geometry. |
| `mesh_activate_gravity_time_initial` | `time_of_birth` | Element inactive before this time. |
| `mesh_activate_gravity_time_strain_settlement` | `-yes` | Combine with `strain_settlement_parameters`. |
| `mesh_activate_gravity_method` | `-method1`/`-method2` | Activation method. |
| `mesh_activate_gravity_stiffness_factor` | factor | Stiffness reduction for method 2 (default 1e-6). |
| `control_mesh_activate_gravity_apply` | indexes | Which `mesh_activate_gravity_*` records apply. |

## Related

- `strain_settlement_parameters` — the settlement creep model combined with
  this record (`mesh_activate_gravity_time_strain_settlement -yes`).
- `force_gravity` — the gravity vector that is gradually activated.

## Estado de implementación

- **Implementado**: `mesh_activate_gravity_time` (gradual gravity activation,
  `mesh_activate_gravity_factor` in mesh.cc, applied in materi.cc). The
  element selection (`_element`/`_element_group`/`_geometry`),
  `_time_initial`, the factor interpolation, and **`_method`** (method 1 and
  2) with **`_stiffness_factor`** are implemented. Method 2 keeps the element
  active with reduced stiffness before activation (the element matrix and
  lhside are scaled by the stiffness factor in materi.cc). Validated with
  `mesh_act_grav` (method 1) and `mesh_act_grav2` (method 2).

## F4 T8 (2026-10-10): aligned with the manual and the Professional

- The records are indexed: `mesh_activate_gravity_time index ...` applies
  to the elements selected by the `_element` / `_element_group` /
  `_geometry` record with the same index (none: all elements), and
  `_method`, `_stiffness_factor`, `_time_initial` are read at that index.
  Before, only index 0 was read (the Professional `mesh_activate1` uses
  index 10: the GNU applied full gravity from the first step).
- Bottom to top: the element start (end) time is interpolated from
  `time_start`/`time_end` and the lowest (highest) vertical coordinate of
  the element within the vertical range of all elements of the record.
- Default `-method2` (manual 6.810). Before its start time a method-2
  element has no gravity and its stiffness (material tangent and stress
  increment) is multiplied by the stiffness factor; between the start and
  end times it has full stiffness and ramping gravity.
- `-method1` keeps the element in the calculation without gravity before
  its start time (the Professional takes it out; not implemented).

Professional `mesh_activate1` (202 nodes, 1000 steps): final `veliy`,
`vely` and `sigyy` = Professional to 3e-4 relative (vely top -9.9504e-6 vs
-9.95e-6; before F4 T8 vely 0) and its target passes. Test
`validation-suite/test-2026/mesh_activate_gravity_time1` (targets =
Professional values, t = 0.5 with the front half way).

## Differences with the Professional version

- The partly activated element is loaded uniformly over its height; the
  Professional nodal values are the exact column solution. GNU off by
  2.2e-9 (0.2 %) on `veliy` in the coarse test and 6 % on `vely` above the
  front; 3e-4 in `mesh_activate1` (DIVERGENCES C46d, class B).
- `-method1` does not take the element out of the calculation.

