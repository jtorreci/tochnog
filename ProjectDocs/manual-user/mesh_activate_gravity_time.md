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
