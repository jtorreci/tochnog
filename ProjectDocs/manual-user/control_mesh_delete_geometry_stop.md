# control_mesh_delete_geometry_stop

## Description

Stops a slow geometry deletion started by `control_mesh_delete_geometry`
together with `control_mesh_delete_geometry_factor`. All remaining delete
factors are destroyed and the partially deleted elements become fully active
again, taking their integration-point strains/stresses from the moment just
before the deletion started (the Professional default,
`global_element_dof_apply -yes`).

This is manual Professional 6.174 (and 6.175 for the geometry restriction).

## Usage

```
control_mesh_delete_geometry_stop <index> -yes
control_mesh_delete_geometry_stop_geometry <index> <geometry_item> <geometry_index>
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index. It does not need to be the same index as the `control_mesh_delete_geometry` / `_factor` records. |
| `-yes` | Stop the deletion and reactivate the elements. `-no` is a no-op. |
| `geometry_item`, `geometry_index` | Optional (`_stop_geometry`): only the elements whose nodes ALL lie inside this geometry entity are reactivated. |

## Behaviour

A `control_mesh_delete_geometry` record with a `_factor` of e.g.
`0. 0.4` deletes the elements in the geometry by 40 % over the timestep: the
elements are not removed from the mesh, they keep an
`element_delete_factor` in (0, 1). `control_mesh_delete_geometry_stop -yes`
at a later control index removes those factors, so the elements resume normal
assembly at the next step.

Elements that were deleted COMPLETELY (a `_factor` end value above
`1 - 1.e-10`, manual 6.170) have their element record removed and can NOT be
reactivated.

## Example

Delete the right bar for 40 % over the timestep at control index 10, then
reactivate it at index 20 and continue:

```
control_timestep                     10  1. 1.
control_mesh_delete_geometry         10  -geometry_line 1
control_mesh_delete_geometry_factor  10  0. 0.4

control_mesh_delete_geometry_stop    20  -yes

control_timestep                     30  1. 1.
```

## See also

- `control_mesh_delete_geometry` (manual 6.166)
- `control_mesh_delete_geometry_factor` (manual 6.170)
