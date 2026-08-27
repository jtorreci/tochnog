# control_print_vtk_other

## Description

`control_print_vtk_other` controls whether "other things" (boundary
conditions, mesh deformation, etc.) are written to the `.vtk` files of
`control_print_vtk` (same `icontrol` index). `-yes` (default) writes
them, `-no` does not.

The Professional manual does not detail the list of "other things". The
GNU implements a **partial subset** (documented in
`manual-developer/control_print_vtk_other.md`):

1. `boundary_condition` — scalar per node: `1.0` when any primary dof
   of the node is bounded (`bounda_unknown`), `0.0` otherwise. A node
   loaded only with `bounda_force` is NOT bounded.
2. `mesh_deformation` — vector per node with the nodal displacement
   (the `dis` dofs); written only when `materi_displacement` is
   initialized.

What is NOT implemented (Professional `gid_other` list): element groups
as cell data, beam/truss vectors, safety slip surfaces, etc. — those
belong to the discarded GiD family.

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_print_vtk` record:

```
control_print_vtk            43  -yes
control_print_vtk_other      43  -no
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `43`      | Index of the control record. Must match a `control_print_vtk` index. |
| switch    | `-yes` (default): boundary_condition and mesh_deformation fields are written. `-no`: they are omitted. |

## Example

```
control_print_vtk            42  -yes
control_print_vtk            43  -yes
control_print_vtk_other      43  -no
```

`tn42.vtk` contains `SCALARS boundary_condition` and (with
`materi_displacement`) `VECTORS mesh_deformation`; `tn43.vtk` has
neither field.
