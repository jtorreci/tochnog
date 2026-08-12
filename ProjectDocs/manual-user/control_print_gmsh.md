# control_print_gmsh

## Description

`control_print_gmsh` writes the mesh and the nodal/element results in the
**Gmsh 2.2** ASCII format (`.msh` file). The `.msh` files can be opened by
Gmsh, FEniCS, Salome, meshio-based tools and other programs that read the
Gmsh format.

This is the Gmsh counterpart of `control_print_gid` / `control_print_vtk`:
use it when the results are to be viewed or post-processed in Gmsh.

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_timestep` record whose results you want to export:

```
control_print_gmsh 20 -yes
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. Must match an active `control_timestep`. |
| switch    | `-yes`: single file `<base>.msh`, mesh written only the first time, results appended per time step. `-separate_index`: file `<base><icontrol>.msh`, mesh+results written on every call. `-separate_sequential`: files `<base>0.msh`, `<base>1.msh`, ... one per call. |

## Options

- `control_print_gmsh_dummy index switch` — write a dummy point element
  (Gmsh type 15, group 1234) in every node, needed by Gmsh to plot vector
  fields at nodes. Default `-yes`; use `-no` to suppress.
- `control_print_gmsh_element_data index switch` — `-yes` (default):
  element data written as `ElementData` (one value per element, averaged
  over the element nodes). `-no`: written as `ElementNodeData` (values at
  each element node).
- `control_print_gmsh_node_method index method` — coordinates used for the
  nodes: `-node` (default, original coordinates), `-node_start_refined`
  (coordinates of the refined/initial mesh), `-node_deformed_mesh`
  (deformed coordinates when `materi_displacement` is active).

## Output

File structure (Gmsh 2.2):

- `$Nodes` / `$Elements` — the mesh (only written the first time for
  `-yes`). Element types: line, triangle, quad, tetra, hexa (linear).
- `$NodeData` blocks — one per dof component, named `node_<dof>` (scalar),
  `node_<dof>_0`, `node_<dof>_1`, ... (vector), `node_<dof>_xx`,
  `node_<dof>_yy`, ... (matrix/stress Voigt components).
- `$ElementData` / `$ElementNodeData` blocks — one per dof component,
  named `element_<dof>...` (see `control_print_gmsh_element_data`).

## Example

```
control_print_gmsh                       20  -yes
control_print_gmsh_dummy                 20  -yes
control_print_gmsh_element_data          20  -yes
control_print_gmsh_node_method           20  -node
control_timestep                         20  0.001 0.04
```

Produces `<base>.msh`; open it in Gmsh or read with
`meshio.read("base.msh")`.
