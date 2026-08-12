# control_print_gmsh

## Implementación

- **Output**: `print_gmsh()` in `print_gm.cc` (same file as `print_gmv`).
  Invoked from the control loop in `top.cc`:
  `db_active_index(CONTROL_PRINT_GMSH, icontrol)` && `ival[0]!=-NO`.
- **Keywords registered** in `database.cc` (data_class CONTROL):
  `control_print_gmsh`, `control_print_gmsh_dummy`,
  `control_print_gmsh_element_data`, `control_print_gmsh_node_method`.
  All data_length 1, `ival[0]` holds the switch/value.
- **New enum**: `CONTROL_PRINT_GMSH` .. `CONTROL_PRINT_GMSH_NODE_METHOD`
  in `tochnog.h` / `tochnog-mod.h` (kept in sync). Also added
  `NODE_DEFORMED_MESH` (keyword `node_deformed_mesh`) for the
  `-node_deformed_mesh` node method option.
- **Format**: Gmsh 2.2 ASCII. `$Nodes` are 1-based (internal 0-based
  renumbering + 1). `$Elements` use Gmsh element types (1=line, 2=tri,
  3=quad, 4=tet, 5=hex) with the same node ordering as `print_vtk`.
  Dummy point elements use type 15 and group 1234.
- **Data**: `$NodeData` / `$ElementData` / `$ElementNodeData` blocks, one
  per dof component. Names follow `node_` / `element_` prefixes (Gmsh
  convention). Uses `dof_label` + `dof_scal_vec_mat` (same detection as
  `print_vtk` / `print_tb`).

## Diseño / decisiones

- Mesh is written only when the target file does not exist yet (for `-yes`),
  so the file grows as a time series. `-separate_index` / 
  `-separate_sequential` always write a fresh file.
- Node data uses the internal 0-based index + 1 (matches renumbering with
  `lowest_node=0`); element data iterates active elements 0..max_element.
- ElementData averages the dof over the element nodes
  (`sum/nnol`). ElementNodeData writes one value per element node.

## Detalles

- `db_version_copy(VERSION_NORMAL, VERSION_PRINT)` +
  `renumbering(VERSION_PRINT, NO, 0, 0, ...)` before reading `NODE_DOF`.
- Matrix dofs (stress) are written as 6 Voigt components xx, yy, zz, xy,
  xz, yz via `stress_indx(kdim,ldim)`.
- Only linear elements (line, tria3, quad4, tet4, hex8) are written;
  higher-order elements are skipped. Truss/beam and tendon elements are
  not written.
- Hardcoded: derived magnitudes (vmises/tresca/principals) are NOT written
  by print_gmsh (they are in print_tabular/print_vtk).

## Pendiente

- Higher-order elements (tria6, quad9, tet10, hex27) not exported yet.
- `-separate_sequential` numbering uses a static counter (resets per run).
- The `$NodeData` string tag uses `node_<dof>`; Gmsh's first string tag is
  the view name.
- meshio 5.3.5 has a known incompatibility reading Gmsh 2.2 `$NodeData`
  integer tags; Gmsh 4.x reads the files correctly.
