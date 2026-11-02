# Global print defaults (F5 W2 V2) — developer

Keywords: `print_gmsh_calculation`, `print_tecplot_calculation`,
`print_vtk_calculation`, `print_gmsh_dummy`, `print_gmsh_node_type`,
`print_vtk_coord`, `print_vtk_group`, `print_vtk_node_type`,
`print_frd_freecad`, `print_frd_prepomax` (manual Professional 6.979, 6.980,
6.991-6.993, 6.998-6.1002).

## Where

- Records: enums `PRINT_GMSH_CALCULATION` .. `PRINT_FRD_PREPOMAX` right after
  `PRINT_GID_CALCULATION` in `tochnog.h` and `tochnog-mod.h` (same order);
  registered in `db_initialize` (`database.cc`) as INTEGER, length 1,
  `no_index` (one loop over the contiguous range). Aliases in `db_number`:
  `print_gmsh_node_method` -> `PRINT_GMSH_NODE_TYPE`, `print_vtk_node_method`
  -> `PRINT_VTK_NODE_TYPE` (older manual names).
- `print_option(global, control, icontrol, default)` (`miscel.cc`): global
  record, then the control record at `icontrol` overrules it; `control<0` or
  `icontrol<0` reads only the global record. Used by `print_vtk`
  (`print_vt.cc`: coord, node method, group), `print_gmsh` (`print_gm.cc`:
  dummy, node method) and `print_frd` (`print_fr.cc`: freecad, prepomax).
- End-of-run prints in `exit_tn` (`miscel.cc`, after `print_gid`):
  `print_gmsh(-1,-YES)`, `print_tecplot(ival,1)`, `print_vtk(-1)`, each
  default `-NO`.
  - `print_vtk(icontrol<0)` writes `<data_file_base>.vtk` (control prints
    keep `tn<icontrol>.vtk`).
  - `print_gmsh`: static `gmsh_yes_printed/gmsh_yes_time` remember the last
    `-yes` print; the end call returns when it equals `TIME_CURRENT`.
  - `print_tecplot(ival, end_of_calculation)` (new defaulted argument in both
    headers): statics `tecplot_printed/tecplot_time`, same skip rule; it
    appends to `tecplot.plt`.
- `print_vtk_group`: after `POINT_DATA`, a `CELL_DATA ncell` block with
  `SCALARS element_group int` from `ELEMENT_GROUP` (VERSION_PRINT, 0 when
  absent), skipping the same empty elements as `CELLS`.
- `print_frd`: freecad `-yes` returns before writing anything unless
  `(materi_displacement||materi_velocity_integrated) && materi_stress &&
  materi_strain_total`. prepomax `-yes` in 2D (`extrude`): node block of
  `2*(max_node+1)` nodes, copy id `inod+1+max_node+2` at z = 1, elements
  rewritten as penta6 (type 2, tria3/tria6 corners) or hex8 (type 1,
  quad4/quad9 corners), result blocks loop over `nnode_out` with the copy
  repeating the values; `100CL` node count = `nnode_out`.

## Latent defects fixed on the way

- `print_vtk`: the stored / `-node_start_refined` coordinate branches wrote
  `coord[idim]` for every `idim<MDIM`, but `NODE` / `NODE_START_REFINED` hold
  `ndim` values: the z of a 2D mesh was uninitialised stack memory (base
  `78ecf63`: `3.14505418671e+100`, `-5.39e-217`, ...). Now `0.0` for
  `idim>ndim-1` (pinned by `print_vtk_calculation1`).
- `print_gmsh`: same for the `-node_start_refined` branch (now `ndim` values
  plus ` 0`); no test (needs a refined mesh).

## Professional measurements (02-08-2026, scratch probes)

- No `.vtk`/`.plt`/`.msh` at the end without `-yes`; `-yes` writes
  `<name>.vtk`, appends to `<name>.plt` / `<name>.msh`; a final state already
  printed by a control print is not repeated (gm2/te2 probes: sets at 0.1,
  0.2, 0.5).
- vtk: `CELL_DATA` with `element_number` and `element_group`;
  `print_vtk_group -no` changes nothing; `print_vtk_coord -no` removes the
  `VECTORS coord` point field (POINTS kept), `control_print_vtk_coord 0 -yes`
  restores it; `print_vtk_node_type -node_deformed_mesh` deforms POINTS,
  `control_print_vtk_node_type 0 -node` overrules.
- gmsh: `print_gmsh_dummy -no` removes the 4 dummy elements, the control
  record overrules; `print_gmsh_node_type -node_deformed_mesh` deforms
  `$Nodes`.
- frd: freecad without `materi_strain_total` -> no file; with it the file
  equals the default one. prepomax extrudes 2D (quad4 hex8, tria3 penta6,
  quad9 hex8 corners with a 0 element count).

## Differences with the Professional version

- `PRINT-CALC-DEFAULT` (class A, owner 2026-11-02): tecplot / vtk end-of-run
  default `-no` = Professional, against the manual `-yes`.
- `VTK-GROUP-NO` (class D): the Professional ignores `print_vtk_group -no`.
- `C99` (class B, open): `tn<index>.vtk` / `tecplot.plt` against
  `<name>.vtk` / `<name>.plt`.

## Hardcoded / pending

- Owner decision: `print_tecplot_calculation`/`print_vtk_calculation`
  default `-NO` against the manual (`PRINT-CALC-DEFAULT`, class A).
- `tecplot.plt`/`tn<i>.vtk` names (C99 open). `print_vtk_coord` follows the
  Professional `coord` field since C98 (2026-11-02).
- Not implemented: frd `PEEQ` name (no post_calcul fields in the GNU frd);
  `element_number` vtk cell field of the Professional.
- The end-skip statics are per run (one input per process).
