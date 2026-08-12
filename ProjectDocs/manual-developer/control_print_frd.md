# control_print_frd

## Implementación

- **Output**: `print_frd()` in `print_fr.cc` (new file, added to the
  makefile as `PRINT_FR_OBJ`). Invoked from the control loop in `top.cc`:
  `db_active_index(CONTROL_PRINT_FRD, icontrol)` && `ival[0]!=-NO`.
- **Keywords registered** in `database.cc` (data_class CONTROL):
  `control_print_frd`, `control_print_frd_freecad`,
  `control_print_frd_prepomax`. All data_length 1.
- **New enums**: `CONTROL_PRINT_FRD`, `CONTROL_PRINT_FRD_FREECAD`,
  `CONTROL_PRINT_FRD_PREPOMAX` in `tochnog.h` / `tochnog-mod.h`.
- **Format**: CalculiX FRD, long format (FORMAT=1). The layout was taken
  from the CalculiX 2.21 source (`frd.c`, `frdheader.c`, `frdselect.c`).
  Node lines are ` -1%10ld%12.5E%12.5E%12.5E`; element lines
  ` -1%10ld%5ld%5s%5ld` + ` -2` + `%10ld` per node.
- **Element types** (frd): TRIA3=7, TRIA6=8, QUAD4=9, QUAD9=10, TET4=3,
  TET10=6, HEX8=1, HEX27=4. Connectivity reuses the same node order as
  `print_vtk` / `print_gmsh` (QUAD4 = 0,1,3,2; HEX8 = 0,1,3,2,4,5,7,6).
- **Result blocks**: per dof, `-4`/`-5` records with FRD names (DISP,
  VELO, STRESS, TOSTRAIN, NDTEMP or Tochnog name truncated to 8 chars).
  Tensor components use CalculiX indices: SXX(1,1) SYY(2,2) SZZ(3,3)
  SXY(1,2) SYZ(2,3) SZX(3,1).
- **Time in 100CL**: formatted like CalculiX (up to 9 decimals, else
  scientific).

## Diseño / decisiones

- Mesh is written only when the target file does not exist yet (for
  `-yes`), so the file grows as a time series.
- `-separate_index` / `-separate_sequential` always write a fresh file.
- Nodal coordinates are the original ones (no displacement added) for the
  `-yes` case; `node_deformed_mesh` is not applied.

## Detalles

- `db_version_copy(VERSION_NORMAL, VERSION_PRINT)` +
  `renumbering(VERSION_PRINT, NO, 0, 0, ...)`; node numbers are +1
  (FRD is 1-based).
- Matrix dofs (stress) read via `stress_indx(kdim,ldim)` — verified
  against the SQLite output of `control_print_tabular`.
- Structural elements (BAR, TRUSS, BEAM, ...) are skipped.

## Pendiente

- `control_print_frd_freecad` / `control_print_frd_prepomax` are
  registered but their switches currently do not change the output (the
  standard names are always written).
- Higher-order element connectivity (TRIA6, QUAD9, TET10, HEX27) uses the
  full node list; exact CalculiX node ordering for these was not verified
  (only linear elements tested).
- `-separate_sequential` numbering uses a static counter (resets per run).
