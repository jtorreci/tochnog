# control_print_vtk_node_method

## Implementación

- Keyword `control_print_vtk_node_method` (data_class CONTROL, type
  INTEGER, data_length 1) registered in `database.cc` with
  `data_required = CONTROL_PRINT_VTK`.
- New enum `CONTROL_PRINT_VTK_NODE_METHOD` in `tochnog.h` /
  `tochnog-mod.h` (kept in sync).
- Read in `print_vtk()` (`print_vt.cc`) with `GET_IF_EXISTS`; default
  `-NODE_START_REFINED` when the record is absent. The `POINTS` loop
  selects the coordinates:
  - `-NODE_DEFORMED_MESH` + `materi_displacement`: stored coordinates +
    the displacement dof (`coord[idim]+node_dof[dis_indx+idim*nder]`).
  - `-NODE_START_REFINED` with an active `NODE_START_REFINED` record for
    the node: the start-refined coordinates (pattern `print_gm.cc:369`
    and `print_g5.cc:210`).
  - otherwise: the stored `node` coordinates.
- `NODE_START_REFINED` is written by the mesh-adjust/projection
  machinery (adjust.cc, area.cc `PROJECT_EXACT`) and by the restart
  copy in top.cc; for plain models it does not exist and the default
  falls back to the stored coordinates.

## Decisión (documentada, manual 6.345)

- Manual default is `-node_start_refined` (gid family 6.307 states it
  explicitly; the vtk section does not repeat it). Since
  `-node_start_refined` falls back to the stored coordinates when no
  `NODE_START_REFINED` record exists, the default is output-identical
  to the previous behavior for the whole test suite (no model with
  refinement writes vtk).
- **Semantic change (documented)**: BEFORE this feature, `print_vt.cc`
  wrote the DEFORMED coordinates whenever `materi_displacement` was
  active, unconditionally. NOW the deformed coordinates are written
  only with `-node_deformed_mesh`; the default (`-node_start_refined`
  -> stored coordinates) does NOT deform the mesh. This matches the
  manual (deformed coordinates are an explicit option, not a side
  effect of `materi_displacement`).
- Verification: `vtk_nodmeth1` test — Total Lagrange model
  (`group_materi_memory 0 -total`, the only memory that allows
  `materi_displacement`, stress.cc check) compressed with
  `vely=-0.01`; `-node` gives stored y=1, `-node_deformed_mesh` gives
  y=0.998 (disy=-0.002).

## Pendiente

- Nothing.
