# control_print_vtk_other

## Implementación

- Keyword `control_print_vtk_other` (data_class CONTROL, type INTEGER,
  data_length 1) registered in `database.cc` with
  `data_required = CONTROL_PRINT_VTK`.
- New enum `CONTROL_PRINT_VTK_OTHER` in `tochnog.h` / `tochnog-mod.h`
  (kept in sync).
- Read in `print_vtk()` (`print_vt.cc`) with `GET_IF_EXISTS`; default
  `-YES`. The block is written at the end of `POINT_DATA` (after the
  post_calcul fields) and contains:
  1. `SCALARS boundary_condition` — `1.0` when ANY primary dof
     (`i < npuknwn`) of the node is bounded (`node_bounded[i] != 0`),
     else `0.0`. Written only when `NODE_BOUNDED` exists.
  2. `VECTORS mesh_deformation` — the nodal displacement
     (`node_dof[dis_indx+idim*nder]`). Written only when
     `materi_displacement` is active.

## Estado: PARCIAL (documentado)

- The Professional manual (6.346) says "other things like boundary
  conditions, mesh deformation etc." without a list. The GNU subset is
  the two fields above. NOT implemented from the gid_other-style list:
  element groups as cell data, beam/truss vectors, safety slip
  surfaces, group plots — those belong to the discarded GiD family
  (control_print_gid_* DESCARTADO in SEGUIMIENTO). The user manual
  marks the feature as partial and lists what is included.

## Detalles / gotchas

- `NODE_BOUNDED` has NO `version_all` (default 0), so it is NOT copied
  to VERSION_PRINT by `db_version_copy`; it is read from
  VERSION_NORMAL, where node indices are **1-based** (node 0 does not
  exist). The print loop uses 0-based VERSION_PRINT indices (after
  `renumbering(..., lowest_node=0)`), so print index `inod` maps to
  NORMAL index `inod+1` — valid while all nodes are active (no deleted
  nodes). Documented limitation: with deleted nodes the mapping shifts.
- `NODE_BOUNDED` is created by `db_set_int` inside `bounda()` (every
  step_start with bounda records) — do NOT use `db_active_index(
  NODE_BOUNDED, 0, ...)` as the existence test (index 0 is never
  active); use `db_max_index` instead.
- `bounda_force` does NOT set `node_bounded` (only `bounda_unknown`
  does, bounda.cc), so force-loaded nodes report `0.0` — this is the
  value discriminator of the test.
- Verification: `vtk_other1` test — default `tn42.vtk` has
  `boundary_condition` (values 1.0/1.0/0.0/0.0: nodes 1,2 fixed,
  nodes 3,4 force-loaded) and `mesh_deformation`; `-no` `tn43.vtk` has
  neither.

## Pendiente

- Full gid_other list (element groups as cell data, etc.) — out of
  scope (GiD family discarded).

## Nota — sección 6.344 del manual Professional

The manual section 6.344 is a typo: it documents
`control_print_gid_group` (element groups plotted as result field)
inside the vtk section, but the keyword belongs to the gid family
(discarded, see SEGUIMIENTO "Descarte de features"). It is NOT in the
convergence inventory/checklist, so it was NOT implemented.
`control_print_vtk_group` does not exist in Professional.
