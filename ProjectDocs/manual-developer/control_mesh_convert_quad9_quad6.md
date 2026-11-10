# control_mesh_convert_quad9_quad6 (and the continuum quad6)

## Files and functions

- `mesh.cc` — `mesh_convert_quad9_quad6()`: called at the start of
  `mesh_convert_tria6_tria3()` (same `step_start()` lifecycle, own control
  index); shares `mesh_delete_unused_nodes()` with the tria6 conversion.
- `database.cc` — `CONTROL_MESH_CONVERT_QUAD9_QUAD6` (INTEGER, CONTROL; enum
  appended at the end in `tochnog.h` / `tochnog-mod.h`, F5 T4 U9a).
- `polynom.cc`, `pol()` — continuum `-QUAD6`: npol 3 along xi, 2 along eta
  (`nnol_eta`), per-direction `npol_dir`; `element_intpnt_supported()`.
- `miscel.cc` — `recovery_tensor_npol()`, `recovery_tensor_nodes_dir()`,
  `recovery_fit_kind()`, the QUAD6 branch of `recovery_fit_shape()`.
- `point_el.cc` — point location (3 x 2 tensor).
- `area.cc`, `area()` — `force_edge`-type side loads: the 3-node sides
  (`border_nodes_quad6_xi`) and 2-node sides (`border_nodes_quad6_eta`) are
  walked in two passes (`goto quad6_second_pass`), gated on GROUP_INTERFACE.
- `calcul_force.cc` — `post_element_force` takes the quad6 by the `pol()`
  path (`pef_pol`); an interface group stays refused.

## Implementation details

**Element.** `name==-QUAD6` means a continuum element only when the group has
no `GROUP_INTERFACE`: `elem.cc` dispatches interface groups before `pol()`.
The tensor path is reused with `nnol_eta=2`: node order xi fastest (3 x 2),
identical to the interface layout (pairs i, i+3).

**Integration rule** (Professional 02-08-2026, measured from
`element_intpnt_iso_coord` / `_method` / `_npoint` / `_h`): default 3 x 2
Gauss; `-minimal` 2 x 1 Gauss (`npoint_per_dir = npol_dir-1`); `-lobatto`
3 x 2 Lobatto; `-reduced` 3 x 2 Lobatto (npoint_per_dir==npol_dir), and with
`-minimal` the 2 x 1 Gauss. The Lobatto -minimal case is 2 x 1 Lobatto (xi at
the ends, eta at 0).

**Recovery** (`recovery_fit_kind`): the default 3 x 2 Gauss and the -reduced
3 x 2 Lobatto use the physical linear fit (kind 1); Lobatto at the nodes the
point value (kind 3); the 2 x 1 Gauss rule (-minimal, -reduced -minimal)
`F = h_n / sum h_n` (kind 4, the prism6 -minimal rule; measured from the
Professional `element_intpnt_dof` -> `node_dof`). The 2 x 1 Lobatto rule keeps
the `h_inol` fallback (its 2 points cannot carry the linear fit; also matches
the Professional on the probe). Every rule compared on a distorted element:
node_dof equal to the Professional within 5e-14 (`quad6_minimal_recovery1`).

**Rank deficiency of -minimal.** 2 x 1 points give rank 6 for 9 deformation
modes: 3 zero-energy modes per element. Accepted as the Professional
(`quad6_minimal1`: the free mid nodes come out equal to the Professional);
a mesh relying on them for stability is ill-posed, like the quad9 -minimal
(rank 12 of 15).

**Converter.** `-x` keeps columns 0 and 2 (quad9 tensor n0 n3 n6 n2 n5 n8 =
`el[1],el[4],el[7],el[3],el[6],el[9]`); `-y` keeps rows 0 and 2. The element
is cleared with `db_delete_index` before the PUT: a shorter record left a
stale tail after the sentinel which `db_len`'s bisection reads as data (a -x
element with length 10 failed in `renumber.cc` "non-existing node"). A
remaining `dir` other than `-x`/`-y` is `db_error` (Professional: the same
two lines, `Error detected for data item`). The Jacobian of a `-x` element is
negative (`pol()` uses `|detj|`, as the Professional: element_volume positive).

Verified against the Professional 02-08-2026 (tests below: every nodal dof,
the integration rule records and the converter's element/node records match
on the probe files; the analytic checks are quad6_shape1 and
quad6_post_element_force1).

## Pending

- Element-class data of a converted element is not resized (as the tria6
  conversion; for mesh set-up, before the calculation).
- Unsupported for the quad6 (no side / tensor data): `tendon.cc`,
  `delete.cc` mesh cut, `geometry_moving.cc`, `mpc.cc` element polishing,
  `post.cc` post_force_edge, `element_normal_sides` (BACKLOG C131).
- `mesh_convert_tria6_tria3()` still shrinks its record in place without
  `db_delete_index`; its test passes, but the stale tail is the same hazard.

## Differences with the Professional version

- `-quad6` in a 3D calculation runs in the Professional; the manual lists
  it in 2D only and the GNU refuses it (DIVERGENCES QUAD6-3D).
