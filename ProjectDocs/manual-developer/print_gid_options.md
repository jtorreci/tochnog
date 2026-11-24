# print_gid_options

GiD output options: `control_print_gid_*` (manual Professional 6.295-6.315) and the global
`print_gid_*` (6.982-6.989). User pages: `manual-user/control_print_gid_*.md` and
`manual-user/print_gid_define.md`.

## Code map

| File | Role |
|------|------|
| `tochnog.h`, `tochnog-mod.h` | enums `CONTROL_PRINT_GID_*` (after `CONTROL_PRINT_GID_EMPTY`), `PRINT_GID_*` (after `PRINT_GID_CALCULATION`, outside the `PRINT_GMSH_CALCULATION..PRINT_FRD_PREPOMAX` loop range); prototypes. |
| `database.cc` | registration. Control records: `data_class CONTROL`, `data_required CONTROL_PRINT_GID`; lists `DATA_ITEM_SIZE` with `fixed_length 0`; normals 3 doubles. Globals: `no_index`, `print_gid_smooth_dof` a list, `print_gid_smooth_n` without index (L9). |
| `miscel.cc` | `print_option()` (single value) and the new `print_option_list()` (name list), same precedence: the control record at `icontrol` overrules the global record. |
| `top.cc` | `print_gid( task, icontrol )` per control index. |
| `miscel.cc` `exit_tn` | `print_gid( -YES, -1 )` at the end of the run: only the global records apply. |
| `print_gi.cc` | dispatcher and every option helper (static state filled once per print). |
| `print_g6.cc` | the GiD writer (`print_gid_6`); ISO-8859 file, edit it with a byte/latin-1 script. |

## Flow in print_gid_6

1. `renumbering( VERSION_PRINT, YES, ... )` now fills `old_node_numbers`: the print version is
   compacted, while `node_bounded` and `mpc_node_number` live only in `VERSION_NORMAL` with the
   old numbers.
2. `print_gid_options_get( icontrol, max_node, old_node_numbers )` reads every option once
   (`print_option` / `print_option_list`; control-only records with a plain `db` read), validates
   the values (error at print time) and, for `control_print_gid_element_mpc -no`, flags the
   print nodes that are the first node of an `mpc_node_number`.
3. `print_gid_smooth()` smooths the selected dofs of `NODE_DOF` in `VERSION_PRINT` in place, so
   the dof fields and `mesh deform` print the smoothed values.
4. Mesh: `MTYPES` 12 (spring2 added); `print_gid_element_selected()` applies
   `_element_group`, `_element_mpc` and the spring2 switch; `print_coordinates()` takes its
   coordinates from `print_gid_node_coord()` (`_node_type`).
5. `print_gid_gauss_points()` declares the Gauss point sets once, when the `.flavia.res` is
   created, if any index (or `print_gid_group`) requests an element result.
6. Dof loop: a field not selected by `_dof` (and a post field not selected by `_dof_calcul`) is
   skipped by putting `outres` in the fail state (`setstate(ios::failbit)`; the stream ignores
   the output) and clearing it after `End Values`; the loop counters still advance.
   `_other -no` skips `mesh deform`.
7. After the post fields: `print_gid_coord_result()` and `print_gid_element_results()`
   (`element_group`, `truss_vector`, `beam_vectors_0..5`).

## Defaults and the existing output

Without any of these records the GiD files are byte-identical to the previous GNU (differential
check on every test-2014/test-2026 input). Hence the records that add content are off unless
they are set: `coord`, `group` (manual default `-yes`) and `spring2` (manual default 1 node) are
written only when the record (or its global form) is given. Divergence rows: parent.

## Measured Professional rules (02-08-2026)

- Smoothing: Jacobi, new = 0.75 old + 0.25 mean(neighbour nodes of the surrounding elements, self
  excluded, each node once), bounded dofs fixed; displacement and velocity_integrated fixed where
  the velocity is bounded. Matches the Professional to 1e-12 for 1, 2 and 10 passes.
- `truss_vector` = N (t x n) and `beam_vectors_k` = f_k (t x n), t unit axis, n unit normal;
  f_k = `-force_moment[k]` at the first point and `force_moment[6+k]` at the second.
- `_element_mpc -no` drops elements containing a slave (first) node of `mpc_node_number`.
- `_mesh_activate_gravity`: Professional default writes only the already activated elements, `-yes`
  every element. Not implemented (the gravity activation helpers in `mesh.cc` are static); the GNU
  always writes every element. Open question.

## Professional-suspect, GNU follows the manual

- L7: named dofs in `_smooth_dof` are smoothed (Professional: nothing). Guard
  `control_print_gid_smooth_dof2`.
- L8: named fields in `_dof_calcul` are written alone (Professional: the field is removed). Guard
  `control_print_gid_dof_calcul1`.

## Tests

27 inputs `validation-suite/test-2026/control_print_gid_*.dat` and `print_gid_*.dat`, checked in
`scripts/build_safe.sh` (block "GiD output options", helpers `gid_opt_run`, `gid_val`,
`gid_near`) after the mesh_gid checks.
