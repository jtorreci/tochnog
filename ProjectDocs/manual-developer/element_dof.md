# element_dof + element_intpnt_dof (developer)

## Implementation (F5 T2 V5, owner decision 2026-11-02, option 2)

- `tochnog.h` / `tochnog-mod.h`: enum `ELEMENT_INTPNT_DOF` appended before
  `LAST_DUMMY`; prototype `element_dof_nodal_set()`.
- `database.cc`: `ELEMENT_INTPNT_DOF` "element_intpnt_dof" takes the former
  `ELEMENT_DOF` properties (DOUBLE, `data_length` `npointmax*nuknwn`,
  `version_all`, class/required ELEMENT) and is now not fixed length.
  `ELEMENT_DOF` "element_dof" is the nodal record: DOUBLE, `data_length`
  `MNOL*nuknwn`, `version_all`, `print_only`, not fixed length.
- Every former `ELEMENT_DOF` site (the integration-point state of
  `options_element_dof`) now uses `ELEMENT_INTPNT_DOF`: `elem.cc` (GET
  VERSION_NORMAL/NEW, PUT VERSION_NEW), `nonloc.cc`, `calcul_force.cc`
  (materi_stress_force / post_element_force sources), `delete.cc` (mesh
  cut stress source), `create.cc`, `mesh.cc`
  (`mesh_duplicate_element_state`), `generate.cc`, `top.cc` and `input.cc`
  (zero initialisation). `ELEMENT_DOF_INITIALISED` keeps its name and pairs
  with the integration-point record.
- `elem.cc`: both PUTs write `npoint*nuknwn` values (the element's points),
  not `npointmax*nuknwn`. `nonloc.cc` keeps separate GET lengths for the
  `ielem` and `jelem` records (a GET overwrites its length argument).
- `miscel.cc`, `element_dof_nodal_set()`: called by `top()` right after the
  step-end `db_version_copy( VERSION_NEW, VERSION_NORMAL )` (serial). Per
  element with `ELEMENT_INTPNT_DOF` and `ELEMENT_INTPNT_NPOINT`, not an
  interface: `pol()` (start-of-step coordinates, shape check disabled as
  in `element_intpnt_set()`; it also sets the thread-local
  `stress_recovery_*` flags the fit reads), then the operator
  `F(inol,ipoint)` from `recovery_simplex_factor()` (new static, shared
  with `sri_stress_recovery_weight()`) or `recovery_fit_factor()`; principal
  dofs (`DOF_PRINCIPAL >= 0`) take `NODE_DOF` VERSION_NORMAL. Writes
  `ELEMENT_DOF` VERSION_NORMAL (`nnol*nuknwn`). Gated by
  `options_element_dof==-YES` (otherwise the integration-point record is
  not maintained).
- The same function trims an `ELEMENT_INTPNT_DOF` record still at the
  `npointmax*nuknwn` length of the zero initialisation: `db_len()` finds the
  first trailing `DBL_MAX` by bisection, so a shorter PUT into a longer
  record is not enough; the record is copied, `db_delete_index()`'d and
  PUT again.
- The same function deletes the (zero) record of the structural elements
  (`-beam`, `-contactspring`, `-spring1/2`, `-truss`, `-trussbeam`, which
  `elem()` dispatches before the integration points); `elem.cc` reads the
  record with GET_IF_EXISTS (zero buffers when absent).
- `create_element()` excludes `ELEMENT_DOF` from the generic copy (the
  nodal values of the old element do not apply to the new one; the next
  step end refills them).

## Hardcoded choices

- Fallback when no fit applies (or it is singular): `F = h_n(g) V_g /
  sum_g h_n(g) V_g`, plain mean `1/npoint` when that sum vanishes.
- SRI quad4/hex8 (GNU-only) take the plain fit, not the Lagrange mix of
  their lumped nodal stress.

## Tests

- `test-2026/element_intpnt_dof1`, `element_dof_nodal1`,
  `element_dof_nodal_1pt1`, `element_dof_nodal_quad9_1`
  (`beam_rec_check` in `scripts/build_safe.sh`), Professional 02-08-2026
  probes intq / two / p1 / q9, every slot within 1e-11 except disx/disy.

## Pending refactorings

- The displacement is stored in the integration-point record (the
  Professional writes 0); dropping it needs a check of every reader.
- Interface and empty elements keep the `npointmax`-point zero record (the
  zero initialisation does not know the rule).
- The old `.dbs` layout (`element_dof` = integration points) is not
  recognised on input.
