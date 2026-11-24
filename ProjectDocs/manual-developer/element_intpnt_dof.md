# element_intpnt_dof (developer)

The integration-point storage of `options_element_dof` /
`global_element_dof_apply`, renamed from `ELEMENT_DOF` in F5 T2 V5; it is
documented together with the new nodal `element_dof` in
[element_dof.md](element_dof.md).

- Enum `ELEMENT_INTPNT_DOF` (`tochnog.h` / `tochnog-mod.h`, before
  `LAST_DUMMY`), registered in `database.cc` ("element_intpnt_dof", DOUBLE,
  `npointmax*nuknwn`, `version_all`, not fixed length).
- Written by `elem.cc` (`npoint*nuknwn`, VERSION_NEW), read by `elem.cc`,
  `nonloc.cc`, `calcul_force.cc`, `delete.cc`; zero-initialised by `top.cc`,
  `input.cc`, `generate.cc`, `create.cc`; copied by `mesh.cc`; trimmed to
  the element's points by `miscel.cc` `element_dof_nodal_set()`.
