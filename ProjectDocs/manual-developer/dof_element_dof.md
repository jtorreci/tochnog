# dof_element_dof (developer)

- Enum `DOF_ELEMENT_DOF` (`tochnog.h` / `tochnog-mod.h`, before `LAST_DUMMY`),
  registered in `database.cc` next to `DOF_TYPE`: INTEGER, length `nuknwn`
  (one per dof slot, like `DOF_LABEL`), `print_only`, `no_index`, class DOF.
- Filled in `input.cc` `input()` right after the `DOF_LABEL` / `DOF_TYPE` /
  `DOF_PRINCIPAL` PUTs, from `dof_type[]`: `-YES` when the type is in the
  `element_intpnt_dof` restore set of `elem.cc` (the
  `if(element_dof_initialised)` block: history variables, stress, elastic /
  plastic / total / intergranular strain, kappa, kappa_shear, cap1 history,
  damage, rho, Maxwell stress, softvar local / nonlocal), `-NO` otherwise.
  It does not depend on `options_element_dof` (Professional probe a_no).
- Hardcoded: the list of dof types is a copy of the `elem.cc` restore
  condition written with dof types instead of index ranges.

Verification: Professional 02-08-2026 probes f5/v5b `a_def`, `a_yes`,
`a_no` (record identical for -yes/-no; `-only` aborts in the Professional),
`g_*` (one extra dof each: elasti/plasti/total strain, work, kappa_shear,
f, rho, history variables, void). Tests `dof_element_dof1/2`
(`scripts/build_safe.sh`, `dof_element_dof_check`).

Pending refactoring: the restore condition of `elem.cc` uses index ranges
(`epe_indx` .. `+MDIM*MDIM`, which reach 3 slots past a 6-component strain
block); a shared helper `dof_is_element_dof(iuknwn)` used by both places
would remove the duplicated list.

Differences: DIVERGENCES DOF-ELEMENT-DOF-SET (class B): the Professional
marks every non-principal dof except the displacement (`work`, `f`,
`void` too).
