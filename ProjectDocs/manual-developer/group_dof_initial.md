# group_dof_initial — developer

## Files and functions

- `database.cc` — `GROUP_DOF_INITIAL` (DOUBLE, `DATA_ITEM_SIZE`),
  `GROUP_DOF_INITIAL_SPECIFIC_NUMBER` (INTEGER, 1),
  `GROUP_DOF_INITIAL_SPECIFIC_VALUE` (DOUBLE, `1+MDIM`); own data class
  (not ELEMENT, so element renumbering does not move them), require
  `GROUP_TYPE`. Enums after `GROUP_AXISYMMETRIC` in `tochnog.h` /
  `tochnog-mod.h`.
- `general.cc` — `general()`: one selection of the initial-dof source for
  the element.
- `top.cc` — step close: `ELEMENT_DOF_INITIAL_APPLIED` is also written for
  every active element whose group has `GROUP_DOF_INITIAL` or
  `GROUP_DOF_INITIAL_SPECIFIC_NUMBER`.

## Implementation details

Unless `ELEMENT_DOF_INITIAL_APPLIED` is set, `general()` checks level 0
(element records, index = element) then level 1 (group records, index =
`element_group`); at each level the full record, then the specific number.
The first found is used: full -> `element_dof_initial_values` (one value per
dof, the last one repeated); specific -> `edis_node_values` (linear field at
the nodes), with the missing `_value` error naming the element or the group.
The inertia substitution below is unchanged.

Measured on the Professional 02-08-2026 (bar of
`element_dof_initial_specific1`, `node_dof 2 = (1 + past)/2`): group
records alone 0.85 / 0.95; element full > group full (0.65); element full >
group specific (0.65); element specific > group specific (0.95); group full
> group specific (0.85); group full > element specific (0.85) — the last
one is DIVERGENCES DOF-INITIAL-GROUP-PRECEDENCE (GNU 0.95).

## History

- F5 T5 U4 (2026-11-17): new.
