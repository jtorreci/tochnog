# groundflow_phreatic_only — developer

## Files and functions

- `database.cc` — record `GROUNDFLOW_PHREATIC_ONLY` (INTEGER, length 1,
  no index, class GROUNDFLOW); enum after `GROUNDFLOW_PHREATICLEVEL_BOUNDA`
  in `tochnog.h` / `tochnog-mod.h`.
- `groundfl.cc` — `groundflow_phreatic_only_apply()`.
- `input.cc` — `input()` calls it right after `end_data`, before the
  element node check and the calculation.

## Implementation details

`-no` or no record: nothing. `-yes`: for every active `GROUP_TYPE` index
not found in any `GROUNDFLOW_PHREATICLEVEL_MULTIPLE_ELEMENT_GROUP` record,
`-GROUNDFLOW` is removed from the type list (the record is deleted and put
again with the shorter length) and every record named `group_groundflow_*`
(loop over `db_name(idat)`, `idat < LAST_DUMMY`) is deleted at that index.
An empty type list stops with the error of the user page. Any other word
than `-yes` / `-no` is `db_error`.

Measured on the Professional 02-08-2026 (coupled column of three quad4,
level 2.5 listing group 1): the `_10.dbs` has `group_type 2 -materi`, no
`group_groundflow_permeability 2` and `element_intpnt_materi_total_pressure
3` = 0; the GNU gives the same database lines and the same total pressures
(`node_dof_calcul`), with the bottom head -3 written as the GNU dynamic
pressure -0.5. The displacements already differ without the record, a
coupled-mechanics difference outside this unit.

## History

- F5 T5 U3a (2026-11-15): new.
