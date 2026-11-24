# groundflow_phreatic_only

## Description

`groundflow_phreatic_only switch` (manual 6.581, no index). With `-yes`,
the groundflow data is removed for the element groups that are not listed
in a `groundflow_phreatic_level_multiple_element_group` record, so only
the groups that have a multiple phreatic level keep their groundflow.
Default `-no`.

## Effect on the calculation

Before the calculation starts, for each group that no multiple level
lists:

- `-groundflow` leaves its `group_type` (`-materi -groundflow` becomes
  `-materi`);
- its `group_groundflow_*` records are deleted.

Those elements then carry no pore pressure (`materi_total_pressure` 0 at
their integration points), as in the Professional 02-08-2026, which edits
the database the same way (the `.dbs` shows the new `group_type`).

A group whose only type is `-groundflow` would be left without a type:
that is an input error, `Error: group : N has only -groundflow and is not
in groundflow_phreatic_level_multiple_element_group`. The Professional
also stops there, with `Error detected for data item : group_type`.

## Example

```
groundflow_phreatic_level_multiple 0 2.5
groundflow_phreatic_level_multiple_element_group 0 1
groundflow_phreatic_only -yes
group_type 1  -materi -groundflow
group_type 2  -materi -groundflow     ( becomes -materi )
```

Self-checking inputs in `validation-suite/test-2026`:
`groundflow_phreatic_only1` (`-yes`), `groundflow_phreatic_only2` (`-no`),
`groundflow_phreatic_only_error1`.

## Differences with the Professional version

- The error for a group left without a type names the group and the
  reason; the Professional only names `group_type`.
- In the coupled test the total pressures are equal; the mechanical
  results already differ without this record (coupled column with a
  non-static multiple level), so the comparison holds for the database
  edit and the pressures only.
