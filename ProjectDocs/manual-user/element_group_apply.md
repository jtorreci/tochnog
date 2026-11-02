# element_group_apply

## Description

`element_group_apply element group_0 group_1 ...` (manual Professional
6.429) with `control_element_group_apply index number` (6.134): at the
start of control index `index` every element with an
`element_group_apply` record gets the group at position `number` (0 =
first) of its record, and keeps it afterwards.

`control_element_group index -yes` (6.133) evaluates the
[area_element_group](area_element_group.md) records again at the start of
control index `index` (the manual advises not to use it at an index with
timesteps).

## Example (manual)

```
element_group_apply 43 1 7 4
element_group_apply 44 21 22 8
control_element_group_apply 10 0   (element 43 group 1, element 44 group 21)
control_element_group_apply 20 2   (element 43 group 4, element 44 group 8)
```

Tests: `validation-suite/test-2026/element_group_apply1`,
`control_element_group1`.

## Differences with the Professional version

- `element_group_apply` / `control_element_group_apply`: none measured.
- `control_element_group`: the Professional 02-08-2026 rejects the record
  ("I do not know : control_element_group"); the GNU follows the manual
  (DIVERGENCES `CTRL-EG`, class D).
