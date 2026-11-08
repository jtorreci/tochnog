# dependency_item and its modifiers

`dependency_geometry`, `dependency_method`, `dependency_number`,
`dependency_type` (manual Professional 6.404-6.408), plus the
`element_group -all` option of `dependency_item`.

## Description

`dependency_item index -group_xxx element_group dof n` with
`dependency_diagram index dof_0 ... dof_n-1 value_0 ...` makes a real
`group_*` record of `element_group` depend on a dof, a post calculation
result, a space coordinate or `-time_current` through a multi linear
diagram. The modifiers, with the same index, refine it:

| record | values | meaning | default |
|---|---|---|---|
| `dependency_item ... -all ...` | element_group `-all` | the diagram is used for every group that has the record itself | — |
| `dependency_method index method` | `-use` / `-multiply` | `-use`: the diagram value replaces the record value; `-multiply`: it multiplies the record value | `-use` |
| `dependency_type index type` | `-cosinus` / `-sinus` / `-tangent` | work on the cosinus/sinus/tangent of the value: with `-use` the arc-function of the diagram value is used; with `-multiply` the diagram value multiplies the cosinus/sinus/tangent of the record value and the arc-function of the product is used | the value itself |
| `dependency_number index number` | integer or `-all` | only value `number` of the record (0 for the first) is dependent | `-all` |
| `dependency_geometry index -geometry_xxx gindex` | geometry | only elements with all their nodes inside the geometry use the dependency (`geometry_method` of that geometry: `-any` / `-average`, see [element_geometry_selection](element_geometry_selection.md)) | no geometry selection |

Diagram rows: the diagram has either one row of data values, used for
every dependent value of the record, or one row per value of the record.
With `dependency_number` the row of that value is used (the only row when
there is one).

The record that is made dependent must be specified itself (for
`-multiply`, and for `-all`, which only acts on groups that have the
record).

## Example

Lower the tangent of the Mohr-Coulomb friction angle of every group, for
the elements around a point, between time 10 and 11 (manual 6.406):

```
geometry_point 10 0. 0. 5.
dependency_item 1 -group_materi_plasti_mohr_coul -all -time_current 2
dependency_number 1 0
dependency_method 1 -multiply
dependency_type 1 -tangent
dependency_diagram 1 10. 11. 1. 0.
dependency_geometry 1 -geometry_point 10
```

Angles of the group records are radians.

`print_group_data` shows the values computed with the diagrams (see
[element_print_group_data_values](element_print_group_data_values.md)).

Errors:

- `dependency_number` beyond the number of values of the record, or not an
  integer / `-all`: `Error detected for data item : dependency_number`.
- A diagram with more than one row and a row count different from the
  number of values of the record: `Error detected for data item :
  dependency_diagram`.
- `-sinus` / `-cosinus` with a value outside [-1, 1]: `Error detected for
  data item : dependency_type`.
- `-multiply` on a record that is not specified for the group: `Error
  detected for data item : dependency_method`.

Tests (validation-suite/test-2026): `dependency_method1`, `_method2`,
`dependency_type1`-`3`, `dependency_number1`, `_number2`,
`dependency_geometry1`, `dependency_number_error1`,
`dependency_diagram_error1`, `dependency_type_error1`.

## Differences with the Professional version

- Time of the `-time_current` diagram: the GNU evaluates it at the end of
  the step, the Professional 02-08-2026 in the middle of the step; with
  `-total_linear` the Professional also sums the stress increments of the
  steps, each with the E of its step (two steps of 1 on a diagram
  1000 -> 300 by t = 1: GNU sigxx 6 = E(2) eps, Professional 9.5). BACKLOG
  C122, DIVERGENCES DEPENDENCY-TIME.
- A `-sinus` / `-cosinus` value outside [-1, 1] is an input error in the GNU;
  the Professional stops with a NaN error (DIVERGENCES DEPENDENCY-TYPE-RANGE).
