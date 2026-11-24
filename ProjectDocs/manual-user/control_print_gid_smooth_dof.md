# control_print_gid_smooth_dof

Manual Professional 6.311, control_print_gid_smooth_n 6.312, print_gid_smooth_dof 6.988 and print_gid_smooth_n 6.989.

## Description

Smooths nodal dof fields in the GiD results. `-all` smooths every dof; names (initialisation
names such as `-materi_stress`, or dof labels such as `-sigxx`) smooth only those dofs. The
smoothing is done `number_of_smoothings` times (default 10). Each pass replaces the value of a
node by 0.75 times itself plus 0.25 times the mean of the other nodes of the elements around it;
a dof bounded at the node keeps its value (a displacement is bounded where its velocity is). The
smoothed displacement is also used for `mesh deform`. Only the printed values change, not the
calculation.

The global records are the defaults of the control records; `print_gid_smooth_n` has no index
(the index of the manual example is an erratum, as in the Professional).

The records apply to the GiD files of `control_print_gid` (`<base>.flavia.msh`, `<base>.flavia.res`). The mesh is written once, at the first GiD print; the results at every print. Overview of the family: [print_gid_options](../manual-developer/print_gid_options.md).

## Syntax

```
control_print_gid_smooth_dof  index  dof_0 dof_1 ...
control_print_gid_smooth_dof  index  -all
control_print_gid_smooth_n    index  number_of_smoothings
print_gid_smooth_dof  dof_0 dof_1 ...
print_gid_smooth_n    number_of_smoothings
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Index of a `control_print_gid` record. |
| `dof_i` | Initialisation name or dof label, or `-all`. |
| `number_of_smoothings` | Number of passes, default 10. |

## Example

```
control_print_gid             10  -yes
control_print_gid_smooth_dof  10  -all
control_print_gid_smooth_n    10  2
```

## Differences with the Professional version

- With `-all` the values equal the Professional 02-08-2026 (1, 2 and 10 passes).
- Professional-suspect (owner decision L7): the Professional smooths nothing when dofs are named.
  The GNU follows the manual and smooths the named dofs. Guard test:
  `control_print_gid_smooth_dof2`.

## Tests

`validation-suite/test-2026/`: `control_print_gid_smooth_dof1`, `control_print_gid_smooth_dof2`, `control_print_gid_smooth_n1`, `print_gid_smooth_dof1` (checked by `scripts/build_safe.sh` in a scratch copy).
