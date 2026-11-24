# mesh_element_group_apply

## Description

`mesh_element_group_apply index group_0 group_1 ...` (manual Professional
6.817) and `control_mesh_element_group_apply index group_0 group_1 ...`
(6.178): only the elements of the listed element groups are evaluated, in
all timesteps (`mesh_`) or in the timesteps of control index `index`
(`control_mesh_`, which overrules `mesh_` in its index). The other
elements keep their state; nodes that only belong to them are dormant (no
boundary conditions, values kept), as for a `group_time` group that is not
alive. The groups of several `mesh_element_group_apply` records are
united. Default: every group is evaluated.

## Usage

```
mesh_element_group_apply          0 0 2
control_mesh_element_group_apply 10 1
```

Tests: `validation-suite/test-2026/mesh_element_group_apply1` / `2`,
`control_mesh_element_group_apply1` / `2`.

## Differences with the Professional version

The Professional 02-08-2026 contradicts the manual (DIVERGENCES
`MESH-EG-APPLY`, class D; the GNU follows the manual):

- it reads `mesh_element_group_apply` without index (all values are
  groups; a second record is an error), so `mesh_element_group_apply 0 1`
  evaluates groups 0 and 1;
- it applies the selection only at the start (control index 0) and then
  deletes the elements of the other groups for the whole calculation;
  `control_mesh_element_group_apply` at a later index has no effect.

`mesh_element_group_apply 0 0` (group 0 only, all timesteps) gives the
same results in both programs.
