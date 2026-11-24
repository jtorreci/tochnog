# group_dof_initial

## Description

The `element_dof_initial` records for every element of a group (manual
6.609 - 6.611); the index is the element group:

- `group_dof_initial index dof_0 dof_1 ...`: past dofs of the elements of
  the group at their birth step (one value per dof, or one value for all).
- `group_dof_initial_specific_number index dof` and
  `group_dof_initial_specific_value index value_0 value_grad_x ...`: past
  value of one dof as the linear field `value_0 + grad . x` at the nodes
  (`value_0` alone: a constant field).

See [element_dof_initial](element_dof_initial.md) and
[element_dof_initial_specific_number](element_dof_initial_specific_number.md)
for the effect on the calculation (birth step only, inertia terms).

## Precedence

When an element has several of these records, one is used:

1. `element_dof_initial`
2. `element_dof_initial_specific_number` / `_value`
3. `group_dof_initial`
4. `group_dof_initial_specific_number` / `_value`

A record given for one element always wins over a record of its group.

A `group_dof_initial_specific_number` without its `_value` stops with
`Error: group_dof_initial_specific_value is missing for group : N`.

## Example

```
group_dof_initial 0 0.7                       ( every element of group 0 )
element_dof_initial_specific_number 1 -temp   ( element 1: this one wins )
element_dof_initial_specific_value 1 0.1 0.8
```

Self-checking inputs in `validation-suite/test-2026`: `group_dof_initial1`
(`group_dof_initial`), `group_dof_initial2` (`_specific`),
`group_dof_initial3`..`7` (precedence), `group_dof_initial_error1`.

## Differences with the Professional version

- The Professional 02-08-2026 ranks the full records before the specific
  ones: `element_dof_initial` > `group_dof_initial` >
  `element_dof_initial_specific` > `group_dof_initial_specific`. So a
  `group_dof_initial` of the group overrides an
  `element_dof_initial_specific` given for one element (0.85 instead of
  0.95 in `group_dof_initial4`). The GNU lets the element record win
  (owner decision 2026-11-17, DIVERGENCES DOF-INITIAL-GROUP-PRECEDENCE).
  Every other combination gives the Professional result.
