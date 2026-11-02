# global_element_dof_apply

## Description

`global_element_dof_apply switch` (manual Professional 6.553; the GNU
also accepts its older name `options_element_dof`):

- `-yes` (default): dofs like strains, stresses and history variables are
  kept in the element integration points (`element_dof`); each element
  uses its own values. They are also averaged to the global nodes, for
  output and plotting.
- `-no`: the elements use the averaged nodal values and keep nothing of
  their own.
- `-only`: as `-yes`, but the non-principal dofs are NOT added to the
  global nodes: the nodal strains/stresses keep their initial values,
  only the integration points hold the state.

Principal dofs (velocities, temperature, pressure, ...) always live on
the nodes.

## Usage

```
global_element_dof_apply switch
```

## Example

```
global_element_dof_apply -only
```

Tests: `validation-suite/test-2026/global_element_dof_apply1` (`-yes`,
= Professional), `_apply2` (`-no`, = Professional nodal stresses),
`_apply3` (`-only`: node sigyy 0, `post_element_force` from the
integration points 3.2967).

## Differences with the Professional version

- `-only`: the Professional 02-08-2026 stops with "Program error detected
  for data item: element_intpnt_dof"; the GNU follows the manual
  (DIVERGENCES `GEDA-ONLY`, class D).
- With `-no` the GNU `post_element_force` result is 0 (it is evaluated from
  the element dofs), the Professional gives the value (BACKLOG C93).
- The `.dbs` record name is `options_element_dof`.
