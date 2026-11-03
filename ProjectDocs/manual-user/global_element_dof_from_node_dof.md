# global_element_dof_from_node_dof

## Description

Input record (manual Professional 6.554):

```
global_element_dof_from_node_dof switch
```

When `global_element_dof_apply` is `-yes` (elements remember their own
strains, stresses, history variables in `element_intpnt_dof`), an element
whose `element_intpnt_dof` record is not in the input file must start its
integration-point values from somewhere:

- `-yes`: they are initialised from the `node_dof` records (interpolated
  to the integration points), e.g. an initial stress field given with
  `node_dof`.
- `-no` (default): they are not initialised from `node_dof`; they start
  at zero. A `control_reset_dof` still reaches them: the reset value is
  evaluated at every integration point of the elements whose nodes are
  all reset, as in the Professional.

The switch has no effect with `global_element_dof_apply -no` (the elements
then always use the nodal values) or when the input holds
`element_intpnt_dof` records (a restart database).

Default: `-no`, as in the manual and the Professional 02-08-2026. Before
this record existed the GNU always behaved as `-yes`: an input that gives
its initial stresses, void ratio or history variables in `node_dof` must
now write `global_element_dof_from_node_dof -yes` (the suite inputs dof1,
hypo1-4, cda1, mhardsoil_*, ... were updated so).

## Why it is useful

It decides whether an initial stress / strain field written in `node_dof`
records is taken over by the elements or is only an initial guess of the
nodal output.

## Parameters

- `switch`: `-yes` or `-no`.

## Example

`validation-suite/test-2026/global_element_dof_from_node_dof1.dat`: a quad4
with `node_dof` stresses (sxx -10..-12, syy -20..-22) and no
`element_intpnt_dof`, `global_element_dof_from_node_dof -no`: the result is
the one of the model without initial stress (node 3 sigxx 13.46).
`global_element_dof_from_node_dof2.dat` (`-yes`): the initial stresses are
taken over (node 1 sigxx -10). `global_element_dof_from_node_dof3.dat`
(record absent) gives the `-no` result. All equal the Professional
02-08-2026.

## Differences with the Professional version

- None known for the default (`-no` since C100, 2026-11-03).
- With `-yes`, a `control_reset_dof` at a later control step reaches the
  integration points only through the first-step node copy (GNU
  behaviour before C100, not compared with the Professional).
