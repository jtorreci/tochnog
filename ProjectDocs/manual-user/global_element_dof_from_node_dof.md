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
- `-no`: they are not initialised from `node_dof`; they start at zero.

The switch has no effect with `global_element_dof_apply -no` (the elements
then always use the nodal values) or when the input holds
`element_intpnt_dof` records (a restart database).

Default in the GNU: `-yes` (the historical GNU behaviour). The manual and
the Professional default is `-no`; see below.

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
taken over (node 1 sigxx -10). Both equal the Professional 02-08-2026.

## Differences with the Professional version

- Default: the manual and the Professional 02-08-2026 use `-no`; the GNU
  keeps `-yes`, its behaviour before the record existed (about 200 suite
  inputs give `node_dof` records in the input). Open owner decision
  ELEMENT-DOF-FROM-NODE-DOF-DEFAULT (DIVERGENCES, class B). Write the
  record explicitly to get the same result in both programs.
