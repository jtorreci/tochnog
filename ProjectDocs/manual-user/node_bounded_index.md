# node_bounded_index

## Description

`node_bounded_index <node> index_0 index_1 ...` (manual Professional
6.879) is a print record: for every principal dof of the node (the
`node_dof` order) the index of the `bounda_dof` record that prescribes
it, or `-nothing` when the dof is not prescribed by a `bounda_dof`. Every
node gets the record (all `-nothing` for a free node). When several
`bounda_dof` records prescribe the same dof, the last applied one is
listed.

It is an output record; you do not give it in the input.

## Example

`node_bounded_index1.dat` (validation-suite/test-2026): `bounda_dof 0 1
-velx -vely` and `bounda_dof 1 4 -velx` give

```
node_bounded_index 1 0 0 -nothing ...
node_bounded_index 4 1 -nothing -nothing ...
```

## Differences with the Professional version

None found (identical records on the quad4, truss and interface probes,
Professional 02-08-2026).
