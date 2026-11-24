# node_geometry_present + print_node_geometry_present

## Description

`node_geometry_present` lists for every node the geometries in which
the node is present — a diagnostic record to check that the geometries
include exactly the nodes you expect (manual Professional 6.886). The
filling is switched on with `print_node_geometry_present -yes`
(6.995).

`print_node_geometry_present_node_type` (6.996) selects the default
node coordinates of the check; the per-geometry `geometry_node_type`
records override it for their own geometry (see
`geometry_node_type`).

Behavioral notes (measured on the Professional):
- The record is filled at the start of every time step with the node
  state of the converged previous step — the last step wins.
- A step without any present geometry leaves the record empty.
- Static models (no moving mesh) keep the default coordinates, so the
  "present" geometries are the ones containing the nodes at their
  start coordinates.

## Usage

```
print_node_geometry_present -yes
print_node_geometry_present_node_type [-node_start_refined | -node | -plus_displacement]

geometry_point <index> ...
geometry_node_type <index> -node
```

The output can be checked with a `target_item`:

```
target_item <n> -node_geometry_present <node> <position>
```

The stored record is a list of `(geometry name, geometry index)` pairs
per node, so position 0 is the name of the first present geometry and
position 1 its index.

## Example

See `node_type_1.dat`: a moving 1D mesh whose nodes cross fixed
`geometry_point`s; the final record marks node 1 present on
`geometry_point 3` (its position at the start of the last step), giving
`node_geometry_present (1,1) = 3`.
