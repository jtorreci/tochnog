# geometry_node_type + geometry_projection_type

## Description

Two per-geometry records that control how a geometry entity (same
index) is evaluated when the code checks whether nodes are located on
it.

`geometry_node_type` (manual Professional 6.540) selects which node
coordinates are used for the check:

| Value | Coordinates used |
|-------|------------------|
| `-node_start_refined` | the (initial) start coordinates (default) |
| `-node` | the current node coordinates (a moving mesh is followed) |
| `-plus_displacement` | the node coordinates plus the nodal displacements |

`geometry_projection_type` (manual Professional 6.543) selects the
meaning of "on the geometry" for filled entities such as a circle:

| Value | Meaning |
|-------|---------|
| `-project_exact` | everything within the tolerance from the geometry edge (default) |
| `-project_inside` | everything inside the geometry (the filled interior) |

The records affect every consumer that evaluates nodes against the
geometry (boundary conditions, mesh refinement/deletion regions, the
`node_geometry_present` fill, ...), exactly like the Professional where
the records are properties of the geometry itself.

## Usage

```
geometry_point <index> <x> ... <radius>
geometry_node_type <index> [-node | -node_start_refined | -plus_displacement]

geometry_circle <index> <x_c> <y_c> <radius> <tolerance>
geometry_projection_type <index> [-project_inside | -project_exact]
```

## Examples

- `node_type_1.dat`: a 1D follow-material mesh (velocity 1) with three
  fixed `geometry_point`s; `geometry_node_type -node` makes the
  presence checks follow the moving nodes (see
  `node_geometry_present`).
- `validation_14.dat`: the hole circle carries
  `geometry_projection_type -project_inside` so the local refinement
  region of the `geometry_set` uses the filled interior of the circle.
