# bounda_geometry_method

## Description

Selects which node type is used to determine whether a node lies inside the
geometry when boundary conditions are imposed on a geometry.

- `-node_start_refined` (default) uses the refined coordinates
  (`NODE_START_REFINED`).
- `-node` uses the original coordinates (`NODE`).

Useful when the mesh has been refined and the refined coordinates differ from
the original ones.

## Usage

```
bounda_geometry_method <iboun> <node_type>
```

## Parameters

| Parameter  | Meaning                                                       |
|------------|---------------------------------------------------------------|
| `iboun`    | Boundary number imposed on a geometry.                        |
| `node_type`| `-node_start_refined` (default) or `-node`.                   |

## Example

Evaluate boundary 0 on the original node coordinates instead of the refined
ones:

```
bounda_geometry_method 0  -node
```
