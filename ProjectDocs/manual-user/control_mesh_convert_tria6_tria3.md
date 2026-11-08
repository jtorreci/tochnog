# control_mesh_convert_tria6_tria3

## Description

Converts quadratic triangles to linear ones (Professional manual 6.161):
with `-yes` each `-tria6` becomes the `-tria3` of its three corner nodes,
and the mid-side nodes that no element uses any more are deleted. Useful
to reuse a quadratic 2D mesh for a linear calculation.

## Usage

```
control_mesh_convert_tria6_tria3 index switch
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control index; the conversion runs when that index is processed (also without `control_timestep`). |
| `switch` | `-yes` converts, `-no` does nothing. |

## Example

```
number_of_space_dimensions 2
end_initia
node 1 0. 0.
...
element 1 -tria6 1 2 3 4 5 6
element 2 -tria6 3 9 7 5 8 6
control_mesh_convert_tria6_tria3 10 -yes
end_data
```

The result is `element 1 -tria3 1 3 6`, `element 2 -tria3 3 7 6` and the
nodes 1, 3, 6, 7 (the same as the Professional 02-08-2026). Test:
`validation-suite/test-2026/convert_tria6_tria3_1.dat`.

## Differences with the Professional version

None found.
