# control_mesh_switch

## Descripción

Permutes the coordinate axes of ALL nodes of the mesh. Useful to reorient or
rotate the mesh (for example to interchange x and y of a model built in the
wrong plane) without regenerating the geometry.

The keyword specifies the new order of the axes. Every axis of the mesh
dimension must be given exactly once (a full permutation). For example
`control_mesh_switch 0 -y -x` makes the new x-axis take the old y-coordinates
and the new y-axis take the old x-coordinates (x and y are interchanged).

## Uso

```
control_mesh_switch <index> -x -y -z
```

The axes must be given in the new order. In a 3D model all three axes must be
given (for example `-x -z -y`); in a 2D model only two (for example `-y -x`).

## Parámetros

| Parámetro | Significado |
|-----------|-------------|
| `index` | Control item index (usually `0`). |
| `-x` | New x-axis takes the coordinates of the old x-axis. |
| `-y` | New x-axis takes the coordinates of the old y-axis. |
| `-z` | New x-axis takes the coordinates of the old z-axis. |

The axis tokens describe, position by position, which OLD axis supplies the
coordinates for the NEW axis. `-y -x` therefore means: new x = old y and
new y = old x.

## Ejemplo

Interchange x and y of a 2D mesh:

```
number_of_space_dimensions 2
end_initia

node 1  0. 0.
node 2  1. 0.
node 3  1. 1.
node 4  0. 1.
element 1  -bar2 1 2

control_mesh_switch 0  -y -x
end_data
```

After this record the model runs with x and y swapped.

## Once per control index (F4 T7)

Applied once per control index, on its first step, also when the same
index carries a multi-step `control_timestep` (as the Professional; before
F4 T7 the axes were swapped on every step). See
[control_mesh_once_per_control](control_mesh_once_per_control.md).

## Two entries: swap (F4 T8)

Manual Professional 6.244: `control_mesh_switch index dir0 dir1` switches
the coordinates `dir0` and `dir1` of all nodes. Since F4 T8 two entries
always mean that swap, in 2D and in 3D: `-x -y` and `-y -x` both
interchange x and y, `-x -z` in 3D interchanges x and z, `-x -x` does
nothing (Professional 02-08-2026 probes). Before F4 T8 the GNU read the
entries as the new axis order (`-x -y` = no change, and a 3D model needed
three entries). Three entries in 3D (`-z -y -x`, the new axis order) are a
GNU extension; the Professional refuses them.

The Professional leaves `node_start_refined` unswitched; the GNU switches
it with the nodes (DIVERGENCES C46e, class D). Test
`validation-suite/test-2026/mesh_switch_pair1` (node 3 at 0 2).

