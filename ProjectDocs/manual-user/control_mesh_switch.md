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
