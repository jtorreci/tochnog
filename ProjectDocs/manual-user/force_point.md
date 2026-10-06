# force_point

## Descripción

`force_point` applies a point force at an arbitrary point in space, not
necessarily on a node. Tochnog locates the element that contains the point
and distributes the force to the nodes of that element using the shape
functions. It is useful for concentrated loads at positions that do not
coincide with a mesh node.

## Uso

Place it in the data part:

```
force_point 0 x y z fx fy fz
```

In 2D, omit the z coordinate:

```
force_point 0 x y fx fy
```

## Parámetros

| Parameter | Meaning                                                                   |
|-----------|---------------------------------------------------------------------------|
| `0`       | Index of the point force record.                                          |
| `x y z`   | Coordinates of the point (ndim values: 2 in 2D, 3 in 3D).                 |
| `fx fy fz` | Force components (ndim values), external load along +x/+y/+z: a positive `fx` pushes the point towards +x. |

The point must lie inside some element of the mesh, otherwise the run stops
with an error.

## Ejemplo

2D model with a concentrated load at (1.0, 0.5):

```
control_geometry
   cartesian
control_time
   end 1.0
   dt 0.1
control_print
   history 0
   step 10
materi_elasti_young
   0 2.e7
materi_elasti_poisson
   0 0.3
materi_density
   0 2500.
materi_velocity
   0 -yes
force_point
   0 1.0 0.5 100. 0.
```

## Sign (fixed 2026-10-06, audit K04 `force18`)

Before 2026-10-06 the GNU applied the force with the opposite sign (a
positive `fx` moved the point towards −x). Now the force is an external
load along its components, as the Professional (manual 6.501). Pro test
`force18` (1D bar, force 1 at x = 0.8, node 1 fixed): node 2 `velx`
+0.8. Test `validation-suite/test-2026/force_point_sign1.dat`.

## Differences with the Professional version

The point is located iteratively (`point_el`, tolerance 1e-6), so the
nodal shares can differ from the exact shape-function values by ~1e-7
(`force_point_sign1`: disx 1.999999902 vs Pro 2).
