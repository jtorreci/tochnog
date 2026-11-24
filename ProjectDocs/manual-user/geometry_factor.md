# geometry_factor

## Description

Spatial weighting factors for loads applied on a geometry entity
(manual Professional 6.527). The factor multiplies the load value of
`bounda_dof`, `bounda_force` and `force_edge_*` records whose entity
(`-geometry_line`, `-geometry_triangle`, `-geometry_quadrilateral`,
`-geometry_point`) has the same index.

This is the Professional name of the GNU record `geometry_bounda_factor`
(the GNU implementation already reads it in `geometry()`); registering
the alias makes Professional input files work unchanged.

## Syntax

```
geometry_factor index factor_0 factor_1 [factor_2 ...]
```

The number of factors depends on the geometry entity:

| Entity | Factors | Variation |
|--------|---------|-----------|
| `geometry_line` | 2 | linear between the start and end of the line |
| `geometry_line` | 3 | parabolic (values at start, middle, end) |
| `geometry_triangle` | 3 | linear across the three corners |
| `geometry_quadrilateral` | 4 | linear across the four corners |
| `geometry_point` | 1 | half sine wave peaking at the point, zero at the tolerance distance |

## Example

From the Professional manual 6.527: a line from (0,0) to (1,0) with
factors 1 and 4, and a temperature boundary of 20. Node 2 at x=0.2
gets temperature `20 * (1 + 0.2*(4-1)) = 20 * 1.6 = 32` and node 3 at
x=0.4 gets `20 * 2.2 = 44`:

```
geometry_line 1 0. 0. 1. 0. 0.01
geometry_factor 1 1. 4.
bounda_dof 0 -geometry_line 1 -temp
bounda_time 0 0. 20. 1.e6 20.
```

## Tests

Corpus Professional tests unlocked (rc=0): `matrix2` (heat, linear
factor per side, middle point target temp 2.5), `temp2` and `matrix4`
(same linear-per-side pattern). Verified against the Professional
binary (user-supplied 25-10-2023): GNU `post_point_dof` = 2.4999998123
vs Professional 2.500000000000e+00 on `matrix2`.

## Notes

- The factors are evaluated with the LOCAL coordinate of the projection
  of the node on the entity (xi for a line), not with the global
  coordinates.
- `delete3`/`delete2` use the related `control_mesh_delete_geometry_factor`
  record, which is a different feature (delete-geometry weighting) and
  is NOT covered by this alias.
