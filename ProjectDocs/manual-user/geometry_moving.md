# geometry_moving

`geometry_moving` with `geometry_moving_parameter`, `geometry_moving_operat`,
`geometry_moving_operat_parameter`, `geometry_moving_operat_time`,
`geometry_moving_n` and `control_geometry_moving` (manual Professional
6.533-6.538 and 6.135).

## Description

Moving excavations in meshes that cannot follow the excavation zones (for
example a tunnel boring machine through complex geology). A geometrical
entity appears at `start_time`, moves through the mesh with a constant
velocity until `end_time` and then stays where it is. The elements it passes
are excavated:

- an element that was completely inside the entity at some moment is
  removed (the nodes no element uses any more go too);
- an element that was only partly inside is excavated for the same part:
  it keeps `1 - part` of its released stress as a load and has no stiffness
  (the `element_delete_factor` behaviour of
  `control_mesh_delete_geometry_factor`, see
  [delete_factor_element_matrix_scaling](delete_factor_element_matrix_scaling.md)).

The excavation of an element is gradual: its delete factor goes linearly
from 1, at the moment the first part of the element enters the entity, to
`1 - part`, at the moment the last part enters. The excavation is
permanent: the part is everything that was inside the entity at any moment,
also when the entity has moved on.

`control_geometry_moving index -initialise` determines, once and at the
start of that control index, which part of every element the entities will
excavate and when. Put it before the `control_timestep` records of the
excavation.

## Records

| Record | Meaning | Default |
|--------|---------|---------|
| `geometry_moving index entity` (6.533) | The entity: `-geometry_point` (a circle / sphere with a radius), `-geometry_triangle`, `-geometry_quadrilateral` (2D, 3D), `-geometry_hexahedral`, `-geometry_tetrahedral` (3D). | - |
| `geometry_moving_parameter index ...` (6.534) | The start coordinates of every vertex (1 for a point, 3 triangle, 4 quadrilateral, 4 tetrahedral, 8 hexahedral; `ndim` values each, in the node order of the similar element), then for `-geometry_point` its radius. | required |
| `geometry_moving_operat index -translate` (6.535) | How the entity moves: `-translate`, a constant velocity. | required |
| `geometry_moving_operat_parameter index vx [vy [vz]]` (6.536) | The velocity, one value per space direction. | required |
| `geometry_moving_operat_time index start_time end_time` (6.537) | The entity exists from `start_time`, moves until `end_time` and then stays. | required |
| `geometry_moving_n index ntime nspace` (6.538) | `nspace` sample points per direction of each element (`nspace`, `nspace^2`, `nspace^3` points for 1D, 2D, 3D elements) are checked at `ntime+1` moments from `start_time` to `end_time`. | `100 5` |
| `control_geometry_moving index -initialise` (6.135) | Determine the excavation of all elements by all `geometry_moving` entities. | - |

The tests use the geometry tests of `geometry()` (a point is inside a
`-geometry_point` at most at the radius; a `-geometry_hexahedral` accepts
the small isoparametric tolerance of `geometry_hexahedral`). Several
entities may be given; an element takes the entity that excavates the
largest part.

The result is printed per element:

- `element_geometry_moving element index` — the entity;
- `element_geometry_moving_part element part` — the excavated part;
- `element_geometry_moving_in_time element t0 t1` — when the first and the
  last part of the element enter (the moment before the sample point is
  found inside);
- `element_geometry_moving_out_time element t0 t1` — when the first and the
  last part leave again (only when some part leaves).

## Usage

```
geometry_moving                  10 -geometry_point
geometry_moving_parameter        10 -1. -1. 0.2
geometry_moving_operat           10 -translate
geometry_moving_operat_parameter 10 1. 1.
geometry_moving_operat_time      10 0. 3.
geometry_moving_n                10 200 20
...
control_geometry_moving          25 -initialise
control_timestep                 30 1.e-3 2.
```

## Example

`validation-suite/test-2026/geometry_moving_ramp1.dat`: the two quad9 of
the Professional `delete3` test with initial stress `sigyy 1`. A
quadrilateral `x 0..2, y 5.2..6` moves down with velocity `0 -1` from time
0 to 1 and sweeps `y 4.2..6`: 0.8 of the height 2 of the top element, a part
0.4. The stress of the bottom element settles at `1 - 0.4 = 0.6`;
`geometry_moving_ramp2` stops at time 0.6, inside the ramp, with `0.76`.

```
geometry_moving 1 -geometry_quadrilateral
geometry_moving_parameter 1  0. 5.2  2. 5.2  0. 6.  2. 6.
geometry_moving_operat 1 -translate
geometry_moving_operat_parameter 1  0. -1.
geometry_moving_operat_time 1  0. 1.
geometry_moving_n 1  8 5
control_geometry_moving 6 -initialise
control_timestep 10 0.2 1.
```

Other examples: `geometry_moving_band1` (a point sweeps a band, part 0.4 =
analytic), `geometry_moving_band2` (complete excavation deletes elements
and nodes), `geometry_moving_simplex1` (tria3), `geometry_moving_hex1`
(3D hex8 and tet4), `geometry_moving1` (the Professional example).

## Differences with the Professional version

- **Sample points** (DIVERGENCES GEOMETRY-MOVING-SAMPLING, class D): the
  manual checks `nspace` points per direction inside each element; the
  Professional 02-08-2026 checks the element nodes only (the part is the
  fraction of nodes ever inside, `nspace` has no effect). The GNU follows
  the manual: a point that passes between the nodes excavates the element
  (`geometry_moving_band1`: Professional nothing, GNU 0.4), and the
  `delete3` stress of `geometry_moving_ramp1` is 0.6 (analytic) where the
  Professional gives 0.667 (3 of the 9 nodes).
- **Record lengths** (GEOMETRY-MOVING-PARAMETERS, class D): the GNU takes
  the manual lengths. The Professional needs an extra trailing value in
  `geometry_moving_parameter` of a `-geometry_quadrilateral` and as many
  `geometry_moving_operat_parameter` values as vertices x dimensions, and
  then gives wrong leave times (`-1.4e20`).
- The manual writes `-geometry_tethrahedral`; both programs use
  `-geometry_tetrahedral` (the Professional rejects the manual spelling).
- `control_geometry_moving` with another value than `-initialise` is an
  error in the GNU; the Professional accepts `-yes` without effect.
- The Professional example `geometry_moving1` stops with threads ("Program
  error ... element_rhside_delete ... allocated in a parallel loop"); it
  runs with `processors_maximum -no`. The GNU runs the input as given.
