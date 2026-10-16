# Nodal stress recovery

The stress reported at a node (`node_dof`, `post_point_dof` and the stress
fields written to `.flavia.res`) is the ABSOLUTE weighted average of the
per-element nodal stresses of the adjacent elements:

```
sigma_node = sum_e (element stress of e at the node) / (number of adjacent elements)
```

This is the Tochnog Professional convention (each element contributes weight
1, the per-element value is a least-squares fit of its integration-point
stresses). Before 2026-09-15 the GNU version accumulated the element stress
INCREMENTS instead, so at a node shared by two elements whose increments
cancel the reported stress kept its previous value.

The same arithmetic mean applies to every element type. Since 2026-10-11
(F4 T6f) this includes the `tria3` (default 1 point and `-maximal`), the
`tria6`, the `tet4` (1 point) and the 1-point `bar2`, which before were
AREA (volume, length) weighted: at a node shared by a small stiff element
and a large soft one the reported stress is now the plain mean of the two
element values, as in the Professional. Example:
`validation-suite/test-2026/tria3_nodal_recovery_mean1.dat` (two tria3 of
area 0.5 and 2 sharing an edge): node 2 `sigxx` 6.4643 = (1.0938 +
11.8348)/2 (before: 9.6866). The strain tensors follow the same rule.
The bar rules with more than one point follow the Professional since
2026-10-16 (see "Bars and quad16" below), and so do the tet4 `-maximal`,
prism15 and prism6 `-minimal` rules (see "The last rules" below).

## quad4 and quad9 per-element values (K21, 2026-10-12)

Since F4 T6f K21 the per-element nodal value of the 2D Gauss rules is the
one of the Professional 02-08-2026:

- `quad4` with 2x2 Gauss points (default) and `quad9` with 2x2 points
  (`group_integration_points -minimal`): a least-squares plane
  `a + b x + c y` through the integration-point values, fitted in the
  PHYSICAL coordinates and evaluated at the node. On rectangles and
  parallelograms this is the same as before; on distorted elements the
  nodal stresses change (up to tens of per cent of the local value).
- `quad4` with 1 point (`-minimal`): the point value, so the node is the
  plain mean of the adjacent element values (before: area weighted).
- axisymmetric `quad4` element with a node exactly on the axis (r = 0):
  the shape-function weighted average of its point values. Elements that
  only come close to the axis (r = 1e-6) use the plane fit.

Example: `validation-suite/test-2026/recovery_tube1_axis1.dat` (laminar
tube flow, axisymmetric quad4): node 3 sigxx -159.58 (before -66.40),
equal to the Professional. Other examples:
`recovery_quad4_distorted_plane1.dat`, `recovery_examp1_quad4_point_mean1.dat`,
`recovery_incnav4_quad9_minimal1.dat`. The strain tensors and the other
recovered integration-point states follow the same rule, also when the
stress is CONVECTED (`convection_apply -yes` on a fixed mesh with
`materi_velocity`; since 2026-10-12, K21 WU1b, no exception).

Convected stress (C72, 2026-10-13): on a mesh fixed in space with
`convection_apply -yes`, the nodal stress is the same absolute recovery
of the integration-point stress, as in the Professional. Neither program
convects the integration-point stress: it stays where it is while the
material flows through the mesh. Example:
`validation-suite/test-2026/convected_stress_recovery1.dat` (sigzz pulse,
uniform velocity 1): the pulse stays at 0.5 / 1.0 / 0.5, equal to the
Professional. Before, the GNU transported the nodal stress at a wrong
speed (about v times the integration-point volume). Use a Lagrangian
mesh (`options_mesh -follow_material`) when the stress history of a
solid matters (BACKLOG C73).

Difference with the Professional version: on the Eulerian `validation_1`
flow the Professional velocity drifts in time although the problem is
steady; the GNU stays steady (DIVERGENCES C72-V1, class B). Other
Professional-side residuals of the convected path are listed in BACKLOG
C72.

## The other Gauss rules and the 3D elements (C70, 2026-10-13)

Since F4 T6f C70 the same per-element rule covers every quad/brick Gauss
rule and the tet10 / prism6 rules, as in the Professional 02-08-2026:

- `quad9` 3x3 (default), `hex8` 2x2x2 (default), `hex27` 3x3x3 (default)
  and 2x2x2 (`-minimal`), `tet10` (5-point rule) and `prism6` (6 points):
  the least-squares linear field `a + b x + c y (+ d z)` through the
  integration-point values, in PHYSICAL coordinates, evaluated at the
  node. Before, quad9 / hex8 / hex27 / tet10 used the same fit in the
  element's natural coordinates, which is identical on parallelograms,
  parallelepipeds and straight-edged tet10 but not on curved or distorted
  elements (typically 2-15 % of the local stress); prism6 used an area
  weighted average.
- `hex8` with 1 point (`-minimal`): the point value, as the quad4.
- axisymmetric `quad9` (3x3 and 2x2) with a node exactly on the axis:
  the shape-function weighted average, as the quad4.

Examples (`validation-suite/test-2026`, targets = Professional):
`recovery_hex8_distorted1.dat`, `recovery_hex8_point_mean1.dat`,
`recovery_hex27_distorted1.dat`, `recovery_tet10_curved1.dat`,
`recovery_prism6_linear1.dat`, `recovery_quad9_curved1.dat`,
`recovery_quad9_axis1.dat`.

Difference with the Professional version: the `prism6` integration
points themselves differ (the Professional samples the triangle corners,
the GNU three interior points; a deliberate divergence decided on
2026-10-13 because the GNU rule is more accurate in bending, BACKLOG
C75), so prism6 nodal stresses agree only where the stress field is
linear. The quad4/quad9/hex8/hex27 Lobatto rules follow the
Professional since 2026-10-13 (BACKLOG C79), the bar and quad16 rules
since 2026-10-16 (below), and the tet4 `-maximal`, prism15 and prism6
`-minimal` rules since 2026-10-16 too ("The last rules" below).

## Bars and quad16 (C67, 2026-10-16)

The 1D continuum bars (`bar2`, `bar3`, `bar4` with
`number_of_space_dimensions 1`) and the `quad16` follow the same rules as
the quad and hex elements, measured on the Professional 02-08-2026 (65
distorted probes of unequal element size, every node):

| Rule | Per-element nodal value |
|------|-------------------------|
| Gauss, any number of points (bar2 `-maximal`, bar3, bar4, quad16 default and `-minimal`) | least-squares line `a + b x` (quad16: plane `a + b x + c y`) in physical coordinates, evaluated at the node |
| `-reduced` (Lobatto points, the factor is not applied) and Lobatto `-minimal` | the same physical linear fit |
| Lobatto `-maximal` (the points are the nodes) | the point value |
| 1 point (bar2 default) | the point value |

The node is the arithmetic mean over the adjacent elements; an
axisymmetric quad16 with a node on the axis uses the on-axis rule of the
quad4/quad9. Before, these rules used the shape-function weights and a
length/area-weighted mean (bar3 Lobatto `-minimal`: the mid nodes were 0).
Truss and beam elements carry no recovered nodal stress (their results are
`element_truss_force` and `element_beam_force_moment`), and the `hex64`
is not accepted by either version.

Examples (`validation-suite/test-2026`, targets = Professional):
`recovery_bar2_gauss2_mean1.dat`, `recovery_bar3_gauss_fit1.dat`,
`recovery_bar3_lobatto_minimal1.dat`, `recovery_bar3_reduced1.dat`,
`recovery_bar4_gauss_fit1.dat`, `recovery_bar4_lobatto1.dat`,
`recovery_quad16_gauss_fit1.dat`, `recovery_quad16_lobatto1.dat`,
`recovery_quad16_lobatto_minimal1.dat`,
`recovery_quad16_axisym_lobatto1.dat`.

Minimal input (bar3, 3 elements, every node prescribed):

```
number_of_space_dimensions 1
materi_velocity
materi_stress
end_initia
node 1 0.25
...
element 1 -bar3 1 2 3
...
group_type 0 -materi
group_materi_viscosity 0 1.0
group_integration_method 0 -lobatto
group_integration_points 0 -minimal
```

## The last rules: tet4 `-maximal`, prism15, prism6 `-minimal` (C67, C76, 2026-10-16)

With these three every element and integration rule of the GNU follows the
Professional recovery (81 Professional probes of the tria3, tria6, tet4,
tet10, prism6 and prism15, plane and axisymmetric, every
`group_integration_points` and `group_integration_method` value; every
node compared). Before, they used the shape-function weights and a
volume-weighted mean; for the prism15 that average was unbounded at the
corner nodes (the serendipity corner functions are negative inside the
element), giving nodal stresses many times the real ones.

| Element and rule | Integration points | Per-element nodal value |
|------------------|--------------------|-------------------------|
| tet4 `-maximal` | 5 (Keast: centroid and 4 interior points), as manual 6.623 and the Professional; before: 4 points at the nodes | least-squares plane in physical coordinates (the tet4 stress is constant, so the point value) |
| prism15 (every value) | 21 (7 triangle points x 3 levels) | least-squares plane `a + b x + c y + d z` in physical coordinates |
| prism6 `-minimal` | 2 (triangle centroid at the two Gauss levels, weight 1/4 each); before: the record was ignored | `sum_g h_n(g) sigma_g / sum_g h_n(g)` (a weighted mean of the two point values) |

The node is the arithmetic mean over the adjacent elements, as for every
other element.

Two integration fixes come with this:

- **prism15**: two of the seven triangle points were wrong (one lay
  outside the triangle), so the 21-point rule did not even integrate a
  linear field exactly. The corrected rule is the standard degree-5
  Dunavant rule, the one the Professional uses. Displacements of an elastic
  prism15 block moved by up to 3.7 % and now equal the Professional to
  3e-11.
- **prism6 `-minimal`** is now honoured. It is an under-integrated rule
  (one spurious mode per element, restrained by the mesh in every case
  tried). It passes the patch test and is more flexible in bending than
  the default 6-point rule: for a cantilever 10x1x1 under an end load
  (Timoshenko), the tip error is -18.0 / -7.1 / -2.8 % on meshes
  10x1x1 / 20x2x2 / 40x2x4, against -34.7 / -13.1 / -4.85 % for the
  default rule. The default, `-normal` and `-maximal` keep the 6-point rule
  (C75).

Examples (`validation-suite/test-2026`, targets = Professional):
`recovery_tet4_maximal1.dat`, `recovery_prism15_fit1.dat`,
`recovery_prism6_minimal1.dat`, `prism6_minimal_cantilever1.dat`,
`prism15_elastic_shear1.dat`.

Minimal input (prism6, 2-point rule):

```
number_of_space_dimensions 3
materi_velocity
materi_displacement
materi_stress
end_initia
...
element 1 -prism6 1 2 3 4 5 6
...
group_type 0 -materi
group_materi_elasti_young 0 1.0
group_materi_elasti_poisson 0 0.3
group_integration_points 0 -minimal
```

Differences with the Professional version (integration rules, not the
recovery; BACKLOG C66, C82, C83):

- `tria6`: the GNU integrates with 6 points at the nodes and mid-edges
  (and refuses axisymmetric groups). The manual (6.623) and the
  Professional use 7 points (default and `-maximal`) or 4 points
  (`-minimal`) and recover them with the physical linear fit. Nodal
  stresses differ by 8-10 % of the scale on curved, distorted meshes
  (BACKLOG C66).
- prism15 `-minimal`: the Professional uses 8 points (the 4-point
  triangle rule with a negative centroid weight x 2 levels). That rule
  gives an indefinite element stiffness (negative eigenvalues), so the GNU
  keeps 21 points; nodal stresses differ by 0.9 % of the scale
  (deliberate divergence, owner decision 2026-10-16; BACKLOG C82).
- tria3 `-maximal` and tet10 `-maximal` (deliberate divergence from the
  MANUAL, owner decision 2026-10-16; BACKLOG C83): the manual (6.623)
  specifies 4 and 10 points without giving them; the GNU and the
  Professional both use the 3 nodes and the 5 Keast points, so the
  results agree. The classic 4-point triangle rule has a negative
  centroid weight (as the prism15 case above).

Axisymmetric tria3 `-maximal` (since 2026-10-16, = Professional; BACKLOG
C83): 3 points at the nodes, recovered as the plane tria3 (point value at
the node, arithmetic mean); an element with a node on the axis (r = 0)
uses the 1-point centroid rule instead, because the hoop term of a point
on the axis is infinite. Before, the GNU integrated every axisymmetric
tria3 with 1 point (nodal stresses 3 % of the scale off). The default and
`-minimal` stay 1 point in both versions.

## When it matters

The difference is visible in models where a node is shared by elements with
opposite stress increments, in particular partially deleted / reactivated
elements (`control_mesh_delete_geometry_factor` +
`control_mesh_delete_geometry_stop`). Example: `delete2.dat` (two `bar2` in
series, right bar deleted 40 % then reactivated, `control_reset_dof -sigxx 1`)
now reports the shared node `2` with `sigxx = 0.8`, matching the Professional;
before it reported `0.6`.

A currently deleted element is EXCLUDED from the recovery: it contributes
neither the numerator nor the denominator, so a node whose only adjacent
element is deleted keeps its previous stress and a shared node is recovered
from the active elements alone (`delete3.dat`: `sigyy = 0.6` on the shared
nodes, `1.0` on the nodes exclusive to the deleted element).

## Scope

The absolute average is applied when the stress is a purely recovered field.
Two cases keep the historic incremental update because the stress there is
not a pure post-processing recovery:

- the stress is CONVECTED (e.g. `strconv1`, a `materi_strain_total` transport
  model with `options_mesh -fixed_in_space`);
- an independent strain unknown exists (`materi_strain_total` and the
  plasti/elasti strain variants).

These cases have no Professional reference (their keyword sets are not
supported by the Professional binary) and are unchanged.

## Residual

`mpc5` still reports its stress one iteration behind the velocity (1.0 at
`k=1`, 0.5 at `k=2`; the Professional reports 0.5 at `k=1`). The test passes
at its own iteration count; closing the one-iteration lag needs a post-solve
constitutive re-evaluation (see the developer manual).
