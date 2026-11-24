# group_integration_method -reduced

`group_integration_method index -reduced` selects the "reduced"
integration rule of the Professional manual (6.620): "a weighted
combination of lobatto and 1-point cell centered integration". The
weighting is set by `group_integration_method_reduced_factor` (6.621).
It applies to the bar, quad and hex elements of the group. Available in
the GNU version since F4 T6f C80 (2026-10-16). Before, the GNU stopped with
"I do not know : -reduced".

## Why it is useful

The factor moves the integration points between the nodes (factor 1,
Lobatto) and the element centre (large factor, 1-point rule). Points
closer to the centre soften the element. This reduces shear locking in
bending and volumetric locking for nearly incompressible materials.
Points too close to the centre leave hourglass modes free (see the
warning below).

## Syntax

```
group_integration_method  index  -reduced
group_integration_method_reduced_factor  index  factor
```

| parameter | meaning |
|---|---|
| `index` | element group |
| `factor` | f > 0. 1 means Lobatto; infinite means the 1-point rule. Default 10 when the record is absent (manual 6.621). A factor <= 0 is an input error. |

## The rule per element

| element / `group_integration_points` | rule |
|---|---|
| quad4 2x2, hex8 2x2x2 (default or `-maximal`) | 2 points per direction at iso coordinates +-1/f, equal weights |
| quad4 / hex8 `-minimal` | 1 centre point (the same as `-gauss`) |
| bar2 default | 1 point (the same as `-gauss`) |
| bar2 `-maximal` | Lobatto at the nodes (the factor is ignored) |
| bar3 / quad9 / hex27 (default or `-maximal`) | Lobatto points -1, 0, 1 (the factor is ignored) |
| bar3 / quad9 / hex27 `-minimal` | 2-point Gauss rule (unlike `-lobatto -minimal`, which uses the corners) |
| tria3, tet4, prism6 | the method is ignored (the same as `-gauss`) |
| bar4 / quad16 | Lobatto (the factor is ignored), `-minimal` Gauss (as bar3; measured on the Professional in C67) |
| hex64 | not accepted by either version |

Nodal stresses: for every `-reduced` bar, quad and hex rule, each element
contributes the linear least-squares fit `1,x(,y,z)` of its point
stresses in physical coordinates. The node takes the arithmetic mean over
its elements. This holds even when the points sit at the nodes
(f = 1, or the bar3/bar4/quad9/quad16/hex27 Lobatto points). It differs from `-lobatto`,
where each element contributes the value of its point at the node.

In an axisymmetric group, an element with a node on the symmetry axis
(radius below 1e-9) is integrated with the Gauss rule. The other elements
keep the reduced points. This holds for any factor, also when the reduced
points are off the axis.

## Warning: hourglass modes at a high factor

With a high factor the quad4/hex8 points are close to the centre and the
element stiffness approaches the 1-point rule. The hourglass modes are
then almost free. The solver gives NO warning.

The analytic study of C80 (plane-strain cantilever, thick cylinder,
point load) found the following. The errors are tip deflection / inner
radial displacement versus the exact solution.

| case | Gauss 2x2 | f = 1 | f = 2 | f = 10 (default) |
|---|---|---|---|---|
| cantilever 20x2 elements | -11 % | -46 % | -2.8 % | +31 % |
| cantilever 40x4 elements | -2.9 % | -18 % | -0.7 % | +6.3 % |
| cantilever 10x1 (one element thick) | -32 % | -77 % | -10 % | +2100 % |
| cylinder 8x16, nu = 0.49 | -4.0 % | -11 % | -3.0 % | -0.13 % |
| cylinder 8x16, nu = 0.4999 | -80 % | -92 % | -75 % | -11 % |

- A mesh one element thick in bending, or a point load, shows the
  hourglass at f = 10: zig-zag displacements next to the load.
- f of about 2 gives the best bending accuracy on square elements.
- A high f removes the volumetric locking as nu tends to 0.5.
- f = 1 is stiffer than Gauss in bending and locks volumetrically.

Check a high-f result for zig-zag displacements. Use several elements
across the thickness, and spread point loads.

## Minimal example

```
number_of_space_dimensions 2
materi_velocity
materi_stress
end_initia
node 1 0. 0.
node 2 1. 0.
node 3 0. 1.
node 4 1. 1.
element 1 -quad4 1 2 3 4
group_type 0 -materi
group_materi_elasti_young 0 1.e3
group_materi_elasti_poisson 0 0.3
group_integration_method 0 -reduced
group_integration_method_reduced_factor 0 2.
...
end_data
```

Tests: `validation-suite/test-2026/reduced_*` (quad4 f = 2 and 10,
hex8 f = 2, quad9 default and `-minimal`, axisymmetric quad4 on the axis,
and the default factor).

## Differences with the Professional version

- **Default factor** (deliberate divergence, class A, project owner
  2026-10-16; DIVERGENCES C80). Without the
  `group_integration_method_reduced_factor` record the GNU uses the
  manual default 10. The Professional 02-08-2026 uses 1 (its results
  without the record equal the f = 1 run). Give the record explicitly to
  get the same results in both versions.
- **bar nodal stresses** (BACKLOG C67): resolved 2026-10-16. The
  bar2/bar3/bar4 and quad16 `-reduced` rules are recovered with the
  physical linear fit, as in the Professional (test
  `recovery_bar3_reduced1`).
- **Very large factor** (for example 1e6): the Professional per-element
  stress is the element average (apparently a fallback when its fit is
  singular). Not reproduced or investigated.
