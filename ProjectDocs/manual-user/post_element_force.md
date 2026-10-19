# post_element_force

## Description

`post_element_force` (Professional manual 6.927-6.935) computes the
normal force, the shear forces and the bending moments of a cross
section. The section forces are the sums of the nodal forces that the
selected elements exert on the section nodes. Use it to read the force
in a pile, a wall or a beam modelled with continuum elements, or the
soil force on a pile shaft or toe.

The cross section must lie on element sides, typically the common side
of two element layers or a boundary of the domain. A section through
the interior of elements is not allowed.

## Input

```
post_element_force index dir_normal_x dir_normal_y dir_normal_z
  dir_shear0_x dir_shear0_y dir_shear0_z dir_shear1_x dir_shear1_y
  dir_shear1_z middle_x middle_y middle_z          (3D, 12 values)
post_element_force index dir_normal_x dir_normal_y dir_shear0_x
  dir_shear0_y middle_x middle_y                   (2D, 6 values)
post_element_force index dir_normal_x middle_x     (1D, 2 values)
```

Any other number of values stops the input with
`Number of data values expected : 6` (12 in 3D, 2 in 1D), as in the
Professional.

| parameter | meaning |
|---|---|
| `dir_normal` | normal direction of the section; the nodal force components along it give `normal_force` |
| `dir_shear0`, `dir_shear1` | shear directions (perpendicular); the components along them give the shear forces |
| `middle` | exact middle point of the section; the moments use the distance from it |

The selection records, with the same index:

- `post_element_force_geometry index geometry_item_name geometry_item_index`:
  only the nodes on this geometry, taken at their INITIAL location.
- `post_element_force_group index group_0 group_1 ...`: only these
  element groups.
- `post_element_force_number index number_0 number_1 ...`: only these
  elements.
- `post_element_force_normal index -yes`: only the elements on the
  positive `dir_normal` side of the section plane. Without it, the
  elements on both sides are used if present: the result is the mean
  of the positive-side forces and the negated negative-side forces, so
  an interior section gives the section force and a load sitting on
  the section nodes is split between both sides. With elements on one
  side only (a boundary section), that side gives the result, negated
  when it is the negative side.
- `post_element_force_force` / `post_element_force_inertia index -yes`:
  add the external / inertia forces (1D, see the developer page).
- `post_element_force_multiply_factor index factor`: scale the result.

## Result

`post_element_force_result index ...` is filled at each time step:

- 3D: five values, `normal_force shear0_force shear1_force moment0
  moment1`.
- 2D: three values, `normal_force shear0_force moment0`.
- 1D: one value, the normal force (the GNU wrote the normal force
  followed by four zeros before 2026-10-17; the Professional writes one
  value).

`moment0` sums the normal force components times the distance along
`dir_shear0` from `middle`; `moment1` uses `dir_shear1`.

## Element types

The section works with every continuum element: `quad4`, `quad9`,
`quad16`, `tria3`, `tria6` in 2D, and `hex8`, `hex27`, `tet4`, `tet10`,
`prism6`, `prism15` in 3D. Each element's forces are integrated with
its own integration rule and integration-point stresses, so the same
cantilever gives the same section forces on any mesh. In an
axisymmetric group (`group_axisymmetric -yes`) the forces are those of
the full ring.

Trusses (`-truss`) and springs (`-spring1`, `-spring2`) with a node on
the section add their axial force to it, so the result is the total
force carried through the section (continuum plus structural
elements). Each one counts on the side of its centroid, like a
continuum element; a `-spring1` on a section node counts on the
positive side. Exclude them with `post_element_force_group` to get the
continuum share alone.

Other elements cannot be part of the selection: beams, truss-beams,
contact springs, bars in 2D/3D, interface element groups and `hex64`
stop the calculation with an error. Restrict the section with
`post_element_force_group` or `post_element_force_number`.

The integration-point stresses come from `element_dof`
(`options_element_dof -yes`, the default for the stress analyses).

## Example

A 2D strip of four `quad4` elements, fixed at x=0 and loaded at x=4
(velx 50 and vely 5 at both end nodes). The section at x=0:

```
geometry_line 1 0. 0. 0. 1. 1.e-5
post_element_force 10 1. 0. 0. 1. 0. 0.5
post_element_force_geometry 10 -geometry_line 1
...
target_item 40 -post_element_force_result 10 0
target_value 40 100. 1.e-4
```

The result is `100 10 -40` (normal, shear, moment), exact by
equilibrium. See `validation-suite/test-2026/post_element_force_*` and
`validation-suite/test-2014/tpef_force10`.

## Differences with the Professional version

- 2D and 3D results equal the Professional 02-08-2026 for every
  element type above, plane and axisymmetric, at boundary and interior
  sections and with a reversed normal (C84, 2026-10-17; corpus
  `force14`/`force15`).
- Trusses and springs in the selection = the Professional 02-08-2026
  (C86, 2026-10-19): a truss on the section gives 100 = truss 45 +
  continuum 55, a spring2 100 = spring 55.385 + continuum 44.615; two
  trusses at an interior section node and a 3D truss match too. The
  spring force itself still differs from the Professional when the
  spring moves off its direction (DIVERGENCES C88), and the section
  carries that spring force.
- The refusal of beams, truss-beams, contact springs and interface
  groups is a GNU check: in the 2D probe the Professional beam carries
  no force, so there is no reference for its contribution; contact
  springs and interfaces were not measured (DIVERGENCES, class B row
  C84-REFUSE).
