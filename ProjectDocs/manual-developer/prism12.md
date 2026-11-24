# prism12 (F5 T4 U9c)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — the enum value `PRISM12` (appended by
  U9b, before `LAST_DUMMY`); prototype of `prism12_shape()`.
- `database.cc` — `name[PRISM12] = "prism12"` (registered by U9b; the
  ELEMENT record machinery is name-generic).
- `polynom.cc` — `prism12_shape(L1, L2, zeta, h, p)`: the product of the
  existing `tria6_shape()` and `interpolation_polynomial(zeta, 2, ...)`;
  `pol()` branch `name==-PRISM12` (14-point rule, `db_error` outside 3D),
  the `volume[]` factor `weight*detj`, and `element_intpnt_supported()`.
- `miscel.cc` — nodal stress recovery: `recovery_fit_kind()` kind 1 (the
  physical linear fit) for 12 nodes / 14 points, and the
  `recovery_fit_shape()` branch with the points of `pol()` and
  `prism12_shape()`.
- `point_el.cc` — `prism15_shape()` evaluates the prism12 through
  `prism12_shape()` for 12 nodes; `point_el()` runs the prism Newton
  (`L1, L2, zeta`) for `-PRISM12` and returns the 12 shape values.
- `calcul_force.cc` — `post_element_force_calculate()` integrates the
  prism12 through `pol()` (the non-tensor path of C84).
- `elem.cc` — in the `GROUP_INTERFACE` branch a `-PRISM12` stops with
  "Error: an interface -prism12 is not available yet; element N" before
  `interface_element()` (whose name dispatch would end in
  `assert( name==-HEX8 )`).
- `print_vt.cc` — VTK cell type 31 (see the VTK section).

The continuum element loop needs nothing else (`elem.cc` dispatches on
`GROUP_INTERFACE`; `area()` gives `nside = 0`, as for the prism15 and
prism18).

## Element definition

Node `6k+t` (k = 0 bottom, 1 top; t the tria6 node) has

    N(6k+t) = h6_t(L1, L2) * l_k(zeta),  l_0 = (1-zeta)/2, l_1 = (1+zeta)/2

with the tria6 functions of `tria6_shape()` (node 0 the `L1` corner,
node 2 the `L2` corner, node 5 the `L3` corner, mid nodes 1, 3, 4).
The Professional 02-08-2026 `element_intpnt_h` of a unit prism (its
-minimal run, the only rule that runs) equals the product to 6e-14 at its
8 points, and its `control_mesh_extrude_element -prism12` writes the base
tria6 then the top tria6.

## Integration rule

14 points, `pol_iso_used` = (L1, L2, zeta): the 7-point tria6 rule of the
`TRIA6` branch (its table and point order) x 2 Gauss points in zeta
(-1/sqrt(3), 1/sqrt(3), weight 1), the triangle points fastest. The tria6
weights (sum 1) are halved, so the weights sum to the natural-domain
measure 1 and `volume = weight*detj` (the PRISM15 / PRISM18 convention).
The rule is the manual 6.623 reading: the tria6 seven-point scheme of
`-maximal`, and as many Gauss points as nodes along zeta. det J is at most
of degree 4 in (L1, L2) and 2 in zeta, so the rule integrates the volume
of any prism12 exactly (`prism12_volume1`: curved base edge and sloped
top, V 2.3002).

Professional 02-08-2026: the default, `-maximal`, `-lobatto` and
`-reduced` stop with "Error detected for data item :
element_intpnt_iso_coord" on a unit prism (PRISM12-DEFAULT, class D). Its
`-minimal` (also with `-lobatto`) uses 8 points, method 0: the 4-point
tria6 rule (centroid weight -27/48) x 2 Gauss levels. A NumPy stiffness
(E 1, nu 0.25) with that rule has 3 negative eigenvalues on a regular
prism of height 0.2 and is indefinite on 71 of 200 randomly distorted
elements (amplitude 0.12, up to 4 negative eigenvalues); the 14-point rule
always has exactly the 6 rigid-mode zeros. The Professional stops on a
fixed-base prism of height 0.2 with "Error with pardiso solver. Pardiso
internal error number -4". The C82 / C132 rule applies: the GNU keeps the
14 points for every `group_integration_points` value (PRISM12-MINIMAL,
class A). The probe files are in the U9c scratch directory
(`pro/p12*.dat`, `eig12.py`, `eig12b.py`).

`element_intpnt_h` is written (12 x 14 = 168 values).

## Recovery

Kind 1 (`F = [1,x_n] (A^T A)^-1 [1,x_g]^T`, physical coordinates over the
14 points), as the prism15 and prism18. There is no Professional
reference for the 14-point rule; on the unit prism of `prism12_gauss1` the
top sigzz is -0.457 against -0.447 of the Professional 8-point run.

## Point location

The prism Newton of `point_el()` over (L1, L2, zeta in [-1,1]) with the
prism12 shape functions; final containment as the prism15. A field
quadratic in x, y and linear in z is interpolated exactly
(`prism12_post_point1`, within the Newton tolerance 1e-7). The
Professional gives other values (PRISM12-POST-POINT, class D:
u = 0.01 x^2 (1+z) at (0.25, 0.25, 0.5), (0.5, 0.25, 0.25),
(0.1, 0.6, 0.75): 5e-3, 7.5e-3, 2e-3 instead of 9.375e-4, 3.125e-3,
1.75e-4).

## Interface

The interface status comes from `GROUP_INTERFACE` (elem.cc), not from the
name. Manual 6.624 lists `-prism12` among the interface elements (two
tria6 sides) and 6.158 converts an interface `-tria6` / `-prism15` to it;
neither is implemented (PRISM12-INTERFACE, class B). The refusal sits in
`elem()`, so a `group_interface -prism12` input stops in the first
element loop.

## VTK

`print_vtk()` writes type 31 VTK_QUADRATIC_LINEAR_WEDGE through
`vtk_quadratic_cell()` (not for an interface group). VTK order: corners
(0 1 2 base, 3 4 5 top), base mid-edge (0,1) (1,2) (2,0), top (3,4) (4,5)
(5,3). GNU to VTK: `0 2 5 6 8 11 1 4 3 7 10 9`. Checked with VTK 9.7.2
reading `prism12_vtk1.vtk` and the curved `prism12_volume1` element: class
vtkQuadraticLinearWedge, 12 points, and the cell's own parametric map
(`EvaluateLocation` at four points, L1 = 1-r-s, L2 = r, zeta = 2t-1)
equal to the GNU `prism12_shape()` map to 4e-16 (script `vtkcheck12.py`
in the U9c scratch directory). The Professional drops the cell
(`CELLS 0 0`; VTK-PRISM12, class A by the owner rule 2026-11-09).

## Differences with the Professional version

- PRISM12-DEFAULT (class D): see the integration rule.
- PRISM12-MINIMAL (class A, the C82 / C132 owner rule): see the
  integration rule.
- PRISM12-POST-POINT (class D): see the point location.
- PRISM12-INTERFACE (class B): see the interface section.
- VTK-PRISM12 (class A, owner rule): see the VTK section.

## Verification

- Professional probes (02-08-2026): `element_intpnt_iso_coord` /
  `element_intpnt_h` of every rule value (only `-minimal` runs), the
  flat-prism pardiso stop, post_point in a quadratic field.
- Analytic: the self-weight sums -0.5 and -2.3002 (volume), the linear
  patch test (two layers, 5 free nodes exact, constant stress at every
  point), the quadratic-xy / linear-z field at the 14 points and at three
  post points, the section statics of `prism12_post_element_force1`.
- RED: the binary before U9c stops with "Error detected for data item :
  element" on every prism12 input and aborts on the interface input.

## Open

- BACKLOG C133: no `area()` side table (face loads), refine / split /
  mpc polishing / geometry_moving sampling do not list the element.
- The interface `-prism12` (PRISM12-INTERFACE).
