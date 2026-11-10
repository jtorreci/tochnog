# prism18 (F5 T4 U9b)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enum value `PRISM18` appended before
  `LAST_DUMMY` (same position in both headers); prototype of
  `prism18_shape()`.
- `database.cc` — `name[PRISM18] = "prism18"` next to `PRISM15` (only the
  word; the ELEMENT record machinery is name-generic).
- `polynom.cc` — `prism18_shape(L1, L2, zeta, h, p)`: the product of the
  existing `tria6_shape()` and `interpolation_polynomial(zeta, 3, ...)`;
  `pol()` branch `name==-PRISM18` (21-point rule), the `volume[]` factor
  `weight*detj`, and `element_intpnt_supported()`.
- `miscel.cc` — nodal stress recovery: `recovery_fit_kind()` kind 1 (the
  physical linear fit) for 18 nodes / 21 points, and the
  `recovery_fit_shape()` branch with the points of `pol()` and
  `prism18_shape()`.
- `point_el.cc` — `prism15_shape()` evaluates the prism18 through
  `prism18_shape()` for 18 nodes; `point_el()` runs the prism Newton
  (`L1, L2, zeta`) for `-PRISM18` and returns the 18 shape values.
- `calcul_force.cc` — `post_element_force_calculate()` integrates the
  prism18 through `pol()` (the non-tensor path of C84).
- `print_vt.cc` — VTK cell (see the VTK section).

The continuum element loop needs nothing else (`elem.cc` dispatches on
`GROUP_INTERFACE`; `area()` gives `nside = 0`, as for the prism15).

## Element definition

Node `6k+t` (k = 0 bottom, 1 middle, 2 top; t the tria6 node) has

    N(6k+t) = h6_t(L1, L2) * l_k(zeta)

with the tria6 functions of `tria6_shape()` (node 0 the `L1` corner,
node 2 the `L2` corner, node 5 the `L3` corner, mid nodes 1, 3, 4) and the
quadratic Lagrange functions of the nodes zeta = -1, 0, 1. Manual 6.413
(element list), its side table 9.16 and the element figure give this
layered order; the Professional 02-08-2026 `element_intpnt_h` of a unit
prism equals the product to 9e-14 (probe in the scratch directory of the
U9b session).

## Integration rule

21 points, `pol_iso_used` = (L1, L2, zeta): the 7-point tria6 rule of the
`TRIA6` branch (its table and point order) x 3 Gauss points in zeta
(-sqrt(0.6), 0, sqrt(0.6), weights 5/9, 8/9, 5/9), the triangle points
fastest. The tria6 weights (sum 1) are halved, so the weights sum to the
natural-domain measure 1 and `volume = weight*detj` (the PRISM15
convention). Professional 02-08-2026: the same coordinates and order for
the default, `-maximal`, `-lobatto` and `-reduced` (method 0).

`-minimal` keeps the 21 points (PRISM18-MINIMAL): the Professional uses 8
points (`tri4` of the TRIA6 branch, centroid weight -27/48, x 2 Gauss
levels). A NumPy stiffness of a randomly distorted prism18 (E 1, nu 0.25)
has 6 negative eigenvalues and 12 zero ones with that rule, 6 zero ones
(the rigid modes) with the 21 points. The Professional itself stops on a
single fixed-base element ("Error with pardiso solver. Pardiso internal
error number -4"). The rule mirrors C82 (prism15); owner confirmation
pending (BACKLOG C132).

`element_intpnt_h` is not written (18 x 21 = 378 values > DATA_ITEM_SIZE),
as for the prism15.

## Recovery

Kind 1 (`F = [1,x_n] (A^T A)^-1 [1,x_g]^T`, physical coordinates over the
21 points). Professional nodal stresses of `prism18_gauss1`,
`prism18_volume1` and the patch test agree to 3e-10 (its iterative solver
noise).

## Point location

The prism Newton of `point_el()` over (L1, L2, zeta in [-1,1]) with the
prism18 shape functions; final containment as the prism15. The
Professional interpolates `post_point` values linearly between the
corner nodes (PRISM18-POST-POINT, class D: u = 0.01 z^2, z 0.25 -> 2.5e-3
instead of 6.25e-4).

## Differences with the Professional version

- PRISM18-MINIMAL (class B, pending the owner): see the integration rule.
- PRISM18-POST-POINT (class D): see the point location.

## Verification

- Professional probes (02-08-2026): element_intpnt_iso_coord /
  element_intpnt_h of every rule value; node_dof and node_rhside of
  `prism18_gauss1`, `prism18_volume1`, `prism18_patch1`; section result
  of `prism18_post_element_force1`.
- Analytic: volume 1.825 of the distorted element, the self-weight sum,
  the linear patch test (5 free nodes exact, constant stress at every
  point), post_point in a quadratic field.
- RED: the binary before U9b stops with "Error in data part." on every
  prism18 input.

## Open

- BACKLOG C133: no `area()` side table (face loads), refine / split /
  mpc polishing / geometry_moving sampling do not list the element.
