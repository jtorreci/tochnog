# prism15 (DEV-B quadratic-prism sprint)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enum value `PRISM15`, **appended at
  the end of the enum** right before `LAST_DUMMY` (same position and
  same order in both headers), keeping every existing value stable.
- `database.cc` — `name[PRISM15] = "prism15"` (registered at the end of
  `db_initialize`; the element record machinery is name-generic, no
  `type`/`data_length` registration needed).
- `check.cc` — early `return 0` for `PRISM15` in `check()`, mirroring
  `TET10`/`HEX27`.
- `polynom.cc` — `pol()` branch `name==-PRISM15` (shape functions +
  integration points) and the `volume[]` factor for PRISM15.
- `point_el.cc` — static helper `prism15_shape()` + `point_el()` branch
  (parametric Newton localization for `post_point`).

## Element definition (verified against the Professional)

Node order (corpus/Professional, confirmed from the Professional
`element_intpnt_h` of `prism15.dat`): base corners 1-3, top corners 4-6,
vertical-edge mids 7-9, base triangle edge mids 10-12 (edges 1-2, 2-3,
3-1), top triangle edge mids 13-15.

Natural coordinates: in-plane area coordinates `L1, L2` (`L3 = 1-L1-L2`,
nodes 1/2/3 of a face are the corners of `L1`/`L2`/`L3`) and `zeta` in
[-1,1] along the prism axis (`z_phys = (1+zeta)/2` for unit height).
Shape functions (serendipity 15-node wedge):

- base corners: `N = L_i (1-zeta) (2 L_i - 2 - zeta) / 2`
- top corners:  `N = L_i (1+zeta) (2 L_i - 2 + zeta) / 2`
- vertical mids: `N = L_i (1 - zeta^2)`
- base edge mids (edge i-j): `N = 2 L_i L_j (1-zeta)`
- top edge mids:  `N = 2 L_i L_j (1+zeta)`

The z-direction is quadratic Lagrange on the vertical edges only; the
`(1-zeta^2)` vertical-mid term is the serendipity correction that makes
the element complete to the quadratic in-plane x quadratic-in-z family
minus the side-face mids of the 18-node prism.

## Integration rule (measured from the Professional)

The Professional `.dbs` of `prism15.dat` stores `element_intpnt_coord`
(21 points per element): 3 Gauss levels along z (`z = 0.1127016654,
0.5, 0.8872983346` of the element height = zeta `-sqrt(0.6), 0,
+sqrt(0.6)`) x 7 in-plane points of the degree-5 Dunavant triangle rule
(area-coordinate groups: centroid; permutations of `(a,b,b)` with
`a=(6+sqrt15)/21`, `b=(9-2sqrt15)/21`; permutations of `(c,d,d)` with
`c=(6-sqrt15)/21`, `d=(9+2sqrt15)/21`). The GNU branch reproduces both
the 21-point layout and the `element_intpnt_h` values of the
Professional exactly (checked with a standalone dev harness).

Weights: triangle weights `9/80` (centroid), `(155+sqrt15)/2400` and
`(155-sqrt15)/2400` (sum 1/2 = triangle area) times zeta Gauss weights
`5/9, 8/9, 5/9` (sum 2). The natural-domain measure (triangle area 1/2
x zeta in [-1,1]) is 1, so the `volume[]` chain uses the plain
`volume = weight*detj` factor for PRISM15 (the weights sum to 1).

## Point localization (`point_el.cc`)

`post_point` needs the parametric inversion of the element. PRISM15 has
its own Newton branch over `(L1, L2, zeta)` (mirroring the generic
hex/quad Newton with `prism15_shape()` instead of the tensor-product
interpolation), with the final containment test `L1, L2, L3 in
[-eps, 1+eps]`, `zeta in [-1-eps, 1+eps]`, plus the standard distance
and convergence checks.

## Gotchas / scope decisions

- The `PRISM6` hand-written block in `polynom.cc` (legacy, closed by
  hand) is **left untouched**: no corpus PASS test feeds a real solid
  PRISM6 through `pol()`, and `prism6.dat` is a pre-existing RUNFAIL in
  the corpus. Its integration weights are not consistent with the
  `volume[]` chain of the other custom elements; fixing PRISM6 (with
  the same measured Professional rule and point_el branch) is a
  follow-up, out of the prism15 blast radius.
- `area.cc` side/border tables were NOT added: `area()` only needs
  element side tables when an area-integral record (edge/face loads,
  convection, ...) targets the element, and none of the prism15 corpus
  paths do. Adding `border_nodes_prism15` (6 sides: 2 triangles + 3
  rectangles, quadratic) is the pending step if face loads on prism15
  are ever needed.
- The enum is appended before `LAST_DUMMY` in **both** headers with the
  same comment, per the append-only database rule (inserting in the
  alphabetical section would renumber every following item).
- No exporters/refine/split changes were needed for the target test:
  the corpus run only requires the solve + `post_point_dof` print,
  which are name-generic.

## Verification

- Standalone dev harness (not committed) checked: Kronecker delta at
  the 15 nodes, partition of unity, linear reproduction (`L1, L2,
  zeta`), analytic gradients vs central differences, and the shape
  values at the Professional integration points match
  `element_intpnt_h` of `prism15.dbs` within 1e-9.
- Corpus `prism15.dat`: rc=0. GNU `post_point_dof 10` sigzz =
  1.000000261 vs Professional 1.000000000 (target 1.0, tolerance 1e-2);
  disz = 0.5000000003 vs 0.5. Both solve the uniform-traction state
  exactly.
- Blast radius: full corpus re-run on the branch — no regressions vs
  the 196-PASS baseline (only prism15 RUNFAIL -> PASS).
- Internal suite: full `build_safe.sh` run green (see its own summary).
