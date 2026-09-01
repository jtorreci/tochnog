# mpc_linear_quadratic

## Description

Automatic tying of non-compatible solution fields between linear and
quadratic elements (manual Professional 6.873).

When a quadratic element (quad9, hex27, bar3, tria6, tet10, ...) borders
a linear element (quad4, hex8, bar2, tria3, tet4, ...) at a common
interface, the EXTRA nodes of the quadratic element (mid-edge, mid-face)
are not attached to the linear element: they dangle and the solution
field becomes non-compatible at the interface. With
`mpc_linear_quadratic -yes` Tochnog generates a multi point constraint
for every dangling node: the node is tied to the linear element with the
shape functions of the linear element evaluated at the node position, so
the quadratic field follows the linear field at the interface.

Typical application (manual): stiff structural parts (beams, sheet
piles, tunnel shells) modeled with quadratic elements, surrounded by
linear soil elements.

## Syntax

```
mpc_linear_quadratic -yes
```

A switch only: no index, no arguments.

## Generated constraints

For every node of a quadratic element that is not attached to any linear
element AND lies inside a linear element, one `mpc_node_number` +
`mpc_node_factor` record pair is generated per principal dof:

- mid-edge node -> 2 masters (the edge endpoints of the linear element),
  factors 0.5, 0.5
- mid-face node -> 4 masters (the face corners), factors 0.25 each
  (3D only)

The generated records are regenerated automatically whenever the mesh
changes (refinement, deletion, splitting, ...).

## Example

```
element 1  -quad9 1 2 3 5 6 7 9 10 11
element 3  -quad4 9 11 13 14

mpc_linear_quadratic -yes
```

The quad9's mid-edge node 10 (1,2) lies on the edge 9-11 of the quad4
and is tied: `velx_10 = 0.5*velx_9 + 0.5*velx_11` (same for `-vely`).
Corpus tests: `mpc3`, `mpc4`.

## Tests

The generated records are verified byte-for-byte against the
Professional (mpc3/mpc4/mpc5). The corpus tests `mpc3`/`mpc4` still
report RUNFAIL: the GNU's staggered mixed u-sigma solve does not reach
the exact homogeneous field (within +-1e-3) that these tying tests
require (mpc3 0.299 vs 0.333, mpc4 0.350 vs 0.333; the no-tie fields are
also off) — the remaining gap is the solver family (DIAG-SOLVE-MIXTO),
not the tie generation. `mpc5` additionally needs the
`control_mesh_delete_geometry_factor` deletion family (element
half-deletion + stress reset semantics) — documented as PENDING, see
SEGUIMIENTO-CONVERGENCIA.md.

## Notes

- The generated records are visible in the .dbs database exactly like
  user records (verified against the Professional 25-10-2023 output:
  identical indices 0..N, identical master lists and factors on
  mpc3/mpc4/mpc5).
- The tolerance on the isoparametric coordinates below which a dangling
  node is considered to be inside the linear element is 1.e-4 (the same
  default as `mpc_element_group_eps_iso`, manual 6.865).
