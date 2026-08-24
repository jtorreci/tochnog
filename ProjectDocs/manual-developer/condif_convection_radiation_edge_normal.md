# condif_convection_edge_normal / condif_radiation_edge_normal (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — 14 enums (2 Professional masters + 6
  companions each), inserted around the legacy CONDIF_CONVECTION /
  CONDIF_RADIATION entries (headers in sync).
- `database.cc` — 14 registrations: masters DOUBLE 2 values
  (`h Tenv` / `alpha_r Tr`), `_geometry` INTEGER variable length (dual
  format: node list or geometry entity — the generic area() path
  handles both: `area[0]>0` -> node list), restriction lists
  cross-required.
- `check.cc` — all 14: `check_ndim(2,3)` + `condif_temperature` unknown
  (same as the legacy masters).
- `area.cc` — `MTYPES 7 -> 9`; `type[7]=CONDIF_CONVECTION_EDGE_NORMAL`,
  `type[8]=CONDIF_RADIATION_EDGE_NORMAL` with their `_GEOMETRY`
  companions in `type_area[]`. New helper pair
  `conv_rad_is_master`/`conv_rad_companion` (5 companions: element,
  element_group, element_side, node, element_node). The element-level
  restrictions run before the side loop (same shape as flux_edge); the
  node-level restrictions are checked in the application branch. The
  application branch itself was extended from
  `RADIATION||CONVECTION` to all four masters: values are read from
  `type[itype]` (each master holds its own record), radiation-vs-
  convection physics discriminated by master. The legacy keywords keep
  their exact previous code path (no restriction support on the legacy
  names; legacy tests convec1/convec2/condif8/condif10 stay green).

## The quad4 border-sides gotcha (found calibrating)

`border_nodes_quad4 = {0,1, 1,3, 3,2, 2,0}`: tochnog connects quad4
local nodes as (1,2) bottom, (2,4) right, (4,3) top, (3,1) left — the
Z node convention (n1,n2 bottom in x-order, n3,n4 top in the SAME
x-order), NOT counter-clockwise BL,BR,TR,TL. An "edge" that is not in
this table silently never fires (the element conduction matrix is
computed from shape functions and is correct for any ordering, which
is why mechanics tests with ccw numbering pass). Tests using area()
features MUST use the Z convention for the relevant edges.

## Verification (suite 76/76)

- `condif_convec`: analytic steady state T=0.5 (h=k=L=1, Tenv=1) with
  the Professional name.
- `condif_rad`: analytic nonlinear steady state T+T^4=1 -> 0.7245,
  10 Newton iterations (`control_timestep_iterations`).
- `condif_convec_el`: `_element` restricted to an element not touching
  the geometry -> T=0 exactly (the A/B: without the restriction the bar
  heats up).
