# Geometry options and entities of the 2026 manual (developer)

## Where

- `tochnog.h` / `tochnog-mod.h`: twelve `GEOMETRY_*` enums after
  `DEPENDENCY_TYPE`; `database.cc`: records of class GEOMETRY (entities
  fixed `ndim+4` / `4*ndim` / `8*ndim`, cylinder part variable);
  `check.cc`: `check_ndim` (circle part 2D, the others 3D).
- `geometry.cc`, `geometry()`:
  - entities `-GEOMETRY_CIRCLE_PART`, `-GEOMETRY_CYLINDER_PART` (start
    vector or coordinate axis; angles with `geometry_angle()` in [0, 2 pi)),
    `-GEOMETRY_TETRAHEDRAL` / `-GEOMETRY_HEXAHEDRAL` through `point_el()`
    (`-TET4` / `-HEX8`, default iso tolerance 1.e-3, the corner order of the
    elements);
  - `geometry_line_eps_iso` replaces `EPS_ISO` of the line test;
    `geometry_quadrilateral_eps_iso` sets `eps_iso` of the two triangles;
  - per-index options in the member loop (so they work for each member of
    a `geometry_set`), before the bounda sine factors: `geometry_boundary`
    (`NODE_BOUNDARY`), `geometry_mpc` (`geometry_node_has_mpc()`: first node
    of an active `MPC_NODE_NUMBER`), `geometry_element_geometry`
    (`geometry_element_geometry_filter()`, scans the elements of the node),
    `geometry_exclude` (nested `geometry()` calls with local outputs; a
    `thread_local` depth counter, limit 20, stops cycles);
  - `node_geometry_present` / `element_geometry_present` entity lists
    include the new entities in alphabetical order.
- `area.cc`, `mesh_boundary_set()`: also generates `node_boundary` when a
  `geometry_boundary` record exists; `top.cc`, `step_start()` calls it before
  the geometry-present scans (idempotent).

## Measured on the Professional 02-08-2026

2D 3x3 quad4 and 3D 2x2x2 hex8 meshes with `print_node_geometry_present`:
line/quadrilateral eps, boundary -yes/-no, exclude (point + line), element
geometry -any/-all/-only, circle part, mpc -yes (slave only) / -no (no
restriction), tetrahedral, hexahedral (prism with a skew face, inner cube),
cylinder part along x/y/z, with start vector, two ranges, reversed axis.
Pro and GNU agree except the start plane of the cylinder part and the two
input checks (user page).

## Pending

- `change_geometry` (change.cc) does not move the new entities.
- `geometry_method` (6.532): F5 T4 U4.
