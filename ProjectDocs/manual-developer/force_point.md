# force_point

## Archivos y funciones

- `force.cc` → `force_point_calculate()` (`force.cc:229`) — implements the
  point-force distribution.
- `point_el.cc` → `point_el()` (`point_el.cc:27`) — locates the element
  that contains a given point and returns the shape-function weights.
- `dof.cc` (`dof.cc:39`) — calls `force_point_calculate()` in the same
  block as `force_gravity_calculate()`.
- `database.cc:1628-1632` — keyword registration:
  `strcpy(name[FORCE_POINT], "force_point")`, type `DOUBLE_PRECISION`,
  `data_length = ndim + MUKNWN`, `fixed_length = 0`,
  `data_class[FORCE_POINT] = FORCE_POINT`.
- Enum `FORCE_POINT` in `tochnog.h` / `tochnog-mod.h`.

## Detalles de implementación

- `force_point_calculate()`:
  1. `db_max_index( FORCE_POINT, max_force, VERSION_NORMAL, GET )` — bail
     out early if no point force or no element exists (`force.cc:241-245`).
  2. For each active `FORCE_POINT` record, read its array with
     `db_dbl( FORCE_POINT, iforce, VERSION_NORMAL )` (`force.cc:249`).
  3. Loop over all elements, gather nodal coordinates from
     `NODE_START_REFINED`, and call `point_el( force_point, coords, weight,
     el[0], nnol )`; break on first hit (`force.cc:252-261`).
  4. If no element contains the point → error + `exit(TN_EXIT_STATUS)`
     (`force.cc:264-267`).
  5. Distribute: for each node of the found element,
     `node_rhside[ipuknwn] -= weight[inol]*force_point[ndim+ipuknwn]`
     for each non-zero force component (`force.cc:269-276`). The force is
     subtracted because it is written directly into `NODE_RHSIDE`.
- Data layout: array is `[ndim coordinates][MUKNWN force components]`.

## Dependencias externas

- `point_el()` from `point_el.cc` (shared with `map.cc`, `contact.cc`,
  `post.cc`). Core database accessors (`db()`, `db_dbl()`,
  `db_active_index()`, `db_max_index()`).

## Parámetros hardcodeados / refactorizaciones pendientes

- The element search is a linear scan over all elements per force point
  (`force.cc:252`); with many elements or forces this is O(nel × nforce).
  A spatial index (e.g. bounding-box or binning) would scale it.
- `point_el()` has known pre-existing bugs for truss/beam elements and for
  some quad4 positions; a point force located in such an element may be
  misdistributed or not found. Fixing `point_el()` fixes this feature too.
- Writes directly into `NODE_RHSIDE` (as `db_dbl` then in-place
  subtraction) — works for the explicit formulation, but the sign and the
  `nuknwn/nder` loop are implicit. Document or abstract the
  "apply nodal rhs contribution" pattern shared with other load types.
- If the force point sits exactly on a node the force is still routed
  through shape functions (weight 1 on that node) — fine, but a direct
  node path would be cheaper and clearer.
