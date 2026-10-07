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
  5. Distribute: for each node of the found element and each direction
     `idim<ndim`, `node_rhside[vel_indx/nder+idim] +=
     weight[inol]*force_point[ndim+idim]` (guarded by the record length).
     `node_rhside` holds external − internal, so an external load is
     ADDED (same sign as the `force_element_edge` loads in `area.cc`).
     K04 `force18` fix (2026-10-06): the old code subtracted the force
     (velx −0.8 vs Pro +0.8) and looped over every principal unknown
     (`nuknwn/nder`), reading past the ndim forces of the record and
     assuming the velocity is the first principal unknown.
- Data layout: array is `[ndim coordinates][ndim force components]`
  (`data_length = ndim + MUKNWN` is only the storage size).

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
- Writes directly into `NODE_RHSIDE`. `node_force` (`dof.cc`
  `parallel_new_dof_before()`) still SUBTRACTS its value (opposite to the
  Professional; `delete.cc` compensates for mesh_cut, BACKLOG) — the
  "apply nodal rhs contribution" pattern should be one helper.
- `point_el()` locates the point iteratively (EPS_SIZE 1e-6): nodal
  weights accurate to ~1e-7 (Pro exact).
- If the force point sits exactly on a node the force is still routed
  through shape functions (weight 1 on that node) — fine, but a direct
  node path would be cheaper and clearer.

## Differences with the Professional version

- **`node_inertia` / `node_mass` layout** (GNU extension / input dialect;
  BACKLOG B12). GNU: `node_inertia` holds ndim values (m (a + g));
  `node_mass` accepted. Professional 02-08-2026: One value per dof; rejects
  `node_mass` next to it. Status: Input dialect. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
