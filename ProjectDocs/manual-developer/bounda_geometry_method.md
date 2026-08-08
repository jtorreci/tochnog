# bounda_geometry_method

## Files and functions

- `bounda.cc` — within `bounda()`: local `bounda_geometry_method` (line 46),
  read of the keyword (lines 175–176, `GET_IF_EXISTS`) and use in the node
  loop, `use_geom` branch (lines 319–321).
- `tochnog.h` — enum `BOUNDA_GEOMETRY_METHOD` (line 141).
- `tochnog-mod.h` — mirror enum `BOUNDA_GEOMETRY_METHOD` (line 134).
- `database.cc` — keyword registration (lines 183–186).

## Implementation details

- Read: `db( BOUNDA_GEOMETRY_METHOD, iboun, &bounda_geometry_method, ddum,
  ldum, ..., GET_IF_EXISTS )`, defaulting to `0`.
- Used only when boundary conditions are applied through a geometry
  (`use_geom` branch). It selects the `node_type` argument passed to
  `geometry(...)` (`geometry.cc:25`):
  - `bounda_geometry_method == 0` → `node_type = NODE_START_REFINED`
    (default).
  - otherwise → `node_type = bounda_geometry_method` (e.g. `-NODE`).
- `geometry()` then decides, using that node type's coordinates, whether the
  node is inside the geometry (`in_geometry`).
- `database.cc`: `type = INTEGER`, `data_length = 1`,
  `data_class = BOUNDA`.

## External dependencies

None. Calls the existing `geometry()` routine.

## Hardcoded parameters / pending refactorings

- The magic default `0` overloads the keyword value: `0` means "use
  `NODE_START_REFINED`" while any other value is passed through verbatim as
  the node type. There is no validation that a `-node_start_refined` value
  (i.e. the `NODE_START_REFINED` enum constant) is actually accepted.
- The node type is read once per boundary, outside the node loop, so all
  nodes of the boundary share the same method.
