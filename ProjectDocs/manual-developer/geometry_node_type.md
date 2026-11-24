# geometry_node_type + geometry_projection_type

## Implementation

- Keywords registered in `database.cc`:
  - `GEOMETRY_NODE_TYPE` (INTEGER, one value, class GEOMETRY, per
    geometry index).
  - `GEOMETRY_PROJECTION_TYPE` (INTEGER, one value, class GEOMETRY, per
    geometry index).
  - `PROJECT_INSIDE` (INTEGER pure name entry for the `-project_inside`
    keyword value).
  Enums appended in `tochnog.h`/`tochnog-mod.h` (same enum, same
  order).
- Consumption in `geometry()` (`geometry.cc`), at the head of the
  per-entity loop:
  - `GEOMETRY_NODE_TYPE[index]` (stored negated, like every keyword
    value) overrides the caller-provided node_type and RE-READS the
    node coordinates (NODE / NODE_START_REFINED records, or NODE +
    displacement dofs for PLUS_DISPLACEMENT) before the membership
    test. Works per member of a `geometry_set` (each member keeps its
    own override).
  - `GEOMETRY_PROJECTION_TYPE[index]` overrides the projection type:
    `PROJECT_INSIDE` joins the "filled interior" semantics of the
    delete/cut projection types. The 8 filled-vs-edge decisions of the
    entity branches now test a local `project_inside` flag (set once
    per entity from the effective projection type).
- Absent records keep the caller-provided node_type/projection_type
  unchanged: the blast radius is limited to models that enter the new
  records (previously none parsed).
- Verified against the Professional binary 25-10-2023:
  - `node_type_1.dat`: `-node` follows the moving 1D nodes; the final
    `node_geometry_present` value equals the Professional (node 1 on
    geometry_point 3).
  - `validation_14.dat`: the hole circle with `-project_inside` runs
    the whole adaptive refinement + delete + restart analysis rc=0.

## Pending

- `-project_inside` semantics for entities without an interior/edge
  distinction (lines, points) are handled by their own branches; only
  the circle/sphere/cylinder-like "filled" decisions use the flag.
