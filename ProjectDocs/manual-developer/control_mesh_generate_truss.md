# control_mesh_generate_truss

## Implementación

- **Output**: `generate_beam_truss()` in `generate.cc`. Called from the
  control loop in `top.cc` for each of the three element kinds:
  ```
  generate_beam_truss( icontrol, BEAM );
  generate_beam_truss( icontrol, TRUSS );
  generate_beam_truss( icontrol, TRUSSBEAM );
  ```
  Each task reads its own keyword (`CONTROL_MESH_GENERATE_BEAM`,
  `CONTROL_MESH_GENERATE_TRUSS`, `CONTROL_MESH_GENERATE_TRUSSBEAM`) and
  returns early if it is not active.
- **Keywords** (data_class CONTROL, type INTEGER, data_length 3):
  `control_mesh_generate_truss`, `control_mesh_generate_beam`,
  `control_mesh_generate_trussbeam`, registered in `database.cc`.
  `ival[0]` = element group, `ival[1..2]` = geometry entity item + index.
- **Checks** (`check.cc`): truss/trussbeam require `materi_velocity`;
  beam additionally requires `materi_displacement` (or
  `materi_velocity_integrated`), `beam_rotation` and ndim in {2,3}.
- **Neighbour logic**: for each node on the geometry, the function scans
  `NODE_NODE` (the list of nodes connected through isoparametric
  elements) and generates a truss/beam to each neighbour that is also on
  the geometry. A `generated_list` avoids duplicates.
- **Loose**: `control_mesh_generate_truss_beam_loose` (with `-yes`) makes
  the generated element use new nodes instead of connecting to existing
  ones.

## Diseño / decisiones

- `task` (TRUSS/BEAM/TRUSSBEAM) selects both the keyword to read and the
  element type to generate, so the neighbour-generation loop is shared.
- `NODE_START_REFINED` is used to know the mesh state before generation.

## Detalles

- Requires `materi_velocity` (truss) or displacement + `beam_rotation`
  (beam). Checked at input time in `check.cc`.
- Works in 1D (truss), 2D and 3D (beam requires 2D/3D).

## Pendiente

- `control_mesh_generate_truss_beam_separate` (per-group separate nodes)
  is not registered.
