# control_mesh_generate_interface

## Implementación

- **Generación**: `generate_interface()` in `generate.cc` (added to the
  makefile, invoked from `step_start` in `top.cc` BEFORE the `any_interface`
  scan so the interface histories are allocated for the generated elements).
- **Keywords** (data_class CONTROL) registered in `database.cc`:
  - `control_mesh_generate_interface` (INTEGER, variable length,
    `fixed_length=0`): the record is a list of triples
    `eg_i eg_a eg_b eg_j eg_c eg_d ...`.
  - `control_mesh_generate_interface_geometry` (INTEGER, length 2, required
    `CONTROL_MESH_GENERATE_INTERFACE`): `geometry_item_name
    geometry_item_index` — restricts generation to the given geometry.
  - `control_mesh_generate_interface_method` (INTEGER, length 2, required
    `CONTROL_MESH_GENERATE_INTERFACE`): `method_select method_generate`.
    `-element_geometry` selects between elements by `element_geometry`
    (instead of `element_group`) and/or generates an `element_geometry`
    record for the interface element.
  - `element_geometry` (INTEGER, length 1, data_class ELEMENT): assigns a
    geometrical set number to an element; elements with the same number
    form a geometry referenced as `-element_geometry N`.
- **New enums**: `CONTROL_MESH_GENERATE_INTERFACE`,
  `CONTROL_MESH_GENERATE_INTERFACE_GEOMETRY`,
  `CONTROL_MESH_GENERATE_INTERFACE_METHOD`, `ELEMENT_GEOMETRY` in
  `tochnog.h` / `tochnog-mod.h` (kept in sync).

## Algoritmo

For each triple `(eg_i, eg_a, eg_b)` in the record:

1. Scan element pairs `(iel, jel)` with `iel` matching `eg_a` and `jel`
   matching `eg_b`. The match is by `element_group` by default, or by
   `element_geometry` when `method_select==-ELEMENT_GEOMETRY`.
2. Count the shared-face node pairs: nodes of `iel` and `jel` with
   coincident coordinates (`NODE_START_REFINED`, tolerance `EPS_COORD`).
3. Choose the generated element type from the number of shared nodes:
   - 2D, 2 shared nodes → `-quad4` `{nA0 nA1 nB0 nB1}`
   - 3D, 3 shared nodes → `-prism6` `{nA0 nA1 nA2 nB0 nB1 nB2}`
   - 3D, 4 shared nodes → `-hex8` `{nA0..nA3 nB0..nB3}`
   - otherwise no interface is generated.
4. The new element is assigned to group `eg_i` (by default), or to
   `element_geometry = eg_i` when `method_generate==-ELEMENT_GEOMETRY`,
   marked with `ELEMENT_MACRO_GENERATE = icontrol`, and `ELEMENT_DOF` /
   `ELEMENT_DOF_INITIALISED` / `NONLOCAL_ELEMENT_INFO` are initialized
   (pattern `generate_spring`).
5. `control_mesh_generate_interface_geometry`: every shared node of `iel`
   must lie inside the geometry (pattern `geometry()` in `adjust.cc`),
   otherwise the pair is skipped.
6. Idempotence: once an interface is generated for `iel`, `ELEMENT_MACRO_GENERATE`
   is set on `iel` and `jel`, so later steps skip the pair. The scan loops
   run over the element set present BEFORE the first generation
   (`max_element_old`), so the freshly generated interfaces are not
   reprocessed.

## Semántica (del manual de Professional)

- The interface element gets `element_group eg_i` if it is between
  `eg_a` and `eg_b`; `eg_i eg_a eg_b eg_j eg_c eg_d ...` generates up to
  several interfaces in one record.
- Interfaces can only be generated between exactly two elements (a common
  face shared by two blocks, e.g. duplicated nodes between two meshes).
- Crossing interfaces are not allowed.
- For connected interfaces, list the touching pairs in ONE record (the
  manual: `control_mesh_generate_interface 10 20 30 31 20 40 41`).

## Detalles

- The generated interface element uses the shared nodes of BOTH sides, so
  the two blocks must have duplicated (spatially coincident) nodes on the
  interface face.
- `control_mesh_generate_interface_geometry` filters by the geometry; the
  element pair must be inside the geometry.

## Pendiente

- Tipos cuadráticos (bar3→quad6, tet10→tria12, hex27→quad18) no
  implementados; se generan solo `quad4`/`prism6`/`hex8`.
