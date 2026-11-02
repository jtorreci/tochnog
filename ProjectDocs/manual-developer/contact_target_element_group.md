# contact_target_element_group: generated node_boundary (developer)

## Where

- `area.cc`, `void mesh_boundary_set( void )` (renamed from
  `contact_node_boundary_set` in F5 T2 V3, see `mesh_boundary.md`; after
  `element_normal_set`): generates `NODE_BOUNDARY` (`-yes`) for the nodes
  of the element sides that no other side-bearing element shares. Runs
  only when `CONTACT_TARGET_ELEMENT_GROUP 0` is active and no
  `NODE_BOUNDARY` record exists.
- `top.cc`, `iteration_start()`: calls it after `element_normal_set()`,
  i.e. serially before `element_loop()` and `parallel_contact()` (db
  allocation is not allowed in the parallel loops).
- Prototype in `tochnog.h` / `tochnog-mod.h`.
- Consumer: `contact.cc` `parallel_contact()`: element targets are looped
  only when `db_max_index(NODE_BOUNDARY) >= 0`, and a target side is used
  when `number_of_side_boundary_nodes==ndim`.

## Implementation details

- Side tables: the static `element_normal_sides()` of area.cc (1D bar2;
  2D tria3/tria6/quad4/quad9/quad16; 3D tet4/hex8/hex27). Node ->
  element compressed lists and the "another element holds all side
  nodes" test are the same as in `element_normal_set()`.
- Elements of other types (truss, beam, spring, interface, ...) are
  ignored: they neither define nor hide boundary sides.
- `mesh_has_changed()` deletes NODE-class records (node_boundary
  included); the next `iteration_start()` regenerates the boundary for
  the new mesh.
- Hand-written `node_boundary` (GNU 2014 inputs) disables the generation,
  so those inputs keep their meaning (sfnet `contac1`, `contac4`).

## Verification

- Professional 02-08-2026 node_boundary sets are equal: `contac1` nodes
  1, 2, 3, 5 (interior node 4 excluded); 2D 3 x 2 quad4 block, all nodes
  except the two interior ones.
- `contac1` velx 0.02/0.01 = Pro (before 0.03/0); `contac4` temperatures
  1.5 and dofs = Pro exactly; `axisym3` velx -9.99981e-4 (Pro -9.99981e-4).
- Tests `contact_element_target_boundary1` (1D) and
  `contact_element_target_boundary2d1` (2D quad4): RED rc 1 old / GREEN
  rc 0 new; the Professional meets both targets.

## Pending (BACKLOG C40)

- No check "target element has all nodes on the boundary" (Pro error).
- `number_of_side_boundary_nodes==ndim` rejects hex8 sides (4 nodes).
- No `element_boundary` record; node_boundary not generated for other
  analyses (the Pro prints both in every run).

## Differences with the Professional version

- **element-target contact residuals** (known difference, not decided;
  BACKLOG C40; AUDIT K05). GNU: Runs a single-layer target (non-physical
  result); hex8 target sides never detected; `node_boundary` only for
  element-target contact, no `element_boundary`. Professional 02-08-2026:
  Stops ("a target element has all nodes on the boundary"); writes
  `node_boundary`/`element_boundary` in every run. Status: Open, LOW. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
