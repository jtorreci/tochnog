# control_mesh_move

## Files and functions

- `mesh.cc` — `void mesh_move( double control_mesh_move[], long int length )` (line 183). Implementation of the linear mesh displacement.
- `top.cc` — called in `step_start` (lines 501–507) when `db_active_index( CONTROL_MESH_MOVE, icontrol, VERSION_NORMAL )` is true.
- `tochnog.h` — enum `CONTROL_MESH_MOVE` (line 241) and prototype (line 1356).
- `database.cc` — keyword registration (lines 662–666).

## Implementation details

- The coefficients are organized as `coeff[idim][0] = constant` and
  `coeff[idim][1+jd] = linear coefficient of axis jd`, with `idim, jd < ndim`.
  Only blocks that fit in `length` are read
  (`idim < ndim && (idim*(ndim+1)+ndim) < length`).
- For every active node (`db_max_index(NODE, ...)` + `db_active_index`) the new
  coordinates are
  `new_coords[idim] = coords[idim] + coeff[idim][0] + sum_jd coeff[idim][1+jd]*coords[jd]`.
- The same transform is repeated on `NODE_START_REFINED` when the node is active there.
- Closes with `mesh_has_changed( VERSION_NORMAL )`.
- `database.cc`: `type = DOUBLE_PRECISION`, `data_length = DATA_ITEM_SIZE`,
  `fixed_length = 0`, `data_class = CONTROL`.

## External dependencies

None. Uses only the internal database API (`db`, `db_active_index`,
`db_max_index`) and `mesh_has_changed`.

## Hardcoded parameters / pending refactorings

- The coefficient block layout `idim*(ndim+1)+1+jd` is implicit in the keyword
  syntax; there is no validation that the given block count matches `ndim`.
- The transform loop is duplicated for `NODE` and `NODE_START_REFINED` (the same
  duplication appears in `mesh_mirror`/`mesh_copy`); a shared node-transform
  helper could be extracted.
