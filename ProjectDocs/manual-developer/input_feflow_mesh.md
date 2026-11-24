# input_feflow_mesh

## Archivos y funciones

- `input.cc` → `input_feflow_read()` (`input.cc:1726`) — reads the FEFLOW
  mesh and fills `node` and `element` records.
- `input.cc:578-597` — the main data-part loop branch `idat==INPUT_FEFLOW_MESH`
  parses the `-yes/-no` switch and calls `input_feflow_read()` when the
  switch is `-YES`.
- `database.cc:2872-2888` — keyword registrations: `input_feflow_fem`,
  `input_feflow_mesh`, `input_feflow_mesh_hydraulic_head`.
- Enums `INPUT_FEFLOW_FEM`, `INPUT_FEFLOW_MESH`,
  `INPUT_FEFLOW_MESH_HYDRAULIC_HEAD` in `tochnog.h`
  (`tochnog.h:631-633`) / `tochnog-mod.h` (must stay in sync).

## Detalles de implementación

- The input file is chosen by `input_feflow_fem` (read with GET_IF_EXISTS,
  default `-YES`): `-NO` opens `feflow.dac`, otherwise `feflow.fem`
  (`input.cc:1739-1744`).
- The file is parsed line by line with `in.getline(line,MCHAR)` and
  lowercased (`input.cc:1752-1753`).
- Sections are detected by substring keywords on each line:
  `coordinates` or `nodes` switches to node reading, `elements` switches to
  element reading (`input.cc:1759-1764`). Other header lines
  (`problem`, `dimension`, `version`, `feflow`, `#`, `*`, `/`) are skipped
  (`input.cc:1766-1767`).
- Elements: first token is the element number, then the node numbers, split
  by `strtok` on `" ,\t"`. The element type is derived from the node count
  (`input.cc:1790-1795`): 2→`-bar2`, 3→`-tria3`, 4→`-quad4`, 6→`-tria6`,
  9→`-quad9`; any other count is skipped. Stored with
  `db( ELEMENT, ielem, edata, ddum, nn=1+nnol, VERSION_NORMAL, PUT )`.
- Nodes: first token is the node number, then up to 3 coordinate tokens,
  stored with `db( NODE, inode, idum, xyz, ndim, VERSION_NORMAL, PUT )`
  (`input.cc:1803-1812`).

## Dependencias externas

None beyond the C++ standard library (`<fstream>`). FEFLOW file parsing is
done manually; no FEFLOW library is linked.

## Parámetros hardcodeados / refactorizaciones pendientes

- Input files `feflow.fem` / `feflow.dac` are hardcoded (`input.cc:1742-1744`).
- `input_feflow_mesh_hydraulic_head` is registered but not implemented; the
  keyword is never acted upon.
- Section detection is substring-based (`strstr`), so a line containing the
  word "elements" anywhere (e.g. a comment) switches the parser state; a
  stricter token/heading match would be more robust.
- The `.dac` file format: real FEFLOW `.dac` files store result data, not a
  mesh. `input_feflow_fem -no` therefore assumes a `.dac` that contains
  `coordinates`/`elements` sections — effectively a renamed mesh file.
- The element-by-node-count dispatch could be a small lookup table instead
  of the `if/else` chain.
