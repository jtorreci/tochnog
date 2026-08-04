# input_gmsh

## Archivos y funciones

- `input.cc` → `input_gmsh_read()` (`input.cc:1474`) — reads the whole Gmsh
  file `tochnog_in.msh` (format 2.2) into the database. Only nodes, elements
  and element groups are read.
- `input.cc:599-618` — the main data-part loop branch `idat==INPUT_GMSH`
  parses the `-yes`/`-no` switch and calls `input_gmsh_read()` when the switch
  is `-YES`.
- `database.cc:2890` — keyword registration:
  `strcpy(name[INPUT_GMSH],"input_gmsh")`.
- Enum `INPUT_GMSH` in `tochnog.h` (`tochnog.h:634`) / `tochnog-mod.h` (must
  stay in sync).

## Detalles de implementación

- The `MeshFormat` section is skipped; the reader scans forward with
  `while ( in >> str )` until it finds `$Nodes`, then `$Elements`
  (`input.cc:1494-1507`).
- `$Nodes` is always read with 3 coordinates, `x y z`:
  `in >> node_tag >> xyz[0] >> xyz[1] >> xyz[2];` and stored with
  `db( NODE, node_tag, idum, xyz, ndim, VERSION_NORMAL, PUT )`
  (`input.cc:1499-1501`). `ndim` limits how many coordinates are used.
- `$Elements` lines have the form
  `elmnr elm-type number-of-tags <tags> node-list`
  (`in >> igeom >> ieltype >> ntag`). The `ntag` tag values are skipped
  (`input.cc:1510`).
- Element type mapping (`input.cc:1512-1517`):
  - 1 → `-BAR2`, 8 → `-BAR3`, 2 → `-TRIA3`, 9 → `-TRIA6`,
    3 → `-QUAD4`, 10 → `-QUAD9`.
  - Any other type aborts with
    "Error: gmsh element type X not supported by tochnog."
- Elements are stored with
  `db( ELEMENT, igeom, elem_data, ddum, ntag=1+nnol, VERSION_NORMAL, PUT )`
  where `elem_data[0]` is the negated element type and `elem_data[1..nnol]`
  the node list (`input.cc:1522-1524`).

## Dependencias externas

None beyond the C++ standard library (`<fstream>`). The file format
assumed is Gmsh 2.2 ASCII; no Gmsh library is linked.

## Parámetros hardcodeados / refactorizaciones pendientes

- The input file name `tochnog_in.msh` is hardcoded in `input.cc:1486`.
- The reader always expects exactly 3 coordinates per node, even for 2D
  meshes; `ndim` is set by the control data before the read.
- The type mapping table is written as a chain of `else if`; it could be
  extracted into a `struct { int gmsh_type; int tn_type; int nnol; }` table
  to make the supported types explicit and easier to extend.
- 3D elements are intentionally unsupported (no Tochnog equivalent); the
  error path could list the recognized 3D types to give a better message.
- Physical groups (first tag) are read but currently not used to build
  groups; only the comment "element group from gmsh physical group" remains.
