# hex20

## Design rationale: Lagrange upgrade over plain serendipity

`-hex20` is accepted as INPUT (the conventional serendipity 20-node
format of Abaqus/GiD/most mesh generators) but is never assembled as a
serendipity element: it is upgraded internally to the complete 27-node
Lagrange `-hex27` (see [quad8](quad8.md) for the full rationale). The
Professional follows the same path: its own suite states textually that
hex20 elements "will be automatically converted to ... hex27 volume
elements" (interface_quad8_hex20.dat).

## Implementación

- **Keyword**: `hex20` registered in `database.cc` (name table only),
  enum `HEX20` in `tochnog.h`/`tochnog-mod.h` between `HEX18` and
  `HEX27` (kept in sync). No element routine handles `-hex20`: the
  element is converted before any assembly.
- **Conversión**: `mesh_convert_hex20()` in `mesh.cc`, called from
  `step_start()` in `top.cc` on EVERY step start, right after
  `mesh_convert_quad8()` and BEFORE `extrude()` /
  `interface_convert()` (the 3D `-quad8` interface split must see the
  hex27 bulk with the shared face centre node already in place). Runs
  on `task==NO` intermediate control steps too (a post_point below the
  timestep would otherwise hit the raw hex20). Idempotent: converted
  elements are `-hex27` and skipped.
- **Permutación de nodos** (the delicate point; NOT guessed): verified
  against the Professional 25-10-2023 .dbs. When the Pro converts the
  mesh of its corpus `interface_quad8_hex20` it writes, for the input
  `element 1 -hex20 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20`
  (corners: base BL,BR,TL,TR = 1,2,3,4 / top = 5,6,7,8; base mid-edges
  BM,LM,RM,TM = 9,10,11,12; top mid-edges = 13,14,15,16; verticals at
  (x0,y0),(x1,y0),(x0,y1),(x1,y1) = 17,18,19,20):
  ```
  element 1 -hex27 1 9 2 10 56 11 3 12 4 17 57 18 60 61 59 19 58 20
                    5 13 6 14 55 15 7 16 8
  ```
  which is EXACTLY the GNU tensor hex27 layout (mesh_extrude /
  border_nodes_hex27 conventions: base plane quad9 tensor
  `[BL,BM,BR,LM,C,RM,TL,TM,TR]`, mid plane, top plane), i.e.:
  - base plane slots 1..9: corners el[1,2,3,4] at 1,3,7,9; mid-edges
    el[9,10,11,12] at 2,4,6,8; slot 5 = new base face centre.
  - mid plane slots 10..18: verticals el[17,18,19,20] at 10,12,16,18;
    slots 11,13,15,17 = new y0/x0/x1/y1 face centres; slot 14 = body.
  - top plane slots 19..27: corners el[5,6,7,8] at 19,21,25,27;
    mid-edges el[13,14,15,16] at 20,22,24,26; slot 23 = new top face
    centre.
- **New nodes**: the 6 face centres = average of the 4 corners of the
  face; the body centre = average of the 8 corners. DEDUPLICATED BY
  COORDINATES (EPS 1e-10) against every active node of the mesh: two
  stacked hex20 sharing a face must end with ONE centre node on the
  shared face (hex20.dat of the corpus stacks two hex20 on z=0..2 and
  the shared z=1 face centre appears once in both hex27 records),
  otherwise the mesh tears apart. New nodes copy the NODE-class records
  of the first corner (mesh_convert_quad8 pattern).
- **Tests**: corpus `hex20.dat` (2 stacked hex20, uniaxial traction,
  post_point sigzz = 1.0) rc=0. The auto-converted run is bit-equal
  (1.5e-15 over the whole node_dof) to the same mesh written as explicit
  `-hex27` records, and GNU node_dof matches the Professional .dbs
  (worst delta 3.8e-10 on displacements, 2.1e-8 on sigzz over the 32
  common nodes; the 13 interior nodes of the GNU have no Pro counterpart
  because the Pro keeps hex20 natively).
- **NOTA Pro**: the Pro binary runs `hex20.dat` standalone rc=0 and
  KEEPS `-hex20` in its .dbs (native serendipity kernel there); the
  auto-conversion to hex27 is what its own suite documents textually for
  the interface files. The GNU always converts.

## Pendiente

- None for the volume conversion. The hex20 nodes are searched once per
  conversion step (idempotent); a large hex20 mesh pays a linear scan
  per new centre node on the first step only.
