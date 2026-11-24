# quad8

## Design rationale: Lagrange upgrade over plain serendipity

`-quad8`/`-hex20` are accepted as INPUT (conventional serendipity
formats, standard in Abaqus/GiD/most mesh generators) but are never
assembled as serendipity elements. They are upgraded internally to the
complete Lagrange `-quad9`/`-hex27`, which is a deliberate improvement
over a plain 8/20-node serendipity implementation:

- Serendipity elements drop the interior nodes, so their polynomial is
  incomplete (they miss the `xi^2*eta^2`-type bubble term of the full
  quadratic); the Lagrange upgrade adds those interior degrees of
  freedom and restores the complete quadratic field on the same input
  mesh.
- This mirrors Tochnog Professional: the corpus suite itself states that
  quad8/hex20 "will be automatically converted to quad9 ... / hex27
  volume elements". GNU node_dof agrees with the Professional `.dbs` to
  ~1e-7/1e-16 on the patch tests, which confirms the Professional
  follows the same conversion path rather than a serendipity kernel.
- The alternative (implementing real serendipity shape functions) would
  touch ~20 files (pol(), point_el, area border tables, exporters,
  split/refine) for a strictly weaker element, with no benefit over the
  reference binary.

## Implementación

- **Keyword**: `quad8` registered in `database.cc` (name table only,
  like the other element names; `strcpy(name[QUAD8],"quad8")`), enum
  `QUAD8` in `tochnog.h`/`tochnog-mod.h` between `QUAD6` and `QUAD9`
  (kept in sync). No element routine handles `-quad8`: the element is
  converted before any assembly.
- **Conversión**: `mesh_convert_quad8()` in `mesh.cc`, called from
  `step_start()` in `top.cc` on EVERY step start (not only `task==YES`:
  a `-quad8` is not a native element, so any evaluation between control
  steps - e.g. a `post_point` at an intermediate control index below the
  timestep - would hit the raw quad8; `interface_bar3_quad8` of the
  corpus exposed it). It runs BEFORE `extrude()` (a quad8 mesh must
  become quad9 first so the quadratic extrusion lifts it to hex27) and
  before `interface_convert()` (the interface split must see the quad9
  bulk). Idempotent: converted elements are `-quad9` and skipped.
- **Permutación de nodos** (el punto delicado):
  - quad8 (Professional, every `-quad8` record of the suite):
    `[BL, BR, TL, TR, BM, LM, RM, TM]`.
  - quad9 (GNU, tensor ordering xi fastest / eta slowest, -1 -> 0 -> +1):
    `[BL, BM, BR, LM, CENTRE, RM, TL, TM, TR]` (corners 0,2,8,6 and
    mid-edge nodes 1,3,5,7 in `border_nodes_quad9` of `area.cc`; same
    ordering as the quad9 records of the suite, e.g. `patch1.dat`).
  - Mapping: `quad9 = { q8[1], q8[5], q8[2], q8[6], new_centre, q8[7],
    q8[3], q8[8], q8[4] }`.
  - Verificado contra el binario del Professional (25-10-2023, .dbs):
    `elasti6` -> element `-quad9 1 5 2 6 9 7 3 8 4` (elasti6 is the
    quad8 `{1..8}` with corners 1,2,3,4) and the same slot permutation
    in `interface_bar3_quad8` (Pro: quad9 `{1,2,3,4,C,5,TL',TM',TR'}` =
    same mapping with the interface edge duplicated).
- **Nodo central**: average of the 4 corners; the new node copies the
  NODE-class records of the first corner (NODE_DOF etc., all zero at
  first step) like `interface_convert`/`mesh_extrude`. Node numbering
  continues from the highest existing node.
- **Skip**: `-quad8` elements of an INTERFACE group are left alone (the
  3D `-quad8` interface family is a facial element: it is lifted to the
  `-hex18` quadratic 3D interface by `interface_convert()`, see
  [hex18](hex18.md)).
- **Tests**: corpus `quad8.dat` (2 quad8, patch test sigxx = 0.1),
  `elasti6.dat` (1 quad8 vertical, force_edge top, targets sigyy(node 8)
  = 1.0 ± 1e-4 and post_force_edge_summed = 10), `interface_bar3_quad8`
  (quad8 + bar3 interface, sigyy = -1.0). GNU node_dof vs Pro .dbs:
  agreement ~1e-7 (sigyy node 8 GNU 1.00000008274 vs Pro
  0.9999998606958).

## Pendiente

- None: `-hex20` (3D serendipity volume, the same Lagrange upgrade
  auto-converted to hex27 by `mesh_convert_hex20`) and the 3D `-quad8`
  interface -> hex18 conversion are implemented (see [hex20](hex20.md)
  and [hex18](hex18.md)).
