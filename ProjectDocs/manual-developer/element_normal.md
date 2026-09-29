# element_normal (developer)

Manual Professional 6.445; consumers 6.728 / 6.740
(`group_materi_plasti_mohr_coul_direct_normal_automatic`,
`group_materi_plasti_tension_direct_normal_automatic`). Implemented
2026-09-30.

## The defect it fixes

Before 2026-09-30 `materi()` (`materi.cc`) computed the "automatic" plane
normal as the cross product of two edges of the element (nodes 0-1 and
0-2):

- in 2D both edges lie in the xy plane, so the normal was (0, 0, +-1), a
  direction that carries no stress: the cut-off was inert. Professional test
  `tension_direct1` (column of 2 quad4, upper element with
  `tension_direct 0` + `_normal_automatic`, top pulled up by 1): GNU purely
  elastic (`disy` 0.51 at nodes 3/4 with 5 iterations, 0.42 with the default
  iterations; `sigyy` 0.48-0.67), Pro/analytic `disy` 0 and all stresses 0;
- the loop read `c0[2]`, `c1[2]`, `c2[2]` from `NODE` records that only hold
  `ndim` coordinates (out-of-bounds read in 1D/2D);
- in 3D it gave the normal of the first element face (independent of the
  mesh boundary).

## Rule (measured on the Professional 02-08-2026)

Probe inputs: 2x2, 3x3, 4x3 and single-row quad4 grids with square and
2 x 1 cells, distorted quad4, tria3 (regular and distorted), quad9/tria6
(straight and curved sides), hex8 (1x1x2, 2x2x2, 3x2x2, 3x3x3, two
distorted meshes with non-planar faces), tet4 (distorted), 1D bar2, a
column split over two groups, trusses along the boundary, an input-given
`element_normal`, `_normal_automatic -no`, `control_mesh_delete_element`,
and a 90 degree rotation with `mesh -follow_material`. The Professional
prints `element_normal` and `node_normal` in the `.dbs`, which made the
rule directly observable:

1. boundary side = side (border table) of an isoparametric element whose
   nodes are not all nodes of another isoparametric element. Groups do not
   matter (`q_col2_grp`); a `-truss` element on the boundary nodes does not
   close the side (`q_truss`);
2. node normal = normalize( sum over boundary sides of
   `integral( N_node n dA )` ), `n` oriented towards the element centre;
   zero for a node without boundary side. Plain length/area weighting is
   identical for straight sides and planar parallelogram faces; the
   distorted hex8 probes (non-planar faces) need the consistent
   shape-function weighting (max difference 8e-3 with area weighting,
   5e-13 with the consistent rule, i.e. print precision);
3. element normal = normalize( sum of the node normals of the element
   nodes ), zero when the sum vanishes. Computed for EVERY element
   (trusses get one too), not only for the elements of the automatic
   groups;
4. generated only when some group has a `_normal_automatic` record (any
   switch value: with `-no` the Professional still uses it); otherwise the
   record is absent, even if the input file gave one;
5. computed from the current node coordinates, once, and again after a mesh
   change (`control_mesh_delete_element` probe: element 5 goes from (0,0) to
   (0,1)); NOT after deformation (`mesh -follow_material` 90 degree
   rotation keeps (-0.707,-0.707); a rotation followed by an element
   deletion gives the rotated normal).

Agreement GNU vs Professional on all probes: 5e-13 (print precision),
except the curved quadratic sides (see below).

## Implementation

- `tochnog.h` / `tochnog-mod.h`: enum `ELEMENT_NORMAL` appended before
  `LAST_DUMMY` (both headers in sync); prototypes `element_normal_set()`,
  `element_normal_needed()` next to `area_node_dataitem()`.
- `database.cc` (after `ELEMENT_MIDDLE`): `element_normal`,
  DOUBLE_PRECISION, `data_length` = ndim, `external` 1 (printed in the
  database and accepted in an input/restart file), `version_all` 0
  (VERSION_NORMAL only), class ELEMENT, required ELEMENT.
- `area.cc` (end of file, reuses the static `border_nodes_*` side tables of
  `area()`):
  - `element_normal_sides()`: side table per element type (bar2 in 1D;
    tria3, tria6, quad4, quad9, quad16 in 2D; tet4, hex8, hex27 in 3D) and
    the side kind (end point, Lagrange line, linear triangle, k x k
    lexicographic face);
  - `element_normal_lagrange()`: 1D Lagrange shape functions (equally
    spaced nodes on [-1,1]) and derivatives;
  - `element_normal_side_weights()`: consistent nodal side normals, 4-point
    Gauss per direction (full-precision constants in the file; exact for
    all supported sides);
  - `element_normal_needed()`: any `GROUP_MATERI_PLASTI_TENSION_DIRECT_
    NORMAL_AUTOMATIC` / `..._MOHR_COUL_DIRECT_NORMAL_AUTOMATIC` record;
  - `element_normal_set()`: returns unless needed and absent; builds a
    node -> side-bearing element list (CSR), tests every side of every
    side-bearing element against the elements of its first node,
    accumulates the oriented consistent weights per node, normalizes
    (relative threshold 1e-12 of the summed contribution sizes), then PUTs
    the normalized sum per element (array_normalize threshold 1e-10, zero
    otherwise; -0 printed as 0).
- `top.cc` `iteration_start()`: calls `element_normal_set()`. This runs
  before every `element_loop()` and outside the parallel element loop
  (allocation there aborts with "Data is allocated in a parallel loop").
- `mesh.cc` `mesh_has_changed()`: `db_delete( ELEMENT_NORMAL,
  VERSION_NORMAL )` (only for VERSION_NORMAL: the record has no other
  version). This is what regenerates the record after mesh changes and
  discards input values at the start (top() calls mesh_has_changed).
- `materi.cc` `materi()`: the `_normal_automatic` branch reads
  `ELEMENT_NORMAL` (GET_IF_EXISTS) into `direct_normal`; the cross product
  and the out-of-bounds reads are gone. A zero normal leaves
  `direct_normal` zero, so `set_stress()` takes the full principal-stress
  cut-off (`materi_direct_full_mc`).

## Hardcoded parameters

- 4-point Gauss rule (exact up to degree 7: quad16 lines need degree 5,
  hex27 faces degree 5 per direction).
- Node-normal threshold: |sum| <= 1e-12 * sum|contributions| -> zero.
- Element-normal threshold: `array_normalize` (|sum| < 1e-10 -> zero).

## Pending / differences

- Zero element normal with `_normal_automatic`: the Professional aborts
  with an assertion, the GNU falls back to the full cut-off (kept on
  purpose, documented in the user manual).
- Curved quad9/tria6 sides: the Professional does not follow the
  consistent rule (probe `q9d_2x2`/`t6d_2x2`: some node normals come out
  axis-aligned, e.g. (0,-1) on a side whose normal is (-0.6,0.8)); not
  reproduced.
- Element types without a border table (tet10, hex20, hex64, prism6, quad8,
  ...) contribute no boundary side.
- `node_normal`, `node_boundary`, `element_boundary` (printed by the
  Professional) and `global_normal` are not implemented.
- A restart from a `.dbs` recomputes the normal from the current node
  coordinates (mesh_has_changed at the start of top()); the Professional
  restart behaviour was not measured.
- `mohr_coul_direct` plane mode with the automatic normal: the normal now
  matches the Professional, but a sheared-column probe still differs in
  the state (backlog C14 item 4, plane return `materi_direct_cutoff`).
- Nodal strain recovery differs in `tension_direct1` (`eptyy` 0.667 vs
  0.5 at the middle nodes): audit cluster K11, unrelated to this record.
