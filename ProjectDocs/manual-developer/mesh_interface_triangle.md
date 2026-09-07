# mesh_interface_triangle family (coordinate / element_group / control)

## Keywords (database.cc, appended at the end of the enum)

- `MESH_INTERFACE_TRIANGLE_COORDINATE` — DOUBLE_PRECISION, variable
  length (9 doubles per triangle), data_class CONTROL. Manual
  Professional 6.856. SINGULAR spelling `mesh_interface_triangle_
  coordinate` (the manual text says "coordinates"; the Professional
  binary 25-10-2023 and the corpus `interface11` use the singular form —
  verified against the Pro .dbs echo).
- `MESH_INTERFACE_TRIANGLE_ELEMENT_GROUP` — INTEGER, length 1, data_class
  CONTROL (6.857): the group of the generated interface elements.
- `CONTROL_MESH_INTERFACE_TRIANGLE` — INTEGER, length 1, data_class
  CONTROL (6.201): `-yes` activates the generation at the control index.
- `CONTROL_PRINT_GID_METHOD` — INTEGER, length 1, data_class CONTROL,
  registered PARSE-ONLY. The corpus `interface11` carries the record
  (GiD output method selector) and the GNU GiD printer is
  self-contained; without the registration the .dat cannot parse (the
  unknown token would abort the parse of the preceding variable-length
  record). No method dispatch.

## Dispatch

`generate_interface_triangle( icontrol )` in generate.cc, invoked from
`step_start()` in top.cc right before `generate_interface()`, BEFORE the
any_interface scan (so the interface histories of the generated elements
are allocated at the first task==YES step). The routine returns
immediately unless `control_mesh_interface_triangle` is active at the
current icontrol with value `-yes`.

## Generation algorithm (generate.cc)

Per control index, over the element set present BEFORE the call (pattern
`generate_interface`: generated elements are never rescanned):

1. For each triangle of the coordinate record: classify every `-tet4`
   vertex by the signed distance to the plane. A tet is cut when it has
   vertices on both sides and NO vertex on the plane (|c| < 1e-9 skips
   the tet), and when the centroid of its cut polygon lies inside the
   current triangle (barycentric test in the plane; boundary counts as
   inside). A tet is cut at most once (first matching triangle wins).
2. Cut points are computed per crossed edge and deduplicated by
   coordinates (1e-10): the point of an edge shared by two tets is
   created once, with two node copies:
   - "base" copy: created first, used by the pieces on the +normal side
     of the plane (n = (P1-P0)x(P2-P0) of the cutting triangle);
   - "duplicate" copy: created on demand when the interface of the first
     tet that needs it is generated, used by the -normal side pieces.
   Node numbering mirrors the Professional .dbs of interface11: base
   nodes 7..12 (edge scan of the crossed tets, canonical tet4 edge
   order), duplicates 13..18 (in the order the interfaces reference
   them).
3. Element passes (Professional .dbs numbering of interface11:
   interfaces 4,5,6, tet4 pieces 7,8, prism6 pieces 9..12):
   - PASS 1: one interface element per cut tet, in tet order. The cut
     polygon is ordered counter-clockwise seen from +n (angle sort
     around the centroid). The element is `-prism6` (3 points) or
     `-hex8` (4 points) with side 1 = the -normal side copies and
     side 2 = the +normal side copies, coincident per slot.
   - PASS 2: the `-tet4` halves of the 1+3 splits (the piece at the
     isolated vertex; face ordered so its e1xe2 points towards the
     apex: CCW+ base copies when the isolated vertex is on the + side,
     reversed duplicate copies when on the - side).
   - PASS 3: the `-prism6` halves (per tet: -side piece first, then
     +side piece). Each wedge connects the original tet face (or a
     subtriangle of it) to the cut section slot-to-slot; a straight
     frustum is never twisted, but the reference can be mirrored
     depending on the input tet node order: `cut_fix_prism_orientation`
     mirrors BOTH end triangles when the reference volume
     ((b-a)x(c-a)).(d-a) is negative (the physical element is
     unchanged).
   - The original tets are deleted (`delete_element`).
4. New nodes are created with `cut_create_node`: NODE-class records are
   copied from the first edge endpoint and NODE / NODE_START_REFINED are
   overwritten with the cut coordinate; NODE_DOF / NODE_DOF_START_
   REFINED are INTERPOLATED between the two edge endpoints at the edge
   parameter (a mid-analysis cut inherits the interpolated kinematic
   state). NODE_ELEMENT / NODE_NODE / NODE_GEOMETRY_PRESENT are not
   copied (rebuilt / recomputed downstream).
5. `mesh_has_changed` at the end (NODE_ELEMENT/NODE_NODE invalidated and
   rebuilt lazily).

## Scope / limitations of this implementation

- Linear `-tet4` volumes only.
- A tet with a vertex exactly on the plane is not cut.
- A cut polygon crossing a triangle boundary of the record is not cut
  (the whole polygon must sit inside ONE triangle of the record; the
  Professional cuts with the union of the triangles).
- Zero-thickness interfaces: both copies share the exact coordinates, so
  the GNU interface frame is e1xe2 of the side-1 face without the
  orientation flip (dir = n.(cm1-cm2) = 0 exactly). The convention
  (side 1 = -normal side, CCW+ ordering, side 2 = +normal side) makes
  compression give a NEGATIVE normal stress, like the Professional.

## State of the physics (interface11, NOT closed)

The cut routine generates the structure that matches the Professional
.dbs of interface11 (element types, numbering, groups, node duplication
— verified element by element), but the target
`element_interface_stress_average 4 0 = -1` is NOT reached (GNU gives
about -13.7 with kn=1e11). Root cause measured: the GNU interface
element is NOT consistent with `-prism6`/`-tet4` volume neighbours. The
bug is PRE-EXISTING (independent of this cut routine) and was never
caught because the corpus tests that involve triangular interfaces
(`interface_tria3_prism6`) only target a BULK stress with a loose
tolerance. Minimal reproductions (GNU vs Professional 25-10-2023):

- Two stacked `-prism6` wedges (triangular cross-section, height 0.5
  each) + one zero-thickness `-prism6` interface, top displaced -1,
  bottom fixed: exact answer sigma = -1 (the Professional solves -1.0);
  the GNU solves sigma = -24, i.e. an interface jump 24x the
  equilibrium one, for every kn (asymptotic). With `-hex8` blocks the
  same setup is exact (-1) in the GNU.
- `interface_tria3_prism6` of the corpus (prism6 interfaces converted
  from tria3 on the shared quad face of two hex8 blocks): the GNU
  interface stress is +1.51 where the Professional gives -1.0 (the
  corpus test passes because its target is a bulk post_point stress).

The interface element machinery (interface.cc) was calibrated on
quadrilateral (hex8-face) interfaces (uniform per-pair weights 1/ns1,
face-area measure). With triangular sides the assembled equilibrium of
the interface + wedge/tet pieces is off by a geometry-dependent factor
(24 in the minimal model). Closing interface11 needs the interface
element to be made consistent with `-prism6`/`-tet4` neighbours (and a
re-validation of the whole interface family: interface1-15, conspr1-7,
patch1, tria3_prism6, quad4_hex8) — that is the natural continuation
sprint; the present cut routine is ready and verified structurally.
