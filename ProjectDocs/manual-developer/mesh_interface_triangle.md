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

## State of the physics (interface11, CLOSED 2026-09-07)

The cut routine generates the structure that matches the Professional
.dbs of interface11 (element types, numbering, groups, node duplication
— verified element by element) AND the physics now converges:
`element_interface_stress_average 4 0` = −1.00000008 (target −1.0,
tol 1e-2, rc=0; the Professional solves −0.99999954).

The blocker was a PRE-EXISTING bug of the GNU volume integration of the
`-prism6` wedge in polynom.cc (independent of this cut routine). It was
never caught because the corpus tests that involve triangular interfaces
only targeted bulk stresses with a loose tolerance or prescribed fields
(the stress reads C·B·u and never exercises the volume integral).
Root cause and fix are documented in polynom.cc (PRISM6 branch and the
volume[] branch): the wedge fell through to the hex8 integration
(weight*8*detj) with weights summing to 1.5, i.e. 24x the physical
reference volume 0.5 — every volume integral of the wedge (stiffness,
mass, gravity, nodal face forces) was 24x too large, so the nodal force
a compressed wedge exerts on its triangular faces was 24x the
consistent load (measured: u=−z, E=1, σ=−1 → reactions ±4 per node
instead of ±A/3 = ±1/6). Zero-thickness prism6 interfaces against
`-prism6`/`-tet4` volume neighbours therefore converged to σ_iface =
−24 instead of −1 (minimal model), while hex8-neighbour interfaces were
exact. The rewrite uses the degree-2 triangle rule × 2-point Gauss in
zeta (weights summing to the reference volume 1/2), the layout the
Professional integrates (its element_intpnt_coord of the wedge
elements shows zeta = 0.2113/0.7887 = 1/2 ± 1/(2√3)), which also fixes
the OBLIQUE wedges of the interface11 cut (their Jacobian varies along
zeta; the old rule with both z-levels in the upper half bent the field
away from u=−z: node_dof of the cut plane −0.626..−0.758 instead of
−0.6 and per-pair σ_iface −0.32..−0.88 instead of uniform −1).

Minimal reproductions (GNU post-fix vs Professional 25-10-2023):

- Two stacked `-prism6` wedges (triangular cross-section, height 0.5
  each) + one zero-thickness `-prism6` interface, top displaced −1,
  bottom fixed: GNU σ_iface = −1.00000008, reactions ±0.1666664/
  +0.1666669 per node (the consistent ±1/6); the Professional solves
  −1.0 with reactions ±1/6. Before the fix the GNU solved σ = −24 with
  reactions ±4 (24x), invariant to the side sign.
- `interface11` of the corpus (3 tets cut by the plane z=0.6): the
  three generated interfaces (elements 4/5/6, the target element 4 is
  the prism6 one) give σ = −1.00000008 / −1.00000286 / −1.00000748
  (Pro −0.99999954 / −0.99999830 / −1.00000026) and the cut-plane
  node_dof is −0.6000014 (u=−z within the solver tolerance; before the
  fix −0.626..−0.758 and σ_el4 = −13.74). rc=0.
- `interface_tria3_prism6` of the corpus (prism6 interfaces converted
  from tria3 on the shared quad face of two hex8 blocks) is UNCHANGED
  by this fix (its interface elements do not go through the volume
  integration of pol()); rc=0 kept. NOTE its own latent bug is still
  open: the GNU interface stress is +1.51 where the Professional gives
  −1.0 — the triangulated quad face (two prism6 interfaces sharing the
  diagonal) distributes the contact force 1/6-1/3-1/3-1/6 over the
  four face nodes while the consistent quad load is 1/4 each, and the
  orientation of the converted zero-thickness interface reports
  compression with the opposite sign. It passes rc=0 only because its
  target checks a bulk post_point stress. Fixing it requires touching
  the tria3→prism6 interface conversion / orientation (out of scope of
  the volume-integration fix of this sprint).
