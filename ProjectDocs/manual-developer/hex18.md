# hex18

## Implementación

The `-hex18` quadratic 3D interface (2 quad9 sides, ns1 = 9 facing
pairs) is supported in two places:

### interface_element() — interface.cc

- **Element type**: `name==-HEX18` -> `nnol = 18` (side 1 = first 9
  nodes, side 2 = last 9, quad9 tensor order on both).
- **Frame**: for a hex18 the side-1 nodes are a full quad9 face, so the
  surface corners are the side-1 nodes 0, 2 and 6 (BL, BR, TL); nodes 1
  (BM) is collinear with BL/BR and the naive `nodes[0],nodes[1],nodes[2]`
  cross product would degenerate. The 3D frame code picks the corner
  indices `ic1=2, ic2=6` when `name==-HEX18` (1,2 otherwise).
- **Integration weights**: per facing pair (i, i+ns1), the 2D product of
  the 1D Gauss-Lobatto rule over the quad9 face: w_ip =
  {1,4,1,4,16,4,1,4,1}/36 in tensor slot order (BL..TR). The linear 3D
  interfaces keep the uniform 1/ns1 (verified: interface_quad4_hex8);
  the quadratic face MUST use the Lobatto product - interface3 of the
  corpus loads the top side with -1/-4/-16 (the Lobatto pattern) and
  every integration point reports the same uniform stress -9. With
  uniform weights the corner intpnt would read -20.25 instead.
- **Measure**: the face area is computed from the quad9 face corners
  (side-1 nodes 0, 2, 6, 8); the measure block now reads its 4 corners
  through an index list (hex18: {0,2,6,8}; other elements: the previous
  {0,1,2,i3} layout). The tensor corner order is a bowtie polygon, so the
  triangle-fan area is used as before.
- Everything else (per-pair spring assembly, gap/tension/Mohr-Coulomb
  state, history records, intpnt output records, print_interface_stress)
  is generic over ns1 and needs no hex18-specific branch.

### interface_convert() — quad8 face -> hex18

`control_mesh_convert` now accepts the 3D facial `-quad8` of an
interface group (`interface_quad8_hex20` family):

1. The quad8 record (corners BL,BR,TL,TR + mid-edge BM,LM,RM,TM, no
   centre) is read into a local copy.
2. The 9th side-1 node is the face centre = average of the 4 corners.
   It is looked up BY COORDINATES over the active nodes: when the
   neighbouring volumes are hex20, the mesh_convert_hex20
   auto-conversion (which runs earlier in the same step_start) already
   created the deduplicated shared-face centre. If not found (no hex20
   volumes), the centre node is created with the NODE-class records of
   the first corner.
3. Side 1 is rewritten to the GNU quad9 tensor order with the SAME slot
   permutation as mesh_convert_quad8: tensor = {q8[1], q8[5], q8[2],
   q8[6], centre, q8[7], q8[3], q8[8], q8[4]}.
4. The 9 side-2 nodes are new copies of side 1 shifted `0.01` along the
   interface normal (normal = cross of the quad8 corners BL,BR,TL -
   the frame machinery of the linear conversions), with NODE /
   NODE_START_REFINED / NODE_DOF / NODE_DOF_START_REFINED /
   NODE_MACRO_GENERATE copied.
5. The element is rewritten `-hex18` [side1 tensor][side2 dups] and the
   neighbouring solids on the OTHER side (all 9 side-1 nodes shared +
   centroid on the +normal side) are reconnected to the duplicates -
   the same reconnection rule as the linear conversions.

The converted mesh matches the Professional .dbs layout exactly: Pro
`element 3 -hex18 5 13 6 14 55 15 7 16 8 68..76` for the corpus file
(55 = the shared face centre, 68..76 = the reconnected side-2 block).

### Record lengths (top.cc)

The interface history/output records are flat arrays of
`db_data_length(idat)` per element, and the db() PUT rejects a length
above it. The database defaults fit ns1 <= 4 (bar2/quad6/hex8); a hex18
(ns1 = 9, intpnt records 9*3 = 27 values) dies with "Length too small of
element_interface_strain_normal". The `any_interface` allocation block
in step_start() now RAISES (never lowers - shrinking the flat record
stride corrupts every read/write, exposed by interface2/10/patch of the
corpus) the lengths of the 9 interface records to the largest interface
element of the model before the db_allocate calls.

## Tests

- corpus `interface3.dat` (standalone hex18, bottom fixed, top loaded
  -1/-4/-16): rc=0, element_interface_intpnt_stress = -9 on ALL 9
  intpnts (target -9 +-1e-5).
- corpus `interface_quad8_hex20.dat` (2 hex20 + quad8 face +
  control_mesh_convert): rc=0, post_point sigzz = -1.005025 vs target
  -1.0+-1e-2 - the same 0.5% offset of the established hex8-interface
  family (interface_quad4_hex8: -1.005025), which comes from the GNU
  convention of shifting the interface duplicates 0.01 along the normal
  (the Pro keeps them coincident; see control_mesh_convert developer
  note).

## Pendiente

- `control_reset_interface(_strain)` zero/reset buffers were enlarged to
  9 slots, but the reset records are exercised only up to ns1=4 in the
  corpus; a hex18 + reset combination has no corpus test yet.
