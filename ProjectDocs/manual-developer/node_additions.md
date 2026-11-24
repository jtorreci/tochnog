# node_* additions (developer, Sprint 8)

## Files and functions

- Enums (tochnog.h / tochnog-mod.h, 1053 in sync): NODE_CONVECTION_APPLY,
  NODE_DYNAMIC_PRESSURE, NODE_FORCE, NODE_INERTIA, NODE_MESH, NODE_SLIDE,
  NODE_STATIC_PRESSURE, NODE_TOTAL_PRESSURE.
- database.cc: registrations. CRITICAL: every input-able NODE-class item
  needs `version_all = 1` — mesh_has_changed DELETES all NODE-class
  records without a version in VERSION_NORMAL at every mesh change, so
  without it the record silently disappears after the input (found the
  hard way: GET_IF_EXISTS never fired; debug probe showed max_index=-1).
  Also fixed the PREEXISTING `data_length[NODE_MASS] = 1` bug (was 1;
  node_damping/stiffness correctly use ndim): node_mass in 2D/3D
  consumed only the first component and the parser choked on the rest —
  node_mass was unusable since the GNU origins.
- dof.cc `parallel_new_dof_before()`: node_force ADDED to the node
  rhs of the ndim velocity equations (`node_rhside += node_force`:
  node_rhside holds external - internal; same convention as
  `force_point_calculate()` and the force_element_edge loads). C28
  (2026-10-06): it was subtracted (Pro 1D probe velx +1, GNU -2). The
  "factor 2" is not a double application: the reversed load compressed
  the unit bar (E = 1) by exactly its length, a singular zero-length
  state (now: `node_force 2 -1.0` GNU -2 / Pro -1, `-0.5` both -0.5,
  `force_point` at the node -1 both -1; out of scope, degenerate). `delete.cc` `mesh_cut()`
  dropped its compensating negation in the same work unit (the stored
  substitute now equals the Professional's record). node_inertia filled per step with
  m*(a+g). GOTCHA: the PUT must go to VERSION_NEW — the end-of-step
  db_version_copy(NEW->NORMAL) overwrites NORMAL with the input values
  (a NORMAL write was silently lost; found with a debug pri). The PUT
  only fires for records the user declared (parallel-loop safe: never
  allocates).
- calcul.cc `calculate_operat()`: the three pressure overrides after
  the phreatic calculation of -total/-static/-dynamic (added ldum to
  the local declarations).
- slide.cc: node_slide membership checked before geometry(); additive
  (in_geometry |= member). Note the preexisting quirk: the SLIDE_FRICTION
  read uses `inod` instead of `islide` (works because inod==0==the usual
  index at that point); left untouched.

## Verification / limitations (suite 85/85)

- node_force_inertia: inertia record ~ m*g (node_force -30 active;
  m=2, g=-10): -19.75 with the old reversed sign (target -19.83 +- 0.2
  calibrated on it), -20.109 since C28 (the downward force adds a
  downward acceleration; the Professional cannot run the input:
  `node_mass` rejected, `node_inertia` wants one value per dof). The
  stale target (BACKLOG B12) was retargeted to -20.11 +- 0.2 on
  2026-10-06. node_mass fix exercised (2D values parse).
- C28 sign (2026-10-06): `test-2026/node_force_sign1` (quad4 bar, disx
  2; old -2) and `test-2026/mesh_cut_node_force_sign1` (Pro mesh_cut_1 +
  target `node_force 6` = +1; old -1): RED old / GREEN new. Probes 1D
  bar2 / 2D quad4 / 3D hex8 equal the Professional 02-08-2026.
- node_pressure: static override 7.5 exact vs calculated -2.0 (A/B).
- node_slide: smoke test (clone of slide_axi with the geometry shifted
  off the node and membership by record). LIMITATION (documented): the
  original slide_axi ALSO reads rhside 0 — its ±0.05 target passed
  trivially; a physically measurable slide-friction verification needs a
  dedicated model (the post-solve node_rhside carries the reaction, not
  the friction term the slide added before the solve).
- node_mesh / node_convection_apply: registered + parsed but WITHOUT
  behaviour (the `mesh` per-node record has no GNU consumer; the
  per-node convection switch would need options_convection surgery).
  Documented as partial.
- node_bounded_index / node_geometry_present: print-only records,
  still PENDIENTE (no consumer value in the GNU).

## The 30-minute debugging trail (for the record)

node_mass data_length → parallel-safe PUT → VERSION_NEW vs NORMAL →
version_all deletion: four distinct failure modes, each found with a
one-line `pri` probe. All documented here so the next node_* item starts
from the right pattern.
