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
- dof.cc `parallel_new_dof_before()`: node_force applied to the node
  rhs (force_point sign convention); node_inertia filled per step with
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

- node_force_inertia: inertia record -19.83 ~ m*g (node_force -30
  active; m=2, g=-10). node_mass fix exercised (2D values parse).
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
