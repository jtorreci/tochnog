# post_calcul_static_pressure_height — developer

Enum + registration (DOUBLE, variable DATA_ITEM_SIZE, no_index 1, class
POST). Consumed in groundflow_phreatic_coord() (groundfl.cc): when no
phreatic level applied to the node (found == 0) and the record exists, the
regions (triples coord_min/coord_max/height_ref along the vertical
coordinate) are searched; on a match with the element-group restriction
(post_calcul_static_pressure_height_element_group, INTEGER no_index 1)
static_pressure = force_gravity[ndim-1]*dens*(height_ref - coord_y),
total_pressure = pres_dof + static_pressure (same pattern as the
phreatic-level branch), found = 1. Group check uses
node_attached_element_groups() (group.cc); -ALL values pass any node.
Verified against the Professional .dbs: ground13 st_pres = -1220 EXACT.
ground13 now passes entirely: the bounda_dof -topres prescription (total
pressure) converts the prescribed p_total into the pres dof value through
the SAME static the machinery adds here (pres_dof = load - static when the
node is inside a level/height region), so the resulting to_pres matches the
Professional (see manual-developer/bounda_dof.md).
