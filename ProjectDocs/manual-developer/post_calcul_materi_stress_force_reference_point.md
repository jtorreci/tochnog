# post_calcul_materi_stress_force_reference_point

## Implementación

- Registered in `database.cc` (data class POST, `no_index = 1`):
  DOUBLE_PRECISION, `data_length = DATA_ITEM_SIZE`,
  `fixed_length = 0` (variable: `ngroups * ndim` values, one point per
  element group).
- Validated in `post_calcul_materi_stress_force_validate()`
  (`calcul_force.cc`): total length must be `ngroups * ndim` (clear
  error otherwise); in 3D the record is REQUIRED (error when absent),
  in 2D a warning is issued and the documented default (0,0) is used.
- Consumed in lots 2/3 to orient the forces/moments consistently
  outwards/inwards in the thickness direction (the far/outer nodes of
  [`post_calcul_materi_stress_force_outer`](post_calcul_materi_stress_force_outer.md)
  are also selected with it).

## Implementación (lote 2, 2D)

- Implemented in `msf_element_contribution_2d` (calcul_force.cc): the
  thickness direction is t = (element centroid - reference_point)
  normalized in-plane; the two end faces are the sides whose exterior
  normals are most perpendicular to t (smallest |n*t|); the plot
  vectors point along t (outward/inward per the reference point).
  Degenerate case (reference point at the element centroid): warning
  + the element is skipped (no forces). Pending for 3D (lot 3).
