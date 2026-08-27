# post_calcul_materi_stress_force_outer

## Implementación

- Registered in `database.cc` (data class POST, `no_index = 1`):
  INTEGER, `data_length = 1` (fixed), default `-no` when absent
  (manual 6.915). Value validated (`-yes`/`-no`) in
  `post_calcul_materi_stress_force_validate()` (`calcul_force.cc`).
- Semantics (consumed in lots 2/3): `-yes` restricts the forces and
  moments to the nodes at the outer sides (furthest from the reference
  point).

## Implementación (lote 2, 2D)

- Implemented in `msf_element_contribution_2d` (calcul_force.cc): with
  `-yes` only the PRIMARY nodes at the maximum distance from the
  reference point receive values; the averaged nodes receive nothing
  (documented decision). Pending for hex27 (lot 3).
