# post_calcul_materi_stress_force_outer

## Implementación

- Registered in `database.cc` (data class POST, `no_index = 1`):
  INTEGER, `data_length = 1` (fixed), default `-no` when absent
  (manual 6.915). Value validated (`-yes`/`-no`) in
  `post_calcul_materi_stress_force_validate()` (`calcul_force.cc`).
- Semantics (consumed in lots 2/3): `-yes` restricts the forces and
  moments to the nodes at the outer sides (furthest from the reference
  point).

## Pendiente

- The outer-node selection lands in lot 2/3.
