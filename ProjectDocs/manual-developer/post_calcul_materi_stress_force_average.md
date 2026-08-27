# post_calcul_materi_stress_force_average

## Implementación

- Registered in `database.cc` (data class POST, `no_index = 1`):
  INTEGER, `data_length = 1` (fixed), default `-yes` when absent
  (manual 6.908). Value validated (`-yes`/`-no`) in
  `post_calcul_materi_stress_force_validate()` (`calcul_force.cc`).
- Semantics (consumed in lots 2/3): only quad9/hex27; `-yes` averages
  the forces/moments of the two opposing end faces on the middle-plane
  nodes. The `-primary` method of
  [`control_print_materi_stress_force`](control_print_materi_stress_force.md)
  skips those averaged nodes via the `msf_node_is_averaged()` hook in
  `print_materi_stress_force.cc` (lot 1: returns 0 for every node, so
  `-all` and `-primary` are identical).

## Pendiente

- The averaging itself lands in lot 2/3 together with the flag that
  marks the averaged (non-primary) nodes.
