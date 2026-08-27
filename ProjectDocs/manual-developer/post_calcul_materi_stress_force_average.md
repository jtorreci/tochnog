# post_calcul_materi_stress_force_average

## Implementación

- Registered in `database.cc` (data class POST, `no_index = 1`):
  INTEGER, `data_length = 1` (fixed), default `-yes` when absent
  (manual 6.908). Value validated (`-yes`/`-no`) in
  `post_calcul_materi_stress_force_validate()` (`calcul_force.cc`).
- **Implemented in lot 2 (2D)**: with `-yes` (default) the quad9
  middle-plane nodes (the 3 nodes NOT on either end face) receive the
  average of the forces/moments of the two opposing end faces
  (`msf_element_contribution_2d` in `calcul_force.cc`). With `-no`
  those nodes receive nothing (0) and are NOT marked as averaged.
- The averaged nodes are flagged in the per-node record
  `POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE` (-YES/-NO) and the
  `-primary` method of
  [`control_print_materi_stress_force`](control_print_materi_stress_force.md)
  skips them through the `msf_node_is_averaged()` hook
  (`print_materi_stress_force.cc`). With `-no` nothing is averaged, so
  `-all` and `-primary` print the same lines (verified by
  msf_quad9_noavg).

## Diseño / decisiones

- The middle plane is defined as the nodes NOT on the two end faces
  (for the canonical beam the mid-line between the faces). The average
  is taken on the SIGNED physical values (nor, she, mom) BEFORE the
  plot components are built, so the averaged node carries the value at
  the element middle (e.g. the moment at mid-element = the mean of the
  two face moments, EXACT - verified by msf_beam2d/msf_quad9).
- `-primary` skips a node when it received ANY averaged contribution
  (mixed primary/averaged nodes in irregular meshes are excluded -
  documented conservative choice).

## Pendiente

- The 3D (hex27) averaging lands in lot 3 together with the 3D
  integration.
