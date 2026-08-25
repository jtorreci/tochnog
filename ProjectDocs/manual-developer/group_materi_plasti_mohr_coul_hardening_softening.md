# group_materi_plasti_mohr_coul_hardening_softening

## Implementación

- **Ley**: new block in `plasti_rule()` in `plasti.cc` (CRLF file), right
  before the classic `group_materi_plasti_mohr_coul` block. It reads the
  group record `GROUP_MATERI_PLASTI_MOHR_COUL_HARDENING_SOFTENING`
  (DOUBLE, length 7) and implements both `GET_YIELD_RULE` and
  `GET_FLOW_RULE` (phi_flow) for `plasti_type == -NONE` (first plastic law
  active) and for the explicit type. The flow gradient is provided by the
  standard driver via central finite differences, so the law only defines
  `f` and `f_flow` — exactly like the classic block.
- **Interpolation**: linear in `ratio = kappa/kappa_shear_crit` clamped to
  `[0,1]`; `kappa` is read from the node dof `new_unknowns[kap_indx]`
  (`materi_plasti_kappa`). At `kappa_shear_crit <= 0` the law falls back
  to the final values (`ratio = 1`).
- **Surface**: same as classic MC — `matrix_eigenvalues` on the stress
  tensor, `f = 0.5*(sig_max-sig_min) + 0.5*(sig_max+sig_min)*sin(phi) -
  c*cos(phi)`. Boundary factor (`plasti_on_boundary`) applies to phi and
  phi_flow as in the classic block.
- **Keywords** (data_class MATERI) registered in `database.cc`:
  - `group_materi_plasti_mohr_coul_hardening_softening` (DOUBLE, length 7,
    required GROUP_TYPE): `phi_0 c_0 phi_flow_0 phi_1 c_1 phi_flow_1
    kappa_shear_crit`.
- **New enum**: `GROUP_MATERI_PLASTI_MOHR_COUL_HARDENING_SOFTENING` in
  `tochnog.h` / `tochnog-mod.h` (kept in sync, same order).

## Física

`kappa` is a node dof accumulated by the standard driver in `stress.cc`
(`new_kappa = old_kappa + sqrt(0.5*inc_epp:inc_epp)`), assembled into the
element RHS in `materi()` (evolution equation `volume*h*(new_kappa -
old_kappa)/dtime`) and resolved by the solver with the lumped inertia from
`general()` (dof type `MATERI_PLASTI_KAPPA`, `inertia = 1.`). The law reads
the nodal `kappa` interpolated to the integration point; with softening the
surface moves DOWN inside the step, so the plastic return needs several
global iterations (`control_timestep_iterations`) to converge to the current
surface — with the default single iteration the stress stays above it
(weak kappa coupling, observed 7.55 -> 81.9 after calibration).

## Validation

Test `mmchs_soft` in `validation-suite/test-2014/` (registered in
`scripts/build_safe.sh`): uniaxial tension, `phi = 0`, `c_0 = 80`,
`c_1 = 20`, `kappa_shear_crit = 0.5`, `dt = 0.02` (20 steps),
`control_timestep_iterations 8`. Measured `kappa = 0.318` -> ratio 0.636
-> `c = 41.8` -> analytic `sig_t = 83.6`; measured `sigxx = 81.9` (98%).
A/B without softening (`c_1 = c_0 = 80`) gives `159.2 = 2*c_0`: the test
discriminates. The earlier single-step rig diverged (7.55) because kappa
jumped past kappa_crit and the surface moved faster than the return.

## Pendiente / gotchas

- GOTCHA (test, not code): the target/post `post_point_dof` basename for
  the kappa dof is `-kap`, NOT `-materi_plasti_kappa` — the latter is the
  initia keyword and `array_member()` over `dof_label` misses it, producing
  an out-of-range read (DBL_MAX) instead of an error.
- CRLF: `plasti.cc` is CRLF; the new block was added preserving CRLF.
