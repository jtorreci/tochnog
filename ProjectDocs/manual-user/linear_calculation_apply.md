# linear_calculation_apply

## Description

`linear_calculation_apply switch` (manual Professional 6.797, data part,
no index, default `-no`): with `-yes` the nonlinearities of the input are
skipped, so the calculation is linear. Useful for testing and problem
search: the printed `-post_node_rhside_ratio` should be very small.

The database is rewritten once, at the start of the calculation:

- any `control_plasti_apply` is deleted and `plasti_apply` is set to
  `-no`; the `group_materi_plasti_*` and `group_interface_materi_plasti_*`
  records are deleted (Professional 02-08-2026, manual silent);
- `mesh` (`options_mesh`) is set to `-fixed_in_space`;
- every `group_*_memory` record is set to `-total_linear` when
  `materi_displacement` is initialised, else to `-updated_linear` (groups
  without a memory record keep the default);
- any `dependency_item` / `dependency_diagram` of a `group_*` item on a
  dof is deleted (a dependency on `-time_current` is kept);
- `group_materi_elasti_hardsoil` -> `group_materi_elasti_young` E50ref
  (the Poisson ratio is kept);
- `group_materi_elasti_young_power` -> `group_materi_elasti_young` E0;
- `group_materi_plasti_hypo_wolffersdorff` -> `group_materi_elasti_young`
  hs and `group_materi_elasti_poisson` 0.2;
- `group_spring_stiffness_nonlinear` -> `group_spring_stiffness` with the
  stiffness at strain 0 (linear interpolation in the table; outside the
  table the nearest end value);
- deleted: `group_interface_gap`, `group_materi_damage*`,
  `group_materi_failure*`, `group_truss_rope`, `contact_*`,
  `groundflow_seepage_*`; `contact_apply -no` is added (Professional).

Not in the GNU (no such record): `group_materi_elasti_polynomial`,
`group_groundflow_nonsaturated`, `group_interface_materi_hardening`,
`group_interface_materi_elasti_stiffness_tangential_diagram`.

## Usage

```
linear_calculation_apply -yes
```

Professional example: `test/tutorial/tutorial_4/tutorial_4.dat` (inside
`start_if test_calculation`). Tests: `validation-suite/test-2026/
linear_calculation_apply1` (Mohr-Coulomb removed: sigyy -10 instead of
the plastic -5.07), `_apply2` (hardsoil: -5), `_apply3` (young_power: -7),
`_apply4` (wolffersdorff: -4.1666667), `_apply5` (the rewritten records in
the `.dbs`).

`tutorial_4` (Professional example, with `test_calculation true`): the
GNU runs it only without `post_calcul ... -force` and without
`control_print_gid` (BACKLOG C97). With `linear_calculation_apply -yes`
the GNU ends with `post_node_rhside_ratio` 4.6e-13, as the manual
expects for a linear calculation; the Professional 02-08-2026 ends with
3.9e-2 (its rewritten database is the one listed above, so the remaining
nonlinearity is elsewhere). The nonlinear run already differs between
the two programs (post_node_result 0 velx -2.29 GNU, -4.17 Professional;
C97), so the linear runs are not compared further.

## Differences with the Professional version

- `group_materi_plasti_hypo_wolffersdorff`: the Professional 02-08-2026
  takes the first value (phi) as Young's modulus instead of hs; the GNU
  follows the manual (DIVERGENCES `LINCALC-WOLFF`, class D).
- `contact_*`: the Professional adds `contact_apply -no` but keeps the
  other contact records; the GNU deletes them as the manual says and adds
  `contact_apply -no` too (DIVERGENCES `LINCALC-CONTACT`, class D; the
  contact is off in both).
