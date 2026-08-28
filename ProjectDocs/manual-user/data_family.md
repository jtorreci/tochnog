# data_activate / data_delete / data_ignore + solver & misc aliases

## data_activate / data_activate_time

Time-gated variant of `control_data_activate` (manual Professional
6.397/6.398): with `-no` the listed data items are de-activated
(their records deleted) from the `data_activate_time` time point on
(default: at the start of the calculation). `-yes` is a no-op
(records are active by default; the GNU de-activation is destructive).

## data_delete / data_delete_time

Time-gated variant of `control_data_delete` (manual 6.399/6.400):
delete records by index, range `-ra ... -ra` or `-all`; elements/nodes
are deleted with all their records. Evaluated from the
`data_delete_time` point on (default: at the start).

## data_ignore

`data_ignore -<item>` (manual 6.401): every record with that item name
is SKIPPED at input time. NOTE: the `data_ignore` record must appear
BEFORE the records it ignores (it is evaluated as the input is read).
Typical use: disable `-print_apply`-style switches to get full output.

## Aliases (db_number)

- `control_solver index <type>` → `control_options_solver`
  (per-timestep solver: `-diagonal`, `-matrix_iterative_bicg`, ...).
- `control_solver_bicg_error index error` →
  `control_options_solver_bicg_error`.
- `axisymmetric index -yes` → `group_axisymmetric`.
- `bounda_print_mesh_dof*` → `print_mesh_dof*`.

## control_solver_bicg_stop / partials

`control_solver_bicg_stop index -no`: the calculation CONTINUES when
the BI-CG solver does not converge (warning; default/`-yes` stops).
Since the A+B solver fix (2026-08-28) this also covers the breakdown
exits (the solver no longer reports breakdown/stagnation as success —
see the solver note in the developer manual): a non-converging or
breakdown solve stops the calculation with RC≠0 by default, and
continues with the current solution (warning) with `-no`. The `index`
must match the control index of the `control_timestep` block of the
model (e.g. `control_solver_bicg_stop 5 -no` for `control_timestep 5
...`); with a non-matching index the record is silently ignored.
Registered without behaviour (documented): `control_solver_bicg_restart`
(no restart in the GNU bicg), `control_solver_matrix_save`,
`control_solver_pardiso_ordering`/`_out_of_core` (PARDISO not compiled
in).

## print_mesh_dof (+ _geometry / _values)

One-shot dump (at the first evaluation) of node coordinates and the
listed dof values (all dofs when none listed) to `print_mesh_dof.dat`;
nodes can be restricted to a geometry via `print_mesh_dof_geometry`.
`print_mesh_dof_values` is registered for compatibility (not used by
the dump).

## Example

```
data_activate 1  -bounda_force -no
data_activate_time 1  0.1
```
(test `dsmall`: the load disappears at t=0.1 and the quasi-static
column relaxes; also exercises the solver/axisymmetric aliases and the
print_mesh_dof dump)

```
data_ignore -bounda_time
```
(test `dignore`: declared BEFORE the records, the load defaults to 0)
