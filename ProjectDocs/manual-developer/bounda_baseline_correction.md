# bounda_baseline_correction — developer

## Files and functions

- `database.cc` — `BOUNDA_BASELINE_CORRECTION` (DOUBLE_PRECISION, 2, no
  index), `BOUNDA_BASELINE_CORRECTION_PARAMETERS` (DOUBLE_PRECISION, 3,
  class BOUNDA); enums after `BOUNDA_TIME_SMC_UNITS` in `tochnog.h` /
  `tochnog-mod.h`.
- `bounda.cc`:
  - `bounda_baseline_table()`: the value of a time-value bounda_time table
    (same interpolation as the time loop of `bounda()`, 0 outside).
  - `bounda_baseline_load()`: p(tau) at time_total; fits and PUTs the
    parameters when the record is missing (normal equations of the
    velocity least squares, Gauss elimination with partial pivoting).
  - `bounda()`: after the bounda_time_until scaling, for a record with an
    acceleration dof label, `bounda_baseline` is computed once per record
    and added to `load` only in the `-acc*` branch (acceleration dof and
    the velocity bound v_new = v_old + a*dt).

## Implementation details

Measured on the Professional 02-08-2026 (1D truss, `bounda_dof 10 1 -accx`,
`bounda_time 10 0 1 .5 1 1 1 2 0`): the parameters use the normalized time
tau and the same record layout; the Professional adds p(tau_k) to the
bounda_time points with time_start < t_k <= time_end and interpolates, so
with (-1 0 0) the applied acceleration is 1-2t on [0,0.5] and the velocity
drifts to 0.2. Its fit does not depend on the time step; [0,2] gives
(-0.8182607915309 -1.375947429918 2.186213875160), not reproduced by any
least-squares of velocities that was tried, and a = t on [0,2] gives
(-2 0 0). GNU follows the manual (owner decision 2026-11-22).

## History

- F5 T5 U9 (2026-11-22): added.

## Hardcoded parameters / pending refactorings

- The fit uses the time step of the first step in the interval; a later
  change of the step size does not refit.
- The fit assumes rest (v = 0) at time_start.
- `bounda_time_smc`: `bounda()` passes the whole SMC table with table = 1.
