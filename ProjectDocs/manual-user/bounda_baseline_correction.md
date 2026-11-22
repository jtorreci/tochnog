# bounda_baseline_correction

## Description

`bounda_baseline_correction time_start time_end` (manual 6.19, no index)
and `bounda_baseline_correction_parameters index c0 c1 c2` (manual 6.20,
index = bounda index): baseline correction of a prescribed acceleration
(`bounda_dof ... -accx/-accy/-accz` with a `bounda_time` table), to
suppress the drift of the velocity.

## Effect on the calculation

- For t in [time_start, time_end] the parabola
  p(tau) = c0 + c1*tau + c2*tau^2, tau = (t - time_start)/(time_end - time_start),
  is added to the prescribed acceleration. Outside the interval the
  acceleration is not changed.
- Only the acceleration dofs of the bounda records are corrected; velocity,
  displacement and other dofs of the same record keep the table value.
- Without `bounda_baseline_correction_parameters` for the bounda index the
  coefficients are determined once, at the first step in the interval:
  they minimize the sum of the squared velocities of the time steps in
  [time_start, time_end] (step size = the time step of that first step,
  velocity v_k = v_k-1 + a_k*dt, starting from rest at time_start). The
  result is written as `bounda_baseline_correction_parameters` in the
  database.
- A `bounda_baseline_correction_parameters` record given in the input is
  used as is (no fit).
- Errors: "Error: bounda_baseline_correction needs at least 3 time steps
  between time_start and time_end for bounda index"; "Error:
  bounda_baseline_correction needs a bounda_time table of time-value pairs
  for bounda index" (bounda_time_file, bounda_time_increment, bounda_sine or
  user tables without given parameters).
- `bounda_time_smc` is not implemented in this build.

Differences with the Professional 02-08-2026: DIVERGENCES
BASELINE-CORRECTION-SIGNAL (it adds p to the bounda_time points after
time_start and its fit gives other coefficients).

## Example

```
bounda_dof  10  1 -accx
bounda_time 10  0. 1.  0.5 1.  1. 1.  2. 0.
bounda_baseline_correction 0. 1.
```

gives `bounda_baseline_correction_parameters 10 -1 0 0`: the velocity
stays 0 up to t = 1.

Self-checking inputs in `validation-suite/test-2026`:
`bounda_baseline_correction1` (constant signal, velocity record untouched),
`bounda_baseline_correction2` (fit over [0,2]),
`bounda_baseline_correction3` (given parameters),
`bounda_baseline_correction_error1` (too few steps).
