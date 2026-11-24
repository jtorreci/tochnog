# timestep_predict_velocity / timestep_iterations_automatic_apply

## timestep_predict_velocity

`-yes` asks the solver to use the PREVIOUS timestep velocities as the
prediction (initial guess) for the current timestep, instead of the
zero velocity that the no-inertia/no-convection default uses (manual
Professional 6.1089). Typically wanted in eulerian calculations.

**GNU status: PARTIAL.** The record is registered and accepted, but
the GNU's staggered scheme does not map to this semantics: its
iterative solve always starts from zero and the result is ADDED to the
unknowns as an increment (`node_dof_new += solution`). Starting the
solve from the previous velocities would double-count them. The
record is therefore behavior-neutral (verified: identical results with
and without it, test `tslv_part`).

## timestep_iterations_automatic_apply

`-no` neglects every `control_timestep_iterations_automatic` record
(manual Professional 6.1090): the adaptive iteration/time-step
machinery is switched off as if the records were absent. The default
`-yes` keeps them active.

Wired in `top.cc` right before the automatic branch: when `-no` is
set, the record is not read and the plain defaults apply (verified:
identical results with the automatic record plus `-no` and without
both, test `tslv_auto`).

## Example

```
control_timestep_iterations_automatic 1 1.e-3 0.05
timestep_iterations_automatic_apply -no
```
