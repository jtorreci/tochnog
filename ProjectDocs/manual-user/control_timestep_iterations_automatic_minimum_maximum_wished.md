# control_timestep_iterations_automatic_minimum_maximum_wished

## Description

`control_timestep_iterations_automatic_minimum_maximum_wished index
minimum_iterations maximum_iterations wished_iterations` (manual
Professional 6.387): the minimum, maximum and wished (preferred) number of
equilibrium iterations of the automatic time stepping
([control_timestep_iterations_automatic](control_timestep_iterations_automatic.md))
with the same index. Default 2 8 4. The maximum must be 2 or larger (and
not below the minimum); otherwise the input is rejected.

The time step grows when a step needs fewer iterations than the wished
number and shrinks when it needs more; a step that reaches the maximum is
retried with a ten times smaller step (once a step has converged).

## Usage

```
control_timestep_iterations_automatic_minimum_maximum_wished index min max wished
```

## Example

```
control_timestep                       0  0.05 1.
control_timestep_iterations_automatic  0  1.e-4 1.e-6 0.5
control_timestep_iterations_automatic_minimum_maximum_wished  0  2 5 4
```

Tests: `validation-suite/test-2026/control_timestep_iterations_automatic_default1`
(no record: at most 8 iterations per step) and
`control_timestep_iterations_automatic_minimum_maximum_wished1` (2 5 4: at
most 5).

## Differences with the Professional version

- Default 2 8 4 = Professional 02-08-2026 (the GNU had 2 20 4 hardcoded
  before F5 W1 U4, 2026-11-02).
- The rule that grows / shrinks the step differs: the GNU multiplies the
  step by 1 + (wished - iterations)/20 and divides it by 10 after a step at
  the maximum; the Professional grows it by 1.25 (also 1.5625) and divides
  it by 4. The step sequences differ (4x4 von Mises probe: 239 GNU steps,
  166 Professional) (BACKLOG C94, DIVERGENCES class B).
- [control_timestep_adjust_minimum_iterations](control_timestep_adjust_minimum_iterations.md)
  is accepted; the GNU keeps the minimum (as with `-no`).
