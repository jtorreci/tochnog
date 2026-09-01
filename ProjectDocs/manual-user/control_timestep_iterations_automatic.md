# control_timestep_iterations_automatic (3 values)

## Description

`control_timestep_iterations_automatic` (manual Professional 6.386)
takes **three values** in the Professional:

```
control_timestep_iterations_automatic <index> ratio_criterium minimal_timestep maximum_timestep
```

| Value | Meaning |
|-------|---------|
| ratio_criterium | the post_node_rhside_ratio criterion that stops the equilibrium iterations |
| minimal_timestep | the lower bound for the automatically reduced timestep |
| maximum_timestep | the upper bound for the automatically reduced timestep |

The GNU previously read only two values (ratio, maximum); the missing
`minimal_timestep` made the maximum land in the second slot, so tests
that give all three values (e.g. the slope family:
`control_timestep_iterations_automatic 200 1.e-3 1.e-6 1.e-3`) read a
wrong maximum and either diverged or never advanced.

## Usage

```
control_timestep                200  1.e-6 1.0
control_timestep_iterations_automatic  200  1.e-3 1.e-6 1.e-3
```

## Notes

- The automatic sub-stepping (`control_timestep_iterations_automatic`
  without `_apply -no`) shrinks the timestep when the equilibrium
  criterion is not met and retries.
- slope_classical_numerical now runs to completion with the correct
  maximum (1.e-3); the collapse time still differs from the target
  (physics, see SEGUIMIENTO).
