# control_timestep_iterations_automatic (variable 2-3 values)

## Description

`control_timestep_iterations_automatic` (manual Professional 6.386)
selects the automatic sub-stepping of the timestep. The record is
**variable-length**: the Professional parser accepts a short list and
stops at the next keyword, so both of these forms are valid:

```
control_timestep_iterations_automatic <index> ratio_criterium minimal_timestep maximum_timestep
control_timestep_iterations_automatic <index> ratio_criterium maximum_timestep
```

| Value | Meaning |
|-------|---------|
| ratio_criterium | the post_node_rhside_ratio criterion that stops the equilibrium iterations |
| minimal_timestep | the lower bound for the automatically reduced timestep (3-value form only) |
| maximum_timestep | the upper bound for the automatically grown timestep |

The slot mapping is length-aware:

- **3-value form** (manual 6.386): `[0]=ratio_criterium`,
  `[1]=minimal_timestep`, `[2]=maximum_timestep`.
- **2-value form** (legacy GNU and every input in the test corpus):
  `[0]=ratio_criterium`, `[1]=maximum_timestep`.

`minimal_timestep` is registered but not consumed by the current GNU
automatic scheme (the scheme only reduces and grows the step; it has no
lower bound of its own).

## Usage

Three values (Professional manual form):

```
control_timestep                200  1.e-6 1.0
control_timestep_iterations_automatic  200  1.e-3 1.e-6 1.e-3
```

Two values (legacy form, used by the whole sfnet corpus):

```
control_timestep                       5  1.e-7 80.0e-6
control_timestep_iterations_automatic  5  0.0001 1.0e-6
```

## Notes

- The automatic sub-stepping (`control_timestep_iterations_automatic`
  without `_apply -no`) shrinks the timestep when the equilibrium
  criterion is not met and retries, and grows it when the step converges
  below the wished iteration count.
- The record must not swallow the following keyword: a short list ends
  at the next keyword token (`control_options_solver`,
  `control_mesh_delete_geometry`, `target_item`, ...), which is what the
  variable-length registration guarantees.
