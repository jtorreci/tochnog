# post_node_rhside_ratio_method

## Description

`post_node_rhside_ratio_method method` (manual Professional 6.952). With `-post_node_rhside_free` the `post_node_rhside_ratio` is the free right-hand side itself instead of the free/fixed ratio.

The records it builds on (manual 6.948-6.950):

- `post_node_rhside_fixed` / `post_node_rhside_free`: per principal dof, the signed sum of `node_rhside` over the prescribed / free nodes (for example the summed reaction force and the summed unbalance force).
- `post_node_rhside_ratio`: per principal dof, the mean absolute `node_rhside` of the free nodes divided by that of the fixed nodes; the mean free value when the fixed one is below 1.e-10; the largest value over the dofs (restricted to `post_node_rhside_ratio_dof_type` when given). With `-post_node_rhside_free` it is the mean free value.

The equilibrium iteration control (`control_timestep_iterations_automatic`, `control_timestep_size_automatic_decrease`) uses the ratio of each iteration; the records written to the `.dbs` are those of the final residual of the step.

## Usage

```
post_node_rhside_ratio_method -post_node_rhside_free
```

## Example

`validation-suite/test-2026/post_node_rhside_ratio_method1.dat`: von Mises plasticity, 2 iterations; ratio 1.23884, the mean |velx| right-hand side of the four free nodes. `post_node_rhside_ratio1.dat`: the same input without the record; ratio 12.9442 (vely: 1.17548 / 0.0908132), fixed -3.98859 0.0908132, free 3.98859 -0.0908132.

## Differences with the Professional version

None measured: the values equal the Professional 02-08-2026 (relative 1e-4, the residual of the plastic iterate). The manual says "size" without defining it; the mean absolute value is the measured Professional rule (BACKLOG C110, resolved 2026-11-05).
