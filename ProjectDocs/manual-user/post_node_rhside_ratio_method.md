# post_node_rhside_ratio_method

## Description

`post_node_rhside_ratio_method method` (manual Professional 6.952). With `-post_node_rhside_free` the `post_node_rhside_ratio` is the free right-hand side itself (summed over the dofs of a dof group, the largest group kept) instead of the free/fixed ratio.

## Usage

```
post_node_rhside_ratio_method -post_node_rhside_free
```

## Example

`validation-suite/test-2026/post_node_rhside_ratio_method1.dat`: von Mises plasticity, 2 iterations; the ratio equals the free velocity right-hand side.

## Differences with the Professional version

The GNU `post_node_rhside_fixed` / `_free` hold the largest absolute value per dof, the manual and the Professional the sum, so the values differ (BACKLOG C110).
