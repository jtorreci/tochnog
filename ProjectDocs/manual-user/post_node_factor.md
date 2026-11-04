# post_node_factor

## Description

`post_node_factor index factor` (manual Professional 6.946): multiplies the result of `post_node` with the same index (`post_node_result`). Default 1.

## Usage

```
post_node 1 -node_dof -sum -geometry_line 2
post_node_factor 1 3.
```

## Example

`validation-suite/test-2026/post_count_data1.dat`: three unit quad4 elements, two element groups, the right edge pulled. `post_node_result 1` velx 0.6 (3 x 0.2).

## Differences with the Professional version

None measured.
