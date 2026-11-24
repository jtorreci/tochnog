# post_force_edge_summed

## Description

`-post_force_edge_summed` is a read-only post result: the TOTAL force
following from the `-force_edge` records, integrated over the edges, per
space direction (`number_of_space_dimensions` values). It is one of the
`post_global` items of the Professional manual (6.936), which the GNU
computes on demand.

## Uso

It can be requested in two ways:

1. As a `control_print` item (prints the stored record to the output):

```
control_print 0  -time_current -node_dof -post_force_edge_summed
```

2. As a `target_item` value for automated checks:

```
target_item 1 -post_force_edge_summed 0 1
target_value 1 10. 1.e-4
```

The record stores one value per space direction (x, y, z): the target
above checks the y component (`1` = second value) against 10 (a uniform
vertical edge load of 1.0 per unit length over an edge of length 10).

## Notes

- The value is recomputed at every step close, so time-dependent loads
  are reflected.
- It is the total APPLIED load of the `force_edge` records (all indices),
  including their geometry/element/node restrictions, time and spatial
  factors - the same integration that distributes the traction to the
  element nodes.

## Related

- `force_edge` (alias of `force_element_edge`) — the edge load records.
- `control_mesh_convert`, `quad8` — quadratic element conversion.
