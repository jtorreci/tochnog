# Data item type check of run-time record names

## What it does

Some records name ANOTHER data item that is read while the calculation
runs:

- `post_node index data_item_name -sum|-average|-moment selector`
- `post_integrate index data_item_name data_item_index number`
- `control_repeat_save index data_item_name data_item_index number ...`
- `control_repeat_until_item index start criterium data_item_name
  data_item_index number`
- `control_print_element index data_item_name`

These features work on floating-point values. The named data item must
therefore be a DOUBLE PRECISION record (for example `node_dof`,
`node_rhside`, `post_point_dof`, `element_truss_force`). Naming an
INTEGER record (for example `node_bounded` or `element_group`) now stops
the calculation with:

```
Error: a DOUBLE PRECISION data item is required here, not the INTEGER item: node_bounded
Error detected for data item : post_node
Error detected for record : 0
```

## Why it is useful

Before 2026-10-17 (F4 T6f, BACKLOG C54) an INTEGER record named there was
read into a one-value scratch buffer: the extra values overwrote memory,
so the result was undefined (garbage, a crash, or an endless
`control_repeat_until_item` loop).

## Example

```
post_node 0 -node_bounded -sum -all      ( stops with the error above )
post_node 0 -node_rhside -sum -all       ( valid )
```

## Differences with the Professional version

The Professional 02-08-2026 returns zeros for `post_node` and
`post_integrate` on an INTEGER record and stops on
`control_repeat_save` (`Error detected for data item : node_bounded`).
The GNU stops with an input error in every case (DIVERGENCES C54-TYPE).
`control_print_element` already refused INTEGER records before the check.
