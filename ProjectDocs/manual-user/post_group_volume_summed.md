# post_group_volume_summed

## Description

Output record (manual Professional 6.937): the summed volume of the elements of each element group (group 0 first), empty elements excluded. The GNU computes it on demand, when a `target_item` or a `control_print` asks for it, like `post_force_edge_summed`.

## Usage

```
target_item 3 -post_group_volume_summed 0 1
```

## Example

`validation-suite/test-2026/post_count_data1.dat`: three unit quad4 elements, two element groups, the right edge pulled. `post_group_volume_summed 2 1`.

## Differences with the Professional version

The Professional fills it in every run (a `post_global` item, default `-yes`); the GNU only when it is asked for, so other database files do not change. Same values.
