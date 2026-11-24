# post_count

## Description

`post_count dataitem_name_0 dataitem_name_1 ...` (manual Professional 6.923). Counts the active indices of each listed data item at every step. The counts go to the output record `post_count_result`, in the same order.

## Usage

```
post_count -element -node -geometry_point
```

## Example

`validation-suite/test-2026/post_count_data1.dat`: three unit quad4 elements, two element groups, the right edge pulled. `post_count_result 3 8 1`.

## Differences with the Professional version

None measured (same counts).
