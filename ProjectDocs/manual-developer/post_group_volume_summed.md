# post_group_volume_summed

## Implementation (F5 T2 V9)

- **Record:** DOUBLE, variable length, no index, class POST.
- **Fill:** `post.cc`, `post_count_data_volume()`, only when `post_item_requested()` finds a `target_item` or `control_print` for it (the helper extracted from `post_force_edge_summed_calculate()`, which now uses it too). ELEMENT_VOLUME summed per ELEMENT_GROUP (default 0), elements with ELEMENT_EMPTY other than -no / -front skipped.

## Verification

`post_count_data1` (2 1 = Professional).

## Pending

- The other Professional `post_global` items (`post_group_summed`, `post_element_volume_summed`, ...) are still not implemented.
