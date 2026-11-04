# post_data

## Implementation (F5 T2 V9)

- **Records:** `post_data` INTEGER (triplets), `post_data_factor` DOUBLE (data_required post_data), `post_data_result` DOUBLE length 1; all indexed, class POST.
- **Fill:** `post.cc`, `post_count_data_volume()`: the record length must be a multiple of 3 (`db_error`). DOUBLE items go through `data_item_value()` (bounda.cc; a negative number is a dof label), INTEGER items are read directly. Factor 1 without `post_data_factor`; a shorter factor record stops with an error.

## Verification

`post_count_data1` (target 6 = Professional), `post_data_factor_error1`.
