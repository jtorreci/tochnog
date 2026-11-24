# post_count

## Implementation (F5 T2 V9)

- **Records:** `database.cc`, `post_count` and `post_count_result` INTEGER, variable length, no index, class POST.
- **Fill:** `post.cc`, `post_count_data_volume()` (called by `post()` before `post_global()`): for every listed item, `db_max_index` and a `db_active_index` count.

## Verification

`post_count_data1` (grep of the .dbs line in `scripts/build_safe.sh`; = Professional).
