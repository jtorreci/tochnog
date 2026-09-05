# post_apply

## Implementation

- Keyword `POST_APPLY` (INTEGER, no_index) registered in
  `database.cc`; enum appended in `tochnog.h`/`tochnog-mod.h` (same
  enum, same order). Parse-only: consumption of `-no` PENDING (the
  per-step post evaluation sites would need the global gate; no corpus
  test exercises `-no`).

## Pending

- Gate the post-processing evaluation points (post_point/post_line/
  post_quadrilateral updates, node_rhside_ratio exempted) on the
  switch.
