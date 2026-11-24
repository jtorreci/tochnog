# post_point_element_group

## Implementation (F5 T2 V9)

- **Record:** INTEGER, length 1 (fixed), indexed, data_required post_point.
- **Consumer:** `post.cc`, file statics `post_point_search_group / _coord / _eps` set by `post()` (`post_point_search_default()` resets them before every point, before the post_line loop and for every quadrilateral); `parallel_post_point()` skips elements of another ELEMENT_GROUP.

## Verification

`post_point_search1`.
