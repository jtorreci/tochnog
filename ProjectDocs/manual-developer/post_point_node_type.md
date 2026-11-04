# post_point_node_type

## Implementation (F5 T2 V9)

- **Record:** INTEGER, length 1 (fixed), indexed; `-node` or `-node_start_refined` (`db_error` otherwise).
- **Consumer:** `post.cc`, `parallel_post_point()` reads the element node coordinates from NODE or NODE_START_REFINED (`post_point_search_coord`). The dof values are interpolated as before.

## Verification

`post_point_node_type1` (mesh -follow_material).
