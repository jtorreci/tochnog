# post_point_eps_iso

## Implementation (F5 T2 V9)

- **Record:** DOUBLE, length 1 (fixed), indexed.
- **Consumer:** `post.cc`, `parallel_post_point()` passes `post_point_search_eps` as the `eps_iso` argument of `point_el()` (default 1.e-3, the same as before).

## Verification

`post_point_search1` (found with 0.1, not with the default).
