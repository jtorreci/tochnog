# post_quadrilateral_element_group

## Implementation (F5 T2 V9)

- **Record:** INTEGER, length 1 (fixed), indexed, data_required post_quadrilateral.
- **Consumer:** `post.cc`, `post()`, quadrilateral loop: sets `post_point_search_group` before its points; points not found are skipped and the weights renormalised (existing behaviour).

## Verification

`post_point_search1` (0.0880956 = Professional).
