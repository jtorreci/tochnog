# post_node_rhside_ratio_method

## Implementation (F5 T2 V9)

- **Record:** INTEGER, length 1, no index; only `-post_node_rhside_free` is accepted (`db_error` otherwise).
- **Consumer:** `post.cc`, `post_node_rhside_fixed_free()`: with the method the ratio of a dof group is its summed free value (tmp1), without the fixed threshold; the largest group is kept as before. The ratio is also read by the iteration control (`top.cc`), so the method changes when an iteration counts as converged, as the manual intends.

## Verification

`post_node_rhside_ratio_method1` (ratio = free velocity sum, check in `scripts/build_safe.sh`).

## Pending

- BACKLOG C110: fixed / free are largest absolute values in the GNU, sums in the manual and the Professional.
