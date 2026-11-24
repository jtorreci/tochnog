# post_node_factor

## Implementation (F5 T2 V9)

- **Record:** DOUBLE, length 1 (fixed), indexed, data_required post_node.
- **Consumer:** `post.cc`, `post()`, post_node loop: after the -average division and before the POST_NODE_RESULT PUT, `array_multiply` with the factor (default 1).

## Verification

`post_count_data1` (0.6 = Professional).
