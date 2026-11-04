# post_strain_volume_initial

## Implementation (F5 T2 V9)

- **Record:** DOUBLE, length 1, indexed, class POST.
- **Fill:** `post.cc`, `post_count_data_volume()`: for each `strain_volume_element` index without the record, the ELEMENT_VOLUME of its element(s), PUT once (kept afterwards).

## Verification

`post_count_data1` (1 = Professional).

## Pending

- `strain_volume_*` and `post_strain_volume_absolute / _relative` are registered but not computed (prescribed volume strain not implemented).
