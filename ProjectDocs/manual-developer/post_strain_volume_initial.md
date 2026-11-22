# post_strain_volume_initial

## Implementation (F5 T2 V9)

- **Record:** DOUBLE, length 1, indexed, class POST.
- **Fill:** `materi.cc`, `strain_volume_setup()` at the end of `input()`: the integrated initial volume of each strain_volume set (F5 T5 U8; before, `post.cc` copied ELEMENT_VOLUME once).

## Verification

`post_count_data1` (1 = Professional).

## Related

- [strain_volume_element](strain_volume_element.md): the volume strain and `post_strain_volume_absolute / _relative` (F5 T5 U8).
