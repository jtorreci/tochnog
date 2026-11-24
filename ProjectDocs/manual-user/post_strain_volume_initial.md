# post_strain_volume_initial

## Description

Output record (manual Professional 6.966): the initial volume of the strain_volume set with the same index (see [strain_volume_element](strain_volume_element.md)). It is computed before the calculation from the input coordinates.

## Usage

```
target_item 5 -post_strain_volume_initial 10 0
```

## Example

`validation-suite/test-2026/post_count_data1.dat`: three unit quad4 elements, two element groups, the right edge pulled. (element 3, volume 1).

## Differences with the Professional version

Same value. `post_strain_volume_absolute` / `_relative` and the volume strain itself are on the [strain_volume_element](strain_volume_element.md) page (F5 T5 U8).
