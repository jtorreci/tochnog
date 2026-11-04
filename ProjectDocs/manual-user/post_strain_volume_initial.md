# post_strain_volume_initial

## Description

Output record (manual Professional 6.966): the initial volume of the elements of `strain_volume_element` (same index). It is taken the first time the element volumes exist (the first step) and kept.

## Usage

```
target_item 5 -post_strain_volume_initial 10 0
```

## Example

`validation-suite/test-2026/post_count_data1.dat`: three unit quad4 elements, two element groups, the right edge pulled. (element 3, volume 1).

## Differences with the Professional version

Same value. PARTIAL: the prescribed volume strain itself (`strain_volume_absolute_time`, `post_strain_volume_absolute` / `_relative`) is not implemented in the GNU yet.
