# post_data

## Description

`post_data index dataitem_name_0 dataitem_index_0 dataitem_number_0 ...` (manual Professional 6.924). At every step each listed entry (data item, its index, the position in the record or a dof label such as `-velx`) is multiplied with its factor of `post_data_factor` and summed into `post_data_result` with the same index. Without `post_data_factor` the factors are 1.

## Usage

```
post_data 0 -node 2 0 -node 3 0 -node_dof 4 -velx
post_data_factor 0 1. 2. 10.
```

## Example

`validation-suite/test-2026/post_count_data1.dat`: three unit quad4 elements, two element groups, the right edge pulled. `post_data_result 0` = 1*x(node 2) + 2*x(node 3) + 10*velx(node 4) = 6.

## Differences with the Professional version

None measured. The manual does not say what happens without `post_data_factor`; the GNU uses 1.
