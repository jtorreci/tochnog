# control_data_copy_index

## Description

Copy one record index of a data item to another index (possibly of a
different item) during the calculation, with an optional multiplication
factor from `control_data_copy_index_factor` (same index) — manual
Professional 6.119/6.120. The user is responsible for applying only
logical copy actions (compatible types; integer copies require factor 1).

## Usage

```
control_data_copy_index <index> <data_item_from> <index_from> <data_item_to> <index_to>
control_data_copy_index_factor <index> <factor>
```

## Example

```
control_data_copy_index          1  -group_materi_elasti_young 0  -group_materi_elasti_young 1
control_data_copy_index_factor   1  2.0
```

The Young's modulus of group 1 becomes the one of group 0 times 2.0.

See also `control_data_copy` (all indices at once).
