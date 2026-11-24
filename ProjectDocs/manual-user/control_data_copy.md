# control_data_copy

## Description

Copy data items during the calculation (manual Professional 6.117/6.118).
`control_data_copy` copies ALL indices of `data_item_from` to
`data_item_to`, with an optional multiplication factor from
`control_data_copy_factor` (same index).

Normally the items should have the same length. As a special option you
can copy `node_inertia` to `node_force` records with a
`control_data_copy_factor` of -1: this substitutes material mass inertia
by static nodal forces for the remainder of the calculation (the
d'alembert principle).

The user is responsible for applying only logical copy actions
(compatible types; integer copies require factor 1).

## Usage

```
control_data_copy <index> <data_item_from> <data_item_to>
control_data_copy_factor <index> <factor>
```

## Example

```
control_data_copy          2  -group_materi_elasti_young  -group_materi_elasti_young
control_data_copy_factor   2  0.5
```

Every active index of the Young's modulus record is overwritten with its
own value times 0.5.

See also `control_data_copy_index` for a single-index copy.
