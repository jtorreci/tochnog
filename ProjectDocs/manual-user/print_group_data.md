# print_group_data

## Description

`print_group_data` plots `group_*` data items (e.g.
`-group_materi_elasti_young`) in the GiD output for isoparametric
elements, and fills the `element_print_group_data` records (manual
Professional 6.990).

## Usage

```
print_group_data -group_materi_elasti_young ...
```

## Status

Implemented: see [element_print_group_data_values](element_print_group_data_values.md).
Real group records are printed with their `dependency_item` diagrams and
element distributions applied (F5 T4 U3). No GiD plot.
