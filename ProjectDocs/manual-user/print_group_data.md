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

Registered (INTEGER, no_index, one item name accepted) and parsed so
the input files of the corpus (distri3) load. The GiD/element-record
writing is PENDING (the corpus test that uses the record has no target
on it).
