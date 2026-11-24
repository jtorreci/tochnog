# control_print_partialname

## Description

`control_print_partialname` is similar to the normal `control_print`
command, but instead of printing the given data items it prints EVERY
record whose name STARTS WITH one of the given data item names (prefix
match). For example `control_print_partialname 10 -element` prints all
records starting with `element` (as opposed to
`control_print 10 -element` which prints only the `element` record).

The output goes to stdout, exactly like `control_print`.

## Uso

```
control_print_partialname 20 -element
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. |
| `data_item_name_0 data_item_name_1 ...` | One or more data item names; all records whose registered name starts with any of them are printed (e.g. `-element` matches `element`, `element_group`, `element_mass`, `element_geometry`, ...). |

## Output

For every matched record, one line per active index with the record
name, the index and all its values — the same format as
`control_print`. The `control_print_filter` record (same index) applies.

## Example

```
control_print_partialname 20 -element
```

Prints every `element*` record of the model to stdout.

## Validation

Test `partialname` (validation-suite/test-2014): the captured stdout
contains the `element*` records (element, element_dof, element_group,
element_mass, ...) and NO `node*` records — the prefix discrimination
is checked with grep.
