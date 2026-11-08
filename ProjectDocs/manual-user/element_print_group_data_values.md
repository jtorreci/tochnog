# print_group_data + element_print_group_data_values

## Description

`print_group_data -group_name_0 -group_name_1 ...` (manual Professional
6.991) lists group records. For every element the print record
`element_print_group_data_values <element> v0 v1 ...` (6.446) then holds
the values of those records for the element group of the element, in the
listed order (a record with more than one value contributes all of
them). A group without the record contributes one 0. Real records are
the values the element uses: the `dependency_item` diagrams (see
[dependency_item](dependency_item.md); `-time_current` at the end of the
step that starts) and the element distributions are applied (F5 T4 U3).

## Usage

```
print_group_data -group_materi_elasti_young -group_materi_elasti_poisson
```

## Example

`print_group_data1.dat` (validation-suite/test-2026) gives
`element_print_group_data_values 1 1000 0.3`.

## Differences with the Professional version

- The Professional 02-08-2026 writes the record as
  `element_print_group_data_value` (singular). The GNU writes the manual
  name and reads both spellings (DIVERGENCES PGD-VALUE-NAME, class D).
- With a `-time_current` dependency the Professional reports the value of
  the middle of the last step (BACKLOG C122); dof dependencies use the
  element node average in the GNU, the last integration point in the
  Professional.
- The manual also mentions the GiD plot of these values; the GNU fills only
  the database record.
