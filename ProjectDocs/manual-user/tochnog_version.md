# tochnog_version

## Description

`tochnog_version` carries the BUILD DATE of the binary as a queryable
record: day, month, year (manual Professional 6.1091). The GNU writes
it at the start of the calculation, parsed from the compiler's
`__DATE__`, so the record always matches the running binary. It can be
targeted with `target_item` like any record.

## Usage

The record is written by the program; users read it (typically from
the `.dbs` dump or through a target):

```
tochnog_version  <day> <month> <year>
```

Target example (value 3 is the year):

```
target_item 10  -tochnog_version 0 2
target_value 10  2026. 0.5
```

## Example

The `tslv_ver` test of the validation suite targets the build year and
passes with a 2026 binary.
