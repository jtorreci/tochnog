# control_check_data

## Description

`control_check_data index switch` (manual Professional 6.111): with
`-yes` the database is checked at the start of control index `index`
(the check of [check_data](check_data.md): every specified data item
whose required item is missing is reported).

## Usage

```
control_check_data index switch
```

Test: `validation-suite/test-2026/control_check_data1` (one check, at
index 1).

## Differences with the Professional version

The GNU prints the result (`check_data: ...`); the Professional
02-08-2026 prints nothing when the database is consistent (DIVERGENCES
`OUTPUT-TEXT`, class C).
