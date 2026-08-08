# bounda_found

## Description

Diagnostic flag that reports whether a `bounda_*` record was actually applied to
any node of the model. `-yes` means the record found at least one node to act
on; `-no` means it was never used.

It does not affect the calculation. Use it to debug boundary conditions that
might not apply to any node (e.g. a wrong boundary number, a mismatched
geometry or a time range that never matches). The result is written to the
database, so it can also be inspected after the run.

## Usage

```
bounda_found <iboun> <yes/no>
```

## Parameters

| Parameter | Meaning                                                             |
|-----------|---------------------------------------------------------------------|
| `iboun`   | Boundary number of the `bounda_*` record to report on.              |
| `yes/no`  | `-yes` to print the flag, `-no` to suppress the report (default).   |

## Example

Report whether boundary 1 was actually applied to any node:

```
bounda_unknown  1  -geometry_point 1 -velx
bounda_time     1  0. 1. 100. 1.
bounda_found    1  -yes
```

`-yes` prints that boundary 1 matched at least one node; `-no` prints that it
was never applied, meaning the prescribed velocity had no effect.

**Status:** partially implemented. The record is accepted in the input and
stored, but the flag value itself is not written to the database yet, so no
report is printed. Only the full implementation produces a diagnostic output.
