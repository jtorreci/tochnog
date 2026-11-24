# print_gmsh_dummy

## Description

`print_gmsh_dummy switch` (manual Professional 6.992): default of
`control_print_gmsh_dummy` for every gmsh file (see
[control_print_gmsh](control_print_gmsh.md)): `-yes` (default) writes a
dummy point element (group 1234) in every node, `-no` does not. A
`control_print_gmsh_dummy` record for a control index overrules it.

## Usage

```
print_gmsh_dummy -no
```

| Parameter | Meaning |
|---|---|
| switch | `-yes` (default) / `-no`: dummy point elements. |

Tests `validation-suite/test-2026/print_gmsh_dummy1/2` (= Professional
02-08-2026: 1 and 5 elements).

## Differences with the Professional version

None for the switch (gmsh layout differences: BACKLOG C27).
