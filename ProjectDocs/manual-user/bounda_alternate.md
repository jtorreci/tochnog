# bounda_alternate

## Description

Omits one of the listed `bounda_dof` indices per iteration, rotating through
the list between successive iterations. Only one of the listed boundaries is
skipped at a time; the next iteration skips the next one, cycling.

Useful for very large runs with limited memory: by alternating the applied
conditions on velocities/pressures, the system is solved in parts instead of
keeping every boundary active at once.

For example, if `bounda_dof` records exist with indices 10, 20 and 30, then
with `bounda_alternate 0 10 20 30` the successive iterations omit: 10, 20,
30, 10, ...

## Usage

```
bounda_alternate <index> <bounda_index_0> <bounda_index_1> ...
```

## Parameters

| Parameter        | Meaning                                                    |
|------------------|------------------------------------------------------------|
| `index`          | Record index; use `0`.                                     |
| `bounda_index_i` | Index of a `bounda_dof` record to alternately omit.        |

## Example

Alternate between omitting `bounda_dof` indices 0 and 10:

```
bounda_alternate 0 0 10
```

With `bounda_dof 0 ...` and `bounda_dof 10 ...` present, successive iterations
omit index 0, then 10, then 0, and so on.

**Status:** the record is accepted and the rotation logic is in place, but the
skip currently does not take effect — see the developer manual for details.
