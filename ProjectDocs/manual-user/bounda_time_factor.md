# bounda_time_factor

## Description

Multiplication factor for the loads specified by a `bounda_time` record
with the same index (manual Professional 6.36). Only the LOAD values are
multiplied; the time points are untouched. Default 1.

Handy when an imported time-load table uses another load definition
(e.g. accelerations as fractions of g).

## Usage

```
bounda_time_factor <index> <factor>
```

## Example

```
bounda_force 1  -geometry_line 1 -vely
bounda_time  1  0.0 5.0  100.0 5.0
bounda_time_factor 1  2.0
```

The effective load at any time is 10.0 (test `bt_factor`: sigyy 20 with
the factor, 10 without).

See also `bounda_time_units` (converts BOTH time and data units).
