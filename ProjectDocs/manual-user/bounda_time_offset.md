# bounda_time_offset

## Description

Sets the initial time offset for the load-only mode of `bounda_time`
(`bounda_time_increment`). The evaluation times start at `offset` instead of 0:
the load applied at time `t` is `bounda_time[k]` with
`k = floor( (t - offset) / increment )`, so `t = offset` corresponds to the
first stored load.

Requires `bounda_time_increment` for the same boundary (`iboun`).

## Usage

```
bounda_time_offset <iboun> <offset>
```

## Parameters

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `iboun`   | Boundary number whose load-only time axis is shifted.         |
| `offset`  | Initial time of the first load, must be >= 0.                 |

## Example

Start the load-only series of boundary 1 at `t = 0.02`:

```
bounda_time          1  1. 2. 3.
bounda_time_increment 1  0.01
bounda_time_offset   1  0.02
```

The first load `1.` is applied at `t=0.02`, `2.` at `t=0.03` and `3.` at
`t=0.04`. Before `t=0.02` the load is zero.
