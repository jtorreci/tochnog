# bounda_time_increment

## Description

Declares that the values of `bounda_time` are load-only (a plain list of
amplitudes) instead of time/load pairs. The evaluation times are generated
uniformly as `offset + k*increment`, where `k` is the integer step index and
`offset` defaults to 0 (see `bounda_time_offset`).

The load applied at time `t` is `bounda_time[k]` with
`k = floor( (t - offset) / increment )`, clamped to the list range. Requires
`bounda_time` for the same boundary (`iboun`).

## Usage

```
bounda_time_increment <iboun> <increment>
```

## Parameters

| Parameter   | Meaning                                                           |
|-------------|-------------------------------------------------------------------|
| `iboun`     | Boundary number whose `bounda_time` values are load-only.        |
| `increment` | Fixed time step between successive loads, must be > 0.           |

## Example

Apply the loads `1.`, `2.`, `3.` every `0.01` time units:

```
bounda_time          1  1. 2. 3.
bounda_time_increment 1  0.01
```

At `t=0.` (k=0) the load is `1.`, at `t=0.01` (k=1) it is `2.` and at `t=0.02`
(k=2) it is `3.`. The load stays at `3.` for any later time (k is clamped to the
last index).
