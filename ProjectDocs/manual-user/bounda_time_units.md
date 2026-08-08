# bounda_time_units

## Description

Converts the units of the values stored in `bounda_time`. All times in the
series are multiplied by `factor_time` and all loads/lengths by
`factor_length`. Useful when the input data is written in other units than
the ones used by the model (e.g. hours instead of seconds).

Requires `bounda_time` for the same boundary (`iboun`).

## Usage

```
bounda_time_units <iboun> <factor_time> <factor_length>
```

## Parameters

| Parameter      | Meaning                                       |
|----------------|-----------------------------------------------|
| `iboun`        | Boundary number whose `bounda_time` is scaled. |
| `factor_time`  | Multiplier applied to the time values.        |
| `factor_length`| Multiplier applied to the load/length values. |

## Example

Scale the time axis of boundary 1 from hours to seconds:

```
bounda_time        1  1. 1. 2. 1. 3. 1.
bounda_time_units  1  3600. 1.
```

The time values `1.`, `2.`, `3.` are interpreted as `3600.`, `7200.` and
`10800.` seconds. The load values (`1.` in each pair) are unchanged.
