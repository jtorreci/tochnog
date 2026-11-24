# bounda_constant

## Description

Keeps the prescribed degrees of freedom (velocities, pressures, temperatures)
constant after the first step. On the first step the value of `bounda_time` is
applied; on the following steps the previous nodal value is kept, so the
prescribed dof does not vary with time.

Useful to prescribe a value that must stay fixed while the surrounding material
evolves (e.g. a fixed pressure or velocity that should not be re-applied from
the load curve every step). Requires `bounda_unknown` and `bounda_time` for the
same boundary (`iboun`).

## Usage

```
bounda_constant <iboun> <on/off>
```

## Parameters

| Parameter | Meaning                                                                  |
|-----------|--------------------------------------------------------------------------|
| `iboun`   | Boundary number of the prescribed-dof record to keep constant.           |
| `on/off`  | `-yes` to activate, `-no` to deactivate.                                 |

## Example

Keep the prescribed velocities of boundary 1 constant after the first step:

```
bounda_unknown  1  -geometry_point 1 -velx
bounda_time     1  0. 1. 100. 1.
bounda_constant 1  -yes
```

Boundary 1 applies the velocity from `bounda_time` on the first step and keeps
it unchanged on all subsequent steps.
