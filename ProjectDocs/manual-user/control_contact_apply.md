# control_contact_apply

## Description

Per-timestep on/off switch for the contact algorithm. If the switch is
`-no`, the contact algorithm is not used; with `-yes` (default) it is.
The switch applies to the `control_timestep` record with the same index
(manual Professional 6.112). This is the per-timestep counterpart of the
global `contact_apply` record; any `-no` (global or per-timestep)
disables the algorithm.

Typical use: enable contact only in the timesteps where contact is
expected (e.g. after a placement stage), keeping the algorithm out of
the way during free phases.

## Usage

```
control_contact_apply <index> <-yes|-no>
```

## Parameters

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `index`   | Index of the `control_timestep` record the switch applies to.  |
| switch    | `-yes` use the contact algorithm (default), `-no` skip it.     |

## Example

```
control_timestep 0  0.01 0.1
control_contact_apply 0  -no
```

A block falling onto a contact face is NOT stopped when the gate is
`-no`: it follows the exact free-fall trajectory (disy = -0.055 after
t=0.1 with g=10), while with the contact active the same node stays
around -0.003 (test `contact_block` vs `contact_ctrl_apply`).
