# bounda_sine

## Description

Imposes the `bounda_dof` or `bounda_force` record with the same `index` as the
sum of one or more sine functions of time. Use it to prescribe a
displacement/force history with a known frequency content (for example, to
excite a model with a frequency spectrum and read the response with
`control_print_history`).

The record follows the Professional manual signature 6.34:

```
bounda_sine <index> <start_time> <end_time> <freq_0> <amp_0> [<freq_1> <amp_1> ...]
```

The prescribed value at time `t` is:

```
sum_i amp_i * sin(2*pi*freq_i*t)
```

and it is applied only while `start_time <= t <= end_time`. Before
`start_time` and after `end_time` the dof is no longer prescribed and becomes
free again (same convention as `bounda_time`: "at all times that a dof is not
prescribed in such way, it is free"). The sine functions start at `start_time`
with value 0.

## Parameters

| Parameter    | Meaning                                                                 |
|--------------|-------------------------------------------------------------------------|
| `index`      | Boundary number; must match a `bounda_dof` or `bounda_force` record.    |
| `start_time` | Time at which the sine starts to be imposed.                            |
| `end_time`   | Time after which the sine is no longer imposed (the dof becomes free).  |
| `freq_i`     | Frequency of the i-th sine term [1/time].                               |
| `amp_i`      | Amplitude of the i-th sine term.                                        |

At least one `(freq, amp)` pair is required.

## Example

Prescribe a 10 Hz sine of amplitude 1e-3 between `t=0` and `t=1e8` on
boundary 20:

```
bounda_dof   20  -bottom -disx
bounda_sine  20  0.0 1.e8 1.e1 1.e-3
```

## Differences with the Professional version

- The GNU version previously read the record as `start_time freq amp ...`
  (with no `end_time`), which shifted every `(freq, amp)` pair by one slot.
  This is fixed: the GNU now uses the Professional layout.
- The Professional special case "a frequency of 0 imposes a constant static
  value equal to the amplitude" is NOT implemented yet (PENDING). In the GNU a
  frequency of 0 produces a zero contribution.
- A record written with the old 3-value GNU form (`start freq amp`) now has
  length 3 and is interpreted as `nfreq = (3-2)/2 = 0`, i.e. it contributes
  nothing. No GNU input in the repository used the old form, so no
  compatibility shim was added.

## See also

- `bounda_dof`, `bounda_force`, `bounda_time`, `control_print_history`
