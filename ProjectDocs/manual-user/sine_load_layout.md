# sine time-load layout (`*_sine` family)

## Description

The `*_sine` records impose a time-dependent factor on a companion load
record (an edge force, a normal flux, a volume heat source, ...) as the
sum of sine functions:

```
<keyword> <index> <start_time> <end_time> <freq_0> <amp_0> [<freq_1> <amp_1> ...]
```

The factor at time `t` is `sum_i amp_i * sin(2*pi*freq_i*t)` and it is
applied only while `start_time < t <= end_time`; outside the window the
factor is 0 (the load is not applied). The layout is the same as
`bounda_sine` (manual Professional 6.34).

## Keywords covered

| Keyword | Professional manual | Companion load |
|---|---|---|
| `force_edge_sine` | 6.465 | `force_edge` |
| `force_edge_normal_sine` | 6.476 | `force_edge_normal` |
| `force_edge_projected_sine` | 6.487 | `force_edge_projected` |
| `groundflow_flux_edge_normal_sine` | 6.567 | `groundflow_flux_edge_normal` |
| `condif_heat_edge_normal_sine` | 6.81 | `condif_heat_edge_normal` |
| `condif_heat_volume_sine` | 6.88 | `condif_heat_volume` |

`force_volume_sine` (6.507) is DIFFERENT: its Professional signature has
NO `end_time` (`index start_time freq_0 amp_0 ...`) and is therefore
correct as-is; do not add an `end_time` to it.

## Example

A 10 Hz sine of amplitude 1e-3 on an edge force, active from `t=0` to
`t=1e8`:

```
force_edge_normal 20 5.0
force_edge_normal_geometry 20 -geometry_line 1
force_edge_normal_sine 20 0.0 1.e8 1.e1 1.e-3
```

## Differences with the Professional version

- The GNU previously read the record as `start_time freq amp ...` (no
  `end_time`), shifting every `(freq, amp)` pair by one slot and ignoring
  the window end. This is fixed: the GNU now uses the Professional
  layout and the `end_time` gate.
- The Professional special case "a frequency of 0 imposes a constant
  static value equal to the amplitude" is NOT implemented for these
  consumers (PENDING). A frequency of 0 contributes 0.
- A record written with the old 3-value GNU form now yields
  `nfreq = (len-2)/2` (e.g. length 3 -> 0) and contributes nothing. No
  GNU input in the repository used the old form, so no compatibility
  shim was added.

## See also

- `bounda_sine`, `force_edge`, `force_edge_normal`, `force_edge_projected`,
  `groundflow_flux_edge_normal`, `condif_heat_edge_normal`,
  `condif_heat_volume`
