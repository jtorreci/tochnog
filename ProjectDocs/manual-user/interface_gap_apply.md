# interface_gap_apply

## Description

`interface_gap_apply switch` (manual Professional 6.796, global): with
`-no`, any `group_interface_gap` record is ignored and every interface uses
the default gap 0 (an interface pulled apart opens). Default `-yes`.

## Usage

```
interface_gap_apply switch
```

| Parameter | Meaning |
|-----------|---------|
| `switch` | `-yes` (default) applies the gap records, `-no` ignores them |

## Example

```
group_interface_gap 10 1.e20
interface_gap_apply -no
```

Self-checking inputs in `validation-suite/test-2026/`:
`interface_gap_apply1.dat` (gap 1.e20 ignored: the pulled interfaces open,
stress 0) and `interface_gap_apply2.dat` (`-yes`: normal 10, shear 2).

## Differences with the Professional version

- The Professional 02-08-2026 rejects the record ("I don't know what to do
  with : interface_gap_apply"); the GNU implements the manual
  (DIVERGENCES IFACE-GAP-APPLY). The `control_interface_gap_apply` record
  the manual mentions (it overrules this one) has no manual entry and is
  not implemented.
