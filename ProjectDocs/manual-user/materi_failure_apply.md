# materi_failure_apply

## Description

`materi_failure_apply switch` (manual Professional 6.801) is the global
(no index) form of `control_materi_failure_apply index switch` (6.143).
With `-no`, any failure data in the input file (`group_materi_failure_*`:
rupture, crunching, damage, plasti_kappa, void_fraction) is ignored: no
element is marked as failed and no element is deleted because of a
failure criterion.

## Usage

```
materi_failure_apply switch
control_materi_failure_apply index switch
```

| Parameter | Meaning |
|-----------|---------|
| `switch` | `-yes` (default) applies the failure criteria, `-no` ignores them |
| `index` | (control form) the `control_timestep` index the switch applies to |

## Precedence

The `control_materi_failure_apply` record of the current timestep overrules
the global `materi_failure_apply`; default `-yes`. Elements that already
failed in an earlier timestep keep their deletion schedule when a later
timestep switches failure off (only the detection of new failures is
disabled).

## Example

```
group_materi_failure_rupture 0 0.005 0.2
materi_failure_apply -no
```

Self-checking inputs: `validation-suite/test-2026/materi_failure_apply1.dat`
(global `-no`: `sigxx` 13.4615, identical to the model without the failure
record; with failure active the GNU gives 8.07692) and
`materi_failure_apply2.dat` (global `-no` overruled by
`control_materi_failure_apply 0 -yes`, GNU only).

## Differences with the Professional version

- The failure switches no longer switch off damage. Up to this change the
  GNU `control_materi_failure_apply -no` also froze `group_materi_damage_*`;
  on the Professional 02-08-2026 `materi_failure_apply -no` leaves the
  damage evolution unchanged, and the GNU now does the same.
- The Professional 02-08-2026 binary accepts `control_materi_failure_apply`
  only WITHOUT index (`control_materi_failure_apply 0 -no` is a parse
  error "I do not know : -no"), although its manual documents an index.
  The GNU follows the manual (indexed), so an unindexed
  `control_materi_failure_apply -no` is not accepted by the GNU; use
  `materi_failure_apply -no` for the global switch.
- The Professional 02-08-2026 stops with "Program error detected for data
  item: element_rhside_delete" as soon as an element fails in the test
  model, so only the `-no` behaviour (identical to the model without
  failure data) could be compared.
