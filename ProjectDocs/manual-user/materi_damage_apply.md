# materi_damage_apply

## Description

`materi_damage_apply switch` (manual Professional 6.798) is the global
(no index) form of `control_materi_damage_apply index switch` (6.140).
With `-no`, any damage data in the input file (`group_materi_damage_mazars`)
is ignored for all timesteps: the damage `dam` does not evolve and the
stresses are those of the model without damage records.

It is meant for testing an input file "just linear" without commenting out
every damage record.

## Usage

```
materi_damage_apply switch
control_materi_damage_apply index switch
```

| Parameter | Meaning |
|-----------|---------|
| `switch` | `-yes` (default) applies damage, `-no` ignores the damage data |
| `index` | (control form) the `control_timestep` index the switch applies to |

## Precedence

When both are given, the `control_materi_damage_apply` record with the index
of the current timestep overrules the global `materi_damage_apply` (measured
on the Professional 02-08-2026, same rule as `interface_gap_apply` /
`dependency_apply`). Without any record damage is applied.

With `-no` in a later timestep the damage reached so far is kept (frozen);
it is not reset. In such a timestep the damage data is ignored: the stress
increments are added without the (1 - d) factor (as in the Professional
02-08-2026; `materi_damage_apply2.dat`: `sigxx` 1.0872657685e7 in both).

## Example

```
group_materi_damage_mazars 0  2.6e-4 1. 15000. 1.2 649. 1.05
materi_damage_apply -no
control_materi_damage_apply 0 -yes   ( damage only in timestep 0 )
control_timestep 0 1.e-5 3.e-4
control_timestep 1 1.e-5 3.e-4
```

Self-checking inputs: `validation-suite/test-2026/materi_damage_apply1.dat`
(global `-no`: `dam` 0, `sigxx` 1.404e7, identical to the model without
damage) and `materi_damage_apply2.dat` (global `-no` overruled in timestep
0: `dam` 0.451188). Both pass on the Professional 02-08-2026.

## Differences with the Professional version

- The damage switches do not affect the element failure criteria and vice
  versa (same as the Professional).
- None known for the switch. (The damaged stress level differed before
  2026-09-29 because of a legacy defect, e.g. `sigxx` 1436.83 vs 85598.3 at
  the end of the `damage1` bar; fixed, see
  `materi_damage_effective_stress.md`.)
