# materi_plasti_tension_apply

## Description

`materi_plasti_tension_apply switch` (manual Professional 6.804, global, no
index): with `-no`, any tension-plasticity data in the input file is ignored
for all timesteps. The per-timestep form is
`control_materi_plasti_tension_apply index switch` (6.150), which overrules
the global record for the timesteps of its index. Useful to check a model
without its tension cut-off without commenting out the records.

## Usage

```
materi_plasti_tension_apply switch
control_materi_plasti_tension_apply index switch
```

| Parameter | Meaning |
|-----------|---------|
| `switch` | `-yes` (default) applies the tension plasticity, `-no` ignores it |
| `index` | (control form) the `control_timestep` index the switch applies to |

## What is ignored

- `group_materi_plasti_tension` (the tension yield surface);
- `group_materi_plasti_tension_direct` (the direct tension cut-off; the
  Mohr-Coulomb direct cut-off keeps its own default tension limit 0).

Plastic strains of earlier timesteps are kept.

## Precedence

The control record with the index of the current timestep overrules the
global record; default `-yes` (measured on the Professional 02-08-2026).

## Example

```
group_materi_plasti_tension_direct 0 1.
materi_plasti_tension_apply -no
```

Self-checking inputs in `validation-suite/test-2026/`:
`materi_plasti_tension_apply1.dat` (tension_direct ignored, sigyy
21.978), `materi_plasti_tension_apply2.dat` (global -no, control 1 -yes:
yield in timestep 1), `control_materi_plasti_tension_apply_yield1.dat` (the
control form also ignores the tension yield surface). All equal the
Professional 02-08-2026.

## Differences with the Professional version

None measured. Since F5 the control form also ignores
`group_materi_plasti_tension` (before it only ignored `_tension_direct`),
as the Professional and the manual ("any tension-plasticity data").
