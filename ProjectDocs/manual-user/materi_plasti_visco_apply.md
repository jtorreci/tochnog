# materi_plasti_visco_apply

## Description

`materi_plasti_visco_apply switch` (manual Professional 6.805, global, no
index): with `-no`, any visco-plasticity data in the input file is ignored
for all timesteps; the plasticity becomes rate independent. The
per-timestep form `control_materi_plasti_visco_apply index switch` (6.151)
overrules it for the timesteps of its index.

## Usage

```
materi_plasti_visco_apply switch
control_materi_plasti_visco_apply index switch
```

| Parameter | Meaning |
|-----------|---------|
| `switch` | `-yes` (default) applies the visco-plasticity, `-no` ignores it |
| `index` | (control form) the `control_timestep` index the switch applies to |

## Example

```
group_materi_plasti_tension 0 1.
group_materi_plasti_visco_power 0 1.e-1 1.
materi_plasti_visco_apply -no
```

Self-checking inputs in `validation-suite/test-2026/`:
`materi_plasti_visco_apply1.dat` (1D bar: sigxx 1 instead of 2.666667) and
`control_materi_plasti_visco_apply1.dat` (control 1 -no). Both equal the
Professional 02-08-2026.

## Differences with the Professional version

None measured.
