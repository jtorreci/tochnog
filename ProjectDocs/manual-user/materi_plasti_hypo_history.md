# materi_plasti_hypo_history

## Description

`materi_plasti_hypo_history` (manual Professional 4.23) is an
initialization option that declares the **eight hypoplasticity history
variables** of the Professional, named `hyhis0` .. `hyhis7`. It is the
per-model name of the generic `materi_history_variables` for the
hypoplasticity laws, but with a FIXED number of variables (8) and with
the Professional basenames.

The eight variables (manual Professional 4.23) are:

| Slot | Name   | Meaning                                      |
|------|--------|----------------------------------------------|
| 0    | hyhis0 | void ratio e                                 |
| 1    | hyhis1 | substep size                                 |
| 2    | hyhis2 | mobilized friction angle                     |
| 3    | hyhis3 | stiffness measure                            |
| 4    | hyhis4 | structure s (sensitivity)                    |
| 5    | hyhis5 | overconsolidation ratio                      |
| 6    | hyhis6 | density index                                |
| 7    | hyhis7 | intergranular rho                            |

The variables share the same dof mechanism as the generic history
variables (`materi_history_variables`): one dof per variable, integrated
in the element loop. The only differences are the fixed count (8) and
the `hyhis` basename, so the dofs are addressable in
`target_item`/`post_point_dof`/`control_reset_dof` with the Professional
names (`-hyhis0`, `-hyhis4`, ...).

## Usage

In the initialization part, together with the mechanical dofs:

```
materi_velocity
materi_displacement
materi_stress
materi_strain_total
materi_plasti_hypo_history
end_initia
```

The initial void ratio is normally set with a reset:

```
control_reset_dof              2  -hyhis0
control_reset_value_constant   2  0.645
```

## Parameters

| Record | Meaning |
|--------|---------|
| `materi_plasti_hypo_history` | Declares the 8 hypoplasticity history variables `hyhis0`..`hyhis7` (no number follows the keyword). |

## Notes

- The kernel mapping (manual Professional 4.23): `hyhis0` is the void
  ratio `e` and `hyhis4` is the structure/sensitivity `s`. The Masin
  clay kernel (`masin_umat`) reads exactly these slots. The legacy
  `materi_history_variables` layout (e at `hisv6`, sensitivity at
  `hisv7`) is still honoured when the generic keyword is used instead.
- Requires `materi_stress` in the initialization part (checked by the
  group self-checks for the hypoplasticity group records).
- The corpus tests hypo1..hypo13 use this keyword; hypo2 and hypo4 pass
  rc=0 with the targets of the Professional.
