# group_materi_elasti_k0

## Description

`group_materi_elasti_k0` (manual Professional 6.650) provides the K0
elastic data: when the record is specified AND
`control_materi_elasti_k0` is set to `-yes`, the K0 parameter is used in
the elastic stress law with `group_materi_elasti_young` or
`group_materi_elasti_young_power` (and `group_materi_elasti_poisson`), or
with `group_materi_elasti_hardsoil`, to determine the Poisson coefficient
consistent with K0:

    nu = K0 / (1 + K0)          (derived from K0 = nu/(1-nu))

For `K0 > 0.95` Tochnog takes `0.95`. The manual recommends the
combination as a convenient way to get "K0 stresses" when imposing gravity
in a geotechnics calculation: after gravity is imposed, drop the
`control_materi_elasti_k0` record so the normal
`group_materi_elasti_poisson` is used in the remaining steps.

## Usage

```
group_type 0  -materi
group_materi_elasti_young 0  1000.0
group_materi_elasti_poisson 0  0.3
group_materi_elasti_k0 0  0.5
group_materi_memory 0  -updated_without_rotation

control_timestep 1  0.1 1.0
control_materi_elasti_k0 1  -yes
```

With `K0 = 0.5`: `nu = 0.5/1.5 = 1/3`. The switch acts on the timesteps
whose `control_timestep` has the SAME index (here 1), like every other
`control_*` record.

## Parameters

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_elasti_k0` | `K0` | Coefficient of earth pressure at rest; `nu = K0/(1+K0)`, truncated at `K0 = 0.95`. |
| `control_materi_elasti_k0` | `switch` | `-yes` activates the K0 override of the Poisson ratio for the timesteps with the same control index. |

## Physics

The override applies when the control is `-yes` for the current control
index and `group_materi_elasti_k0` exists:

- `group_materi_elasti_young` / `group_materi_elasti_young_power`: the K0
  Poisson ratio replaces `group_materi_elasti_poisson`; it is also applied
  when no `group_materi_elasti_poisson` record is given (measured on the
  Professional).
- `group_materi_elasti_hardsoil`: the K0 Poisson ratio replaces both
  `nu50` and `nuur` of the hardsoil record.

`group_materi_elasti_young_polynomial` is not affected (not listed in the
manual). In a laterally confined column or oedometer the resulting stress
ratio is exactly `sigma_lateral/sigma_axial = nu/(1-nu) = K0`.

## Validation

- `mk0.dat` / `mk0_off.dat` (validation-suite/test-2014): confined
  compression, `E = 1000`, base `poisson 0.3`, `K0 = 0.5`, control index 0:
  `-yes` gives `sigma_xx = -0.75`, `sigma_yy = -1.5` EXACT; `-no` gives
  `-0.5769` / `-1.3462`.
- `elasti_k0_index1.dat` (test-2026, control index 1; Professional test
  `elasti7` uses index 20): gravity column, `K0 = 0.1`,
  `sigma_xx = -0.05`, `sigma_yy = -0.5` (Professional
  `-0.04999999884` / `-0.5`).
- `elasti_k0_nopoisson1.dat` (test-2026): young_power without any poisson
  record, same targets.
- `elasti_k0_hardsoil1.dat` (test-2026): hardsoil with `nu50 = nuur = 0.2`,
  same targets (without K0 the ratio would be 0.25).

## Notes

- Requires `materi_stress` and `materi_velocity` in the initialization
  part (`group_materi_elasti_young_power` additionally requires
  `materi_strain_total` in the GNU input check; the Professional does not).
- Up to 2026-09-29 the GNU only honored `control_materi_elasti_k0` with
  index 0 (it was silently ignored for any other control index).
