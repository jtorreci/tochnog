# group_materi_elasti_k0

## Description

`group_materi_elasti_k0` (manual Professional 6.650) provides the K0
elastic data: when the record is specified AND
`control_materi_elasti_k0` is set to `-yes`, the K0 parameter is used in
the elastic stress law **with `group_materi_elasti_young` or
`group_materi_elasti_young_power` and `group_materi_elasti_poisson`** to
determine the Poisson coefficient consistent with K0:

    nu = K0 / (1 + K0)          (derived from K0 = nu/(1-nu))

For `K0 > 0.95` Tochnog takes `0.95`. The manual recommends the
combination as a convenient way to get "K0 stresses" when imposing gravity
in a geotechnics calculation: after gravity is imposed, drop the
`control_materi_elasti_k0` record so the normal
`group_materi_elasti_poisson` is used in the remaining steps.

## Uso

```
group_type 0  -materi
group_materi_elasti_young 0  1000.0
group_materi_elasti_poisson 0  0.3
group_materi_elasti_k0 0  0.5
group_materi_memory 0  -updated_without_rotation

control_materi_elasti_k0 0  -yes
```

With `K0 = 0.5`: `nu = 0.5/1.5 = 1/3`.

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_elasti_k0` | `K0` | Coefficient of earth pressure at rest; `nu = K0/(1+K0)`, truncated at `K0 = 0.95`. |
| `control_materi_elasti_k0` | `switch` | `-yes` activates the K0 override of the Poisson ratio. |

## Física

The override applies ONLY when all of the following hold: the control is
`-yes`, `group_materi_elasti_k0` exists, AND `group_materi_elasti_young`
(or `_young_power`) AND `group_materi_elasti_poisson` records exist (the
manual's combination; the `group_materi_elasti_hardsoil` combination is
NOT implemented — see the developer manual). In a laterally confined
oedometer the resulting stress ratio is exactly
`sigma_lateral/sigma_axial = nu/(1-nu) = K0`.

## Validation

- `mk0.dat` / `mk0_off.dat` (validation-suite/test-2014): confined
  compression (oedometer), `E = 1000`, `eps_zz = -0.001`, base
  `poisson 0.3`, `K0 = 0.5`.
  - with `control_materi_elasti_k0 0 -yes`: `nu = 1/3`,
    `sigma_xx = -0.75` and `sigma_yy = -1.5` EXACT
    (`sigma_xx/sigma_yy = 0.5 = K0`);
  - with `control_materi_elasti_k0 0 -no`: the base poisson `0.3` is
    used, `sigma_xx = -0.5769` and `sigma_yy = -1.3462` — A/B
    discrimination.

## Notas

- Requires `materi_stress` and `materi_velocity` in the initialization
  part.
- The `control_materi_elasti_k0` keyword was registered in Sprint 9 as a
  partial (no hook); the hook is now wired in `set_stress()`.
