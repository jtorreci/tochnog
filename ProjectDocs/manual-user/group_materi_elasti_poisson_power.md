# group_materi_elasti_poisson_power

## Description

`group_materi_elasti_poisson_power` (manual Professional 6.653) makes the
Poisson ratio a **power law of the pressure state** (theory section 2.2.2
of the Professional manual):

    nu = nu0 + nu1 * (p/p1)^alpha,      with the condition  nu <= nu2

where `p` is the pressure, `p = -(sig11 + sig22 + sig33)/3` (positive in
compression, same sign convention as `group_materi_elasti_young_power`).
It combines with `group_materi_elasti_young` or
`group_materi_elasti_young_power` (the young modulus is taken from those
records) and replaces the constant `group_materi_elasti_poisson`.

## Uso

```
group_type 0  -materi
group_materi_elasti_young 0  1000.0
group_materi_elasti_poisson 0  0.3
group_materi_elasti_poisson_power 0
                        0.2     ( nu0 )
                        0.1     ( nu1 )
                        0.5     ( nu2: nu capped at nu2 )
                        1.      ( p1, must be > 0 )
                        1.      ( alpha )
```

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_elasti_poisson_power` | `nu0 nu1 nu2 p1 alpha` | Power-law Poisson ratio: `nu = nu0 + nu1*(p/p1)^alpha` with the cap `nu <= nu2`. `p1` must be positive. |

## Física

The pressure-dependent Poisson ratio follows the same structure as the
`group_materi_elasti_young_power` law (E = E0 + E1(p/p1)^alpha) and is
evaluated with the CURRENT stress state each time the elastic stiffness is
computed. In a confined compression (oedometer) the lateral/axial stress
ratio is exactly `nu/(1-nu)` at every state, so the power law directly
controls the K0-like lateral stress growth with pressure.

## Validation

- `mpower.dat` (validation-suite/test-2014): oedometer with `E = 1000`,
  `eps_zz = -0.0012`, `nu0 = 0.2`, `nu1 = 0.1`, `nu2 = 0.5`, `p1 = 1`,
  `alpha = 1`. The self-consistent analytic fixed point
  (`nu = 0.2 + 0.1*p`, `p = E*eps/(3(1-2*nu))` -> `E*eps = 1.2`) is
  `nu = 0.4`, `p = 2`, `sigma_yy = -2.5714`, `sigma_xx = -1.7143`.
  The code is incremental (each step uses the stiffness of the previous
  step's stress, starting from `nu = nu0` at zero stress), so the state
  stays below the fixed point: `sigma_xx = -0.4389`,
  `sigma_yy = -1.4101` (Professional 02-08-2026: -0.4417 / -1.4123).
  `control_timestep_iterations 8`.
- `elasti_poisson_power_replace1.dat` (validation-suite/test-2026): the
  same model in the Professional dialect, targets -0.4389 / -1.4101.

**Change 2026-10-06 (audit K06).** Before this date the power law did not
replace the constant Poisson ratio: the stiffness was the SUM of the
elastic matrix of `young`/`poisson` and the one of `young`/`nu(p)`, about
twice the Professional stiffness (`mpower` gave -1.4389 / -3.2924; its
targets were recalibrated). Results of inputs with this record change
accordingly.

## Notas

- Requires `materi_stress` and `materi_velocity` in the initialization
  part (same check as the other elastic records).
- The base `group_materi_elasti_poisson` record is overridden when the
  power record exists (the elastic matrix is rebuilt with `nu(p)` only).
- With `group_materi_elasti_young_power` the modulus is the power-law
  modulus of that record.
