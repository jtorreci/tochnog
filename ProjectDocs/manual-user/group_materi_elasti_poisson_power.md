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
  The incremental code (the stiffness is evaluated with the previous-step
  stress) lands at `nu ~ 0.43`, `p ~ 2.06` -> measured
  `sigma_xx = -1.4389` and `sigma_yy = -3.2924` (84% of the fixed point;
  window targets that exclude the constant-poisson base
  `-0.5769/-1.3462`). `control_timestep_iterations 8` keeps the
  within-step coupling tight.

## Notas

- Requires `materi_stress` and `materi_velocity` in the initialization
  part (same check as the other elastic records).
- The base `group_materi_elasti_poisson` record is read only as a
  fallback/default and is overridden when the power record exists.
