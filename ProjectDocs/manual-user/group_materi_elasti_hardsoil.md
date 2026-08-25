# group_materi_elasti_hardsoil

## Description

`group_materi_elasti_hardsoil` (manual Professional 6.649) is the
elastic part of the **Hardening-Soil (HS) model** (theory section
"Hardening-Soil model"). It replaces the constant Young modulus with a
power law in the **minor principal stress** (the confining pressure):

```
first loading:      E50 = Eref_50 * ((sig3 + c*cot(phi))/(sigmaref_50 + c*cot(phi)))^m
                    nu  = nu50
unload/reload:      Eur = Eref_ur * ((sig3 + c*cot(phi))/(sigmaref_ur + c*cot(phi)))^m
                    nu  = nuur
```

The manual orders the principal stresses `sig3 > sig2 > sig1` with `sig1`
the **largest compressive** stress. In this code the sign convention is
tension positive, so the largest compression is the **smallest algebraic
eigenvalue** and the manual's `sig3` (the least compressive, the
confining stress in a triaxial test) is the **largest algebraic
eigenvalue** of the stress tensor:

```
sig3_manual = -(largest algebraic eigenvalue) + c*cot(phi)
```

Taking the smallest eigenvalue would give the AXIAL stress (the manual's
`sig1`) and the stiffness would depend on the axial load instead of the
confinement — that is not the HS model (the tests verify the correct
mapping analytically).

The cohesion terms `c*cot(phi)` are read from
[`group_materi_plasti_hardsoil`](group_materi_plasti_hardsoil.md); if
the plastic record is absent, `c = 0` is assumed (the cohesion term
vanishes identically, `cot(phi)` is never evaluated).

The loading/unloading switch uses the maximum |p| history dof
([`materi_plasti_hardsoil_history`](materi_plasti_hardsoil_history.md)):
if the current pressure of the step (estimated as `p_old + dp`, with
`dp` the elastic pressure increment of the step evaluated with the
first-loading tangent) is SMALLER than the maximum at the start of the
step, the material is unloading/reloading and `Eur`/`nuur` apply;
otherwise it is first loading and `E50`/`nu50` apply. The `dp` estimate
catches the first unloading step (same decision logic as
`group_materi_elasti_stress_pressure_history_factor`).

## Uso

```
group_materi_elasti_hardsoil 0
                        1000.    ( Eref_50 )
                        100.     ( sigmaref_50 )
                        0.3      ( nu50 )
                        0.5      ( m )
                        3000.    ( Eref_ur, typically 3*Eref_50 )
                        100.     ( sigmaref_ur )
                        0.2      ( nuur )
```

The data check requires the initia `materi_plasti_hardsoil_history`
(the maximum |p| history used by the loading/unloading switch).

## Parámetros

| # | Parameter | Meaning |
|---|-----------|---------|
| 1 | `Eref_50` | Reference secant stiffness at `sigmaref_50` (from the 50% strength triaxial test) |
| 2 | `sigmaref_50` | Reference confining stress for `Eref_50` |
| 3 | `nu50` | Poisson ratio for first loading (e.g. 0.3 drained) |
| 4 | `m` | Stress exponent (e.g. 0.5) |
| 5 | `Eref_ur` | Reference unloading/reloading stiffness (typically `3*Eref_50`) |
| 6 | `sigmaref_ur` | Reference confining stress for `Eref_ur` |
| 7 | `nuur` | Poisson ratio for unloading/reloading (e.g. 0.2 drained) |

The reference bases `sigmaref_50 + c*cot(phi)` and
`sigmaref_ur + c*cot(phi)` must be positive (data error otherwise).

## Validation

- `mhardsoil_elast` (plane-strain uniaxial, `c=10`, `phi=30 deg`,
  `m=0.5`): `E50 = 1000*sqrt(17.32/117.32) = 384.23` and
  `sigyy = -1.098901*384.23*0.0004 = -0.1689` EXACT (A/B
  `mhardsoil_elast2` with `sigmaref=200` gives `-0.1241`).
- `mhardsoil_unload` (confined oedometer, 4 load + 2 unload steps,
  `m=0`, `Eur = 3*E50`): `sigyy = -0.5385 + 1.346154*3000*0.0002 =
  +0.2692` EXACT (A/B `mhardsoil_unload_flat` with `Eur = E50` ends at
  `-0.2692`), the same numbers as `msph` with the stress factor.

**Clamp**: if the base `sig3 + c*cot(phi)` is not positive (a very
tensile `sig3` beyond the cohesion) the power law is not evaluated and
`E = Eref` is used (the manual does not specify the clamp; this avoids
zero or negative stiffness). With `m = 0` the modulus is `Eref` at all
stresses.
