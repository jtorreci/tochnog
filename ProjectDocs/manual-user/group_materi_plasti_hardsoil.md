# group_materi_plasti_hardsoil

## Description

`group_materi_plasti_hardsoil` (manual Professional 6.703) is the
plastic part of the **Hardening-Soil (HS) model**. The yield function
reads

```
f = q/(E50*(1 - q/qa)) - 2*q/Eur - gamma_p
```

with

- `q` the equivalent shear stress
  `q = sqrt(0.5*[(s11-s22)^2 + (s22-s33)^2 + (s33-s11)^2] +
  3*(s12^2 + s23^2 + s31^2))`;
- `E50`, `Eur` the power-law moduli of
  [`group_materi_elasti_hardsoil`](group_materi_elasti_hardsoil.md)
  evaluated at the CURRENT minor principal stress (coupled
  elasticity-plasticity, one step lag);
- `qa = qf/Rf` the asymptotic shear stress;
- `gamma_p` the equivalent plastic shear strain (hardening variable).

The manual gives `qa = qf/Rf` but not the explicit `qf`. The standard
Schanz derivation is used: in triaxial `q = sig1 - sig3` and at
Mohr-Coulomb failure
`sig1 = sig3*(1+sin(phi))/(1-sin(phi)) + 2c*cos(phi)/(1-sin(phi))`,
hence

```
qf = sig1 - sig3 = 2*sin(phi)*(sig3 + c*cot(phi))/(1 - sin(phi))
```

with `sig3` the minor principal stress in the manual (compression
positive) convention (`sig3_manual = -(largest algebraic eigenvalue) +
c*cot(phi)` in this tension-positive code, see the elastic page).

The flow rule is **associative**: the driver computes the flow direction
by central finite differences of `f` (the law only defines `f`). The
parameter `psi` (dilatancy angle) is accepted for input compatibility;
the non-associative HS flow potential is not implemented (see the
developer manual).

The hardening variable `gamma_p` is the accumulated plastic strain size
from [`materi_plasti_kappa`](materi_plasti_kappa.md)
(`kappa = int sqrt(0.5*deps_p:deps_p)`, the same measure the manual
calls "equivalent plastic shear strain") plus the optional extra initial
contribution of
[`control_materi_plasti_hardsoil_gammap_initial`](control_materi_plasti_hardsoil_gammap_initial.md).

The manual warns that the model "requires sufficient small timesteps";
it is suited for monotonic loading.

## Uso

```
group_materi_elasti_hardsoil 0
                        1000.    ( Eref_50 )
                        100.     ( sigmaref_50 )
                        0.3      ( nu50 )
                        0.       ( m )
                        3000.    ( Eref_ur )
                        100.     ( sigmaref_ur )
                        0.2      ( nuur )
group_materi_plasti_hardsoil 0
                        0.523599 ( phi, radians = 30 deg )
                        10.      ( c )
                        0.       ( psi, radians )
                        0.9      ( Rf )
```

The data check requires `materi_stress`, `materi_velocity`,
`materi_plasti_kappa` (the hardening variable), `materi_strain_plasti_hardsoil`
and `materi_plasti_hardsoil_history` to be initialized.

## Parámetros

| # | Parameter | Meaning |
|---|-----------|---------|
| 1 | `phi` | Maximum friction angle (radians, `0 < phi < pi/2`) |
| 2 | `c` | Cohesion |
| 3 | `psi` | Maximum dilatancy angle (radians, accepted; the flow is associative, see Description) |
| 4 | `Rf` | Failure ratio, `0 < Rf` (typically 0.9) |

## Validation

- `mhardsoil_gp0` (initial deviatoric stress `sigxx = -2` with the
  `gammap_initial` control): the yield function value at the initial
  state `f = q/(E50*(1-q/qa)) - 2*q/Eur = 2/(1000*(1-2/38.49)) -
  4/3000 = 0.0007763` is stored EXACTLY as the extra gamma_p: with the
  control the state does not relax (`sigxx = -2.0` EXACT, `kappa = 0`).
- `mhardsoil_gp0_off` (same, without the control): the return relaxes
  the deviatoric stress (`sigxx -2 -> -1.225`) and `kappa` grows to
  `0.000293`; the end state satisfies `f ~ 0` (it sits on the expanding
  surface): `0.687/(1000*(1-0.687/38.49)) - 2*0.687/3000 - 0.000293 =
  -0.000061`.
- `mhardsoil_plast` (initial `sigxx = -4`): `kappa = 0.0006` = 2x the
  `sigxx = -2` case: the hardening scales with the initial deviatoric
  stress.
