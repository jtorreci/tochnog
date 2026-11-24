# group_materi_plasti_cap1

## Description

`group_materi_plasti_cap1` (manual Professional 6.691) is the first cap
model of the Professional version for **permanent plastic deformations
under high pressures in granular materials**. It is intended to be used
**in combination with shear plasticity models** (Drucker-Prager,
Mohr-Coulomb, ...): the standard driver evaluates all plastic laws at
each integration point and activates the one with the **largest yield
function value**, so the cap limits the mean pressure p while the shear
law governs the deviatoric strength.

The cap surface in the (p, q) plane is

```
f = q^2/M^2 + p*(p* - p*c) = 0
p*  = p + c*cot(phi)
p*c = pc + c*cot(phi)
```

with `p = -(sig11+sig22+sig33)/3` (positive in compression),
`q = sqrt(0.5*[(s11-s22)^2 + (s22-s33)^2 + (s33-s11)^2] + 3*(s12^2 +
s23^2 + s31^2))`, `phi` the Coulomb friction angle (radians), `c` the
cohesion and `pc` the history variable of the model (an ellipsoid/tear
shape in the p-q plane that translates along the p axis as pc hardens).
The flow rule is associative (`f_flow = f_yield`); the flow direction is
provided by the standard driver via central finite differences.

## Uso

In the data part, on the material group, together with the shear
plasticity laws (optional):

```
group_materi_plasti_cap1 0
                        0.5    ( phi, radians )
                        10.    ( c )
                        1.2    ( M, tangent of the Critical State Line )
                        0.2    ( lambda_star, compression index )
                        0.02   ( kappa_star, swelling index )
                        833.3  ( K_ref, bulk modulus at p_ref )
                        100.   ( p_ref )
                        0.     ( m, exponent of the hardening law )
```

The history variable `pc` comes from the initialization option
[`materi_plasti_cap1_history`](materi_plasti_cap1_history.md) and its
initial value is given via `node_dof` (see that page for the exact
record layout). Without the initia the data check fails.

## Parámetros

| # | Parameter | Meaning |
|---|-----------|---------|
| 1 | `phi` | Coulomb friction angle (radians) |
| 2 | `c` | Cohesion |
| 3 | `M` | Tangent of the Critical State Line (typically `M = 6*sin(phi)/(3 - sin(phi))`, but it is READ from the record, not computed) |
| 4 | `lambda_star` | Compression index (e.g. 0.15) |
| 5 | `kappa_star` | Swelling index (e.g. 0.03) |
| 6 | `K_ref` | Bulk modulus at the stress `p_ref` (typically `E_ref/(3(1-2*nu))`) |
| 7 | `p_ref` | Reference pressure of the hardening law |
| 8 | `m` | Exponent of the hardening law (e.g. 0.6; m = 0 gives linear hardening) |

Hardening (manual theory cap1, rate form):

```
eps_p_cv_dot = (lambda_star/kappa_star - 1)/K_ref * (p_ref/p*c)^m * pc_dot
```

i.e. the cap plastic volumetric strain rate drives the history variable
`pc`; `pc` only grows while the cap produces volumetric plastic
compression and stays constant during unloading/reloading.

## Validation

- `mcap1.dat` (registered in `scripts/build_safe.sh`): single hex8,
  isotropic compression, `phi = 0.5 rad`, `c = 10`, `M = 1.2`,
  `lambda* = 0.2`, `kappa* = 0.02`, `K_ref = 833.3`, `p_ref = 100`,
  `m = 0`, `pc_0 = 100`. The cap activates at `p = pc = 100`, hardens to
  `pc = 110.0` after 20 plastic steps (discrete fixed point
  `dpc = (5.0 - dpc)/9 -> 0.5/step`) and the 5-step elastic unload ends
  at `sigxx = -84.9996` (analytic `-85.0`, error 0.0005%).
- `mcap1_elast.dat`: elastic twin without the cap -> `sigxx = -175.0000`
  (exact) — discriminates surface AND hardening.
- `mcap1_comb.dat`: cap1 + druck_prag (phi 30 deg, c 10). On the
  isotropic path `f_dp < 0` always, so cap1 dominates (max-f selection)
  and the response is identical to `mcap1`.

## Notas

- `phi` must be in the open interval (0, pi/2): `cot(phi)` diverges at
  `phi = 0` (the manual expects `phi > 0`); both `sin(phi) = 0` and
  `cos(phi) = 0` raise a data error.
- `M != 0` is required.
- The hardening update uses the **trace of the plastic strain increment
  of the step** (clamped to >= 0, i.e. only compression hardens). When
  combined with shear plasticity, the volumetric plastic strain of the
  shear mechanism also contributes to the trace — for an exact cap-only
  hardening the cap should be the dominant mechanism in compression
  (see the developer manual for the design decision).
