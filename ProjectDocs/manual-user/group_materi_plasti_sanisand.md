# group_materi_plasti_sanisand

## Description

SANISAND — "Simple ANIsotropic SAND" model of Dafalias & Manzari (2004,
*J. Engng. Mechanics* ASCE 130(6):622-634) for sands. It is a rate-independent
plasticity model with a conical yield surface that opens around the hydrostatic
axis, kinematic hardening of the back stress `alpha`, a fabric tensor `z`
that evolves during plastic loading, and memory of the back stress at stress
reversal points (`alpha_sr`).

Key features:
- Critical state line in the `e:p` plane (`e0`, `lambda`, `xi`) and in the
  `q:p` plane (`M_c`, `M_e`).
- Bounding surface and dilatancy surface (state parameter `psi`).
- Fabric-dilatancy coupling: the dilatancy coefficient `A_d` grows with the
  accumulated fabric.
- Handles drained and undrained (pore-water bulk modulus `bulk_w`) loading.

The model is formulated in the same pressure units as `p_a`; the rest of the
input file must be consistent.

Requires `materi_stress`, `materi_strain_total`, a velocity formulation, and
`group_materi_memory -updated_without_rotation`. It uses **36 history
variables** (`materi_history_variables 36`).

## Usage

```
materi_history_variables 36

group_materi_plasti_sanisand <element_group>
                        p_a e0 lambda xi M_c M_e mm G0 nu
                        h0 c_h n_b A0 n_d z_max c_z bulk_w ptmult e
```

## Parameters

| Param | Meaning |
|-------|---------|
| `p_a` | Atmospheric pressure (reference). |
| `e0` | Void ratio on the critical state line at p = 0. |
| `lambda` | CSL slope in the `e:p` plane. |
| `xi` | CSL exponent in the `e:p` plane. |
| `M_c` | CSL slope in `q:p`, triaxial compression (or `phi_c` in degrees if > 5). |
| `M_e` | CSL slope in `q:p`, triaxial extension (or `phi_e` in degrees if > 5; 0 implies `phi_e = phi_c`). |
| `mm` | Opening of the yield surface cone. |
| `G0` | Shear modulus constant. |
| `nu` | Poisson's ratio. |
| `h0` | Plastic modulus constant. |
| `c_h` | Plastic modulus constant. |
| `n_b` | Plastic modulus constant (bounding surface). |
| `A0` | Dilatancy constant. |
| `n_d` | Dilatancy constant. |
| `z_max` | Maximum fabric index. |
| `c_z` | Fabric evolution constant. |
| `bulk_w` | Pore water bulk modulus (undrained); 0 for drained/consolidation. |
| `ptmult` | Mean-stress shift `p_t = ptmult * p_a` (stabilisation). |
| `e` | Initial void ratio. |

Typical values (Nevada sand, from the reference): `p_a=101, e0=0.834,
lambda=0.019, xi=0.7, M_c=1.25, M_e=0, mm=0.01, G0=125, nu=0.05, h0=7.05,
c_h=0.968, n_b=1.1, A0=0.704, n_d=2.1, z_max=4.0, c_z=600, bulk_w=0,
ptmult=1.0, e=0.70`.

## Example

```
materi_history_variables 36
...
node_dof -ra -from 1 -to 4 -ra
         0. 0.
         0. 0.
         -100. 0. 0. -100. 0. -100.  (kPa, initial stress)
         0. 0. 0. 0. 0. 0.
         0. 0. 0. 0. 0. 0. 0.7 0. ...  (36 history values; hisv[6]=e0)
...
group_materi_memory  0  -updated_without_rotation
group_materi_plasti_sanisand 0
                        101.0 0.834 0.019 0.7 1.25 0.0 0.01 125.0 0.05
                        7.05 0.968 1.1 0.704 2.1 4.0 600.0 0.0 1.0 0.70
```

This reproduces the regression test
(`validation-suite/test-2014/hyposanisand1.dat`): Nevada sand, laterally
confined biaxial compression.

## Notes

- The back stress `alpha` is initialised automatically so that the yield
  surface cone is centred about the current stress state; the fabric `z` and
  the reversal memory `alpha_sr` start at zero.
- **Validation status**: the void ratio matches the reference Fortran UMAT
  to 0.001% and the early triaxial steps (1-8) to <0.12%. The late
  deviator stress (step 20) differs by ~6% — this is the reproducibility
  limit between two numerically equivalent implementations (the same C code
  recompiled with different optimisation flags varies MORE than the
  C-vs-Fortran difference), NOT a model error. See the developer manual for
  the full analysis.
