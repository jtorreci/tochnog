# group_materi_plasti_hypo_wolfersdorff (and variants)

## Description

Hypoplasticity for granular materials (the Karlsruhe school, von
Wolffersdorff / Gudehus). A rate-type constitutive model for sands and
gravels with a single nonlinear tensorial relation between stress rate and
strain rate, without a yield surface. The state is the stress tensor and the
void ratio.

All variants share the kernel `hypo.c` (pure-C port of the original Fortran
`hypo.f`). They are selected by which keyword is active:

| Keyword | hypo_type | Purpose |
|---------|-----------|---------|
| `group_materi_plasti_hypo_wolfersdorff` | 0 | Standard von Wolffersdorff model (8 params). |
| `group_materi_plasti_hypo_lowangles` | 1 | Wolfersdorff variant with two extra exponents (10 params). |
| `group_materi_plasti_hypo_cohesion` | — | Cohesion parameter, subtracted from the normal stresses (stabilisation at free surfaces). |
| `group_materi_plasti_hypo_intergranularstrain` | — | Intergranular strain extension (Niemunis-Herle) for small-strain stiffness. |
| `group_materi_plasti_hypo_pressuredependentvoidratio` | — | Switch to initialise the void ratio from the pressure (instead of a fixed value). |

Requires `materi_stress`, `materi_strain_total`, a velocity formulation, and
`materi_history_variables 4` (void ratio + intergranular strain if used).

## Usage

```
materi_history_variables 4

group_materi_plasti_hypo_wolfersdorff <element_group>
                        phi_c h_s n e_d0 e_c0 e_i0 alpha beta

group_materi_plasti_hypo_lowangles <element_group>
                        phi_c h_s n e_d0 e_c0 e_i0 alpha beta rval powxi

group_materi_plasti_hypo_cohesion <element_group>  c

materi_strain_intergranular
group_materi_plasti_hypo_intergranularstrain <element_group>
                        R m_R m_T beta_r chi

group_materi_plasti_hypo_pressuredependentvoidratio <element_group>  -yes|-no
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `phi_c` | Critical state friction angle (degrees). |
| `h_s` | Granular hardness (kPa). |
| `n` | Exponent of the barotropic factor. |
| `e_d0` | Minimum void ratio at zero pressure. |
| `e_c0` | Critical void ratio at zero pressure. |
| `e_i0` | Maximum void ratio at zero pressure. |
| `alpha` | Exponent of the pyknotropic factor. |
| `beta` | Exponent of the barotropic factor. |
| `rval` | Low-angles exponent (lowangles variant). |
| `powxi` | Low-angles exponent power (lowangles variant). |
| `c` | Cohesion (kPa); subtracted from the normal stresses before the law evaluation, and the linear contribution is used if the pressure falls below `-3*c`. |
| `R` | Intergranular strain radius. |
| `m_R` | Intergranular strain stiffness factor (loading). |
| `m_T` | Intergranular strain stiffness factor (unloading). |
| `beta_r` | Intergranular strain evolution exponent. |
| `chi` | Intergranular strain exponent. |

## Example

```
group_materi_plasti_hypo_wolfersdorff 0
                        30.   (deg)
                        5800.e3 (kPa)
                        0.28
                        0.84
                        0.53
                        1.00
                        0.13
                        1.05
```

This reproduces the regression test `hypo1.dat` (Karlsruhe sand, biaxial
compression, `sigyy=-863` at the end).

## Notes

- The variants can be combined: `cohesion` + `intergranularstrain` +
  `pressuredependentvoidratio` are additive records on top of the base model
  (wolfersdorff or lowangles).
- The kernel `hypo.c` is a pure-C port (validated by hypo1-4), reentrant
  after P4-F.
- Regression tests: `hypo1.dat` (wolfersdorff), `hypo2/3/4.dat`
  (intergranular strain), `hypo_cohesion.dat`, `hypo_lowangles.dat`,
  `hypo_pdvr.dat` (pressure-dependent void ratio).
