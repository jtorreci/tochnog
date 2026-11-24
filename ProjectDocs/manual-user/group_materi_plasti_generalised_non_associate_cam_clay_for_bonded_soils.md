# group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils

Generalised non-associated Cam clay for bonded soils and soft rocks (manual Professional 6.700, with
the initialisations 4.21 and 4.38).

## Description

The Professional manual names the model, "a modification of the 'Milan' model of Prof. Roberto Nova",
but gives no parameters (6.700 reads `index ...`), and the Professional 02-08-2026 stops on every input
we tried. The GNU therefore defines the model from the published Nova bonded-soil model:

- R. Nova, "Elastoplastic models for soils and soft rocks", 16th ICSMGE, 2005, pp. 3651-3656 [A].
- M. Arroyo, R. Castellanza, R. Nova, "Compaction bands and oedometric testing in cemented soils",
  Soils and Foundations 45(2), 2005, pp. 181-194 [B].

Stresses are compression positive here. With the mean stress p and the deviatoric stress q
(q = sqrt(3 J2)):

- **Bonding.** Two history variables with stress dimension: p_s, the stress memory of the unbonded
  soil (like the Cam clay preconsolidation pressure), and p_m, the bonding. The tensile strength is
  p_t = k p_m.
- **Shifted stress.** p* = p + p_t, eta = q/p*, and the size of the surfaces is p_c = p_s + p_m + p_t.
  The bonded surface therefore runs from p = -p_t to p = p_s + p_m on the p axis.
- **Yield surface and plastic potential.** Both have the Lagioia, Puzrin & Potts (1996) shape, each
  with its own (mu, alpha, M). The shape integrates the stress-dilatancy rule
  d = mu (M - eta) (alpha M/eta + 1):

  ln(p*/p_c) + I(eta) = 0,  I(eta) = int_0^eta e/Q(e) de,
  Q(e) = (1-mu) e^2 + mu (1-alpha) M e + mu alpha M^2.

  The yield surface uses (mu_f, alpha_f, M_f); the plastic potential uses (mu_g, alpha_g, M_g), so the
  flow is non-associated. Plastic volume change stops at eta = M_g (the critical state).
- **Hardening and debonding.** With the plastic volumetric strain eps_v^p (compression positive) and
  the plastic deviatoric strain eps_s^p:

  dp_s = rho_s p_s (deps_v^p + xi_s deps_s^p)
  dp_m = -rho_m p_m (|deps_v^p| + xi_m deps_s^p)

  Plastic compaction and shear (xi_s) harden p_s; any plastic volume change and shear (xi_m) destroy
  the bonding.
- **Elasticity.** The ordinary GNU elasticity records of the group (for example
  `group_materi_elasti_young` and `group_materi_elasti_poisson`).

## Syntax

```
group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils  index
   mu_f alpha_f M_f  mu_g alpha_g M_g  rho_s xi_s  rho_m xi_m  k
```

| Value | Meaning | Range |
|---|---|---|
| mu_f, alpha_f, M_f | Shape of the yield surface | mu_f > 0, alpha_f >= 0, M_f > 0 |
| mu_g, alpha_g, M_g | Shape of the plastic potential | mu_g > 0, alpha_g >= 0, M_g > 0 |
| rho_s, xi_s | Hardening of p_s by plastic compaction and by plastic shear | rho_s >= 0 |
| rho_m, xi_m | Debonding of p_m by plastic volume change and by plastic shear | rho_m >= 0 |
| k | Tensile strength p_t = k p_m | k >= 0 |

The record takes exactly 11 values. A negative or zero value outside these ranges stops with
`Error detected for data item`.

The record needs the initialisations `materi_stress`, `materi_strain_plasti` and
[materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils_history](materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils_history.md),
whose dofs `gencamhis0` (p_s) and `gencamhis1` (p_m) must be given initial values, for example with
`control_reset_dof`; p_s must be positive. The plastic strain of this model alone can be added with
[materi_strain_plasti_generalised_non_associate_cam_clay_for_bonded_soils](materi_strain_plasti_generalised_non_associate_cam_clay_for_bonded_soils.md).

## Example

```
materi_strain_plasti
materi_stress
materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils_history
end_initia
...
group_materi_elasti_young    0  10000.
group_materi_elasti_poisson  0  0.25
group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils 0
  0.9 0.4 1.1  0.9 0.4 1.2  10. 0.5  5. 1.  0.1

control_reset_dof             1 -sigxx -sigyy -sigzz
control_reset_value_constant  1 -100.
control_reset_dof             2 -gencamhis0
control_reset_value_constant  2 100.
control_reset_dof             3 -gencamhis1
control_reset_value_constant  3 50.
```

## Behaviour to know

- **Step size.** The GNU return mapping takes the flow direction at the elastic trial stress, so the
  flow rule is met to the first order of the step. In a constant-volume shear test the stress ratio at
  the critical state is 1.125, 1.186 and 1.195 for strain steps of 5e-4, 1e-4 and 5e-5 (M_g = 1.2).
  Use small steps.
- **Isotropic tension.** The element yields at p = -p_t with dilatant plastic flow, which softens both
  p_s and p_m.

## Verification

Tests in `validation-suite/test-2026`:

- `materi_plasti_gencam_bonded_iso1`: isotropic compression against the closed form
  p = p_s + p_m, p_s = p_s0 exp(rho_s eps_v^p), p_m = p_m0 exp(-rho_m eps_v^p).
- `materi_plasti_gencam_bonded_tension1`: isotropic extension, failure at the tensile strength
  k p_m.
- `materi_plasti_gencam_bonded_shear1`: constant-volume shear with shear hardening and debonding; the
  final state lies on the yield surface by an independent quadrature of I(eta).
- `materi_plasti_gencam_bonded_error1` to `_error5`: input errors.

The closed forms of I(eta) were checked against numerical quadrature on 156 parameter sets (largest
relative difference 5e-12).

## Differences with the Professional version

GENCAM-BONDED-DEFINITION: the Professional 02-08-2026 registers the three names but does not document
the parameters and never ran the model in our probes (fewer than 12 values: data error; 12 to 20
values: "Program error, assertion failed"). The GNU record, its 11 values and the behaviour above are
a GNU definition from the literature, not a convergence with the Professional.

## Future work

- **Lode-angle dependence (13 values).** Separate compression and extension slopes M_c and M_e for the
  yield surface and the potential, as in [B].
- **The full model of [B].** Also the four-parameter hyperelasticity of Borja, Tamagnini & Amorosi
  (1997), so the parameter sets of [B] Table 2 can be used as published.
