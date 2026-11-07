# group_materi_elasti_borja_tamagnini

## Description

Borja Tamagnini nonlinear (hyper)elasticity for soils (manual Professional
6.641; Borja, R.I., Tamagnini, C. 1998, "Cam-Clay plasticity, part III",
Comput. Methods Appl. Mech. Engrg. 155, 73-95, eq. 3.4-3.12). The bulk
stiffness grows exponentially with the elastic volumetric compression and
the shear modulus grows with the pressure, as in Cam-Clay type soils. The
law can be combined with plasticity (it acts on the elastic strain).

With `eps_v` the elastic volumetric strain, `e` the deviatoric elastic
strain, `eps_s = sqrt(2/3 e:e)` and `Omega = -eps_v / khat`:

```
p_stress = -pr exp(Omega) (1 + 3 alpha eps_s^2 / (2 khat))   (mean stress)
mu       = G0 + alpha pr exp(Omega)
sigma    = p_stress I + 2 mu e
```

The stress at zero elastic strain is `-pr` (the reference pressure). As the
other hyperelastic laws, the stress changes by sigma(new strain) -
sigma(old strain), so initialise the stresses at `-pr` (for example with
`control_reset_dof` / `control_reset_value_constant`) to obtain the total
stress of the paper.

## Syntax

```
group_materi_elasti_borja_tamagnini index G0 alpha khat pr
```

| Parameter | Meaning |
|---|---|
| `G0` | shear modulus at zero pressure (mu_0) |
| `alpha` | coupling of the shear modulus with the pressure |
| `khat` | elastic compressibility (larger than 0) |
| `pr` | reference pressure (positive in compression) |

Requires `materi_stress` and `materi_strain_elasti`.

## Example

```
group_materi_elasti_borja_tamagnini 0 100. 2. 0.01 10.
control_reset_dof 0 -sigxx -sigyy -sigzz
control_reset_value_constant 0 -10.
```

Tests: `validation-suite/test-2026/group_materi_elasti_borja_tamagnini1`-`4`, `_error1`
(analytic targets).

## Differences with the Professional version

- The Professional 02-08-2026 implements a different, incremental
  "modified" law that does not follow the reference paper (DIVERGENCES
  BORJA-TAMAGNINI, class D, owner decision 2026-11-07: the paper is
  followed). Isotropic compression eps_v -0.03 from zero stress: GNU
  -190.855 for any step size, Professional -30 / -53.54 / -69.24 with 1 / 4
  / 20 steps. Pure shear: GNU modulus G0 + alpha pr, Professional
  G0 + alpha / khat.
