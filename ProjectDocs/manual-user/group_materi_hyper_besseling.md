# group_materi_hyper_besseling (and the hyperelastic family)

## Description

Hyperelastic (large-strain, total-Lagrange) rubber/foam material models. The
strain-energy function `W` is evaluated from the right Cauchy-Green tensor `C`
and the stress is obtained by central differences of `W` with respect to `C`,
so every model below is available with a fully consistent tangent
(`group_materi_hyper_stiffness -yes`, the default).

All models share the same code path and can be combined (the energy
contributions `W` of the active models are added). The following keywords are
implemented and verified:

| Keyword | Parameters | Strain energy contribution |
|---------|------------|----------------------------|
| `group_materi_hyper_besseling` | `K_1 K_2 alpha` | `K_1*(J1-3)^alpha + K_2*(J2-3)` |
| `group_materi_hyper_blatz_ko` | `G beta` | `G/2*(I1-3 + 2/beta*(J^-beta - 1))` |
| `group_materi_hyper_mooney_rivlin` | `K_1 K_2` | `K_1*(J1-3) + K_2*(J2-3)` |
| `group_materi_hyper_neohookean` | `K_1` | `K_1*(J1-3)` |
| `group_materi_hyper_reduced_polynomial` | `K_1 K_2 ...` | `sum K_i*(J1-3)^i` |
| `group_materi_hyper_volumetric_linear` | `K` | `K/2*(J-1)^2` |
| `group_materi_hyper_volumetric_murnaghan` | `K beta` | `K/beta*(1/(beta-1)*J^-beta+1)*J` |
| `group_materi_hyper_volumetric_ogden` | `K beta` | `K/beta*(1/beta*(J^-beta-1)+ln(J))` |
| `group_materi_hyper_volumetric_polynomial` | `K_0 K_1 ...` | `sum K_i/2*(J-1)^(2*(i+1))` |
| `group_materi_hyper_volumetric_simotaylor` | `K` | `K/2*((J-1)^2 + ln(J)^2)` |

Where `I1,I2,I3` are the invariants of `C`, `J = sqrt(I3)` is the volume
ratio, `J1 = I1/I3^(1/3)` and `J2 = I2/I3^(2/3)` are the isochoric invariants.
The energy is truncated at zero (`W<0 -> W=0`) to protect the tangent in
extreme compression.

Useful for rubber, foam, and other highly-compressible or nearly-incompressible
materials under large deformations (geotextiles, sealing layers, foundations
on soft rubber-like subgrades).

Requires `materi_stress`, `materi_strain_elasti`, `materi_strain_total`, and
`group_materi_memory -total` (or `-total_linear`) with `materi_displacement`.

## Usage

```
group_materi_hyper_besseling <element_group>  K_1 K_2 alpha
group_materi_hyper_blatz_ko    <element_group>  G beta
group_materi_hyper_mooney_rivlin <element_group>  K_1 K_2
group_materi_hyper_neohookean  <element_group>  K_1
group_materi_hyper_reduced_polynomial <element_group>  K_1 K_2 ...
group_materi_hyper_volumetric_linear <element_group>  K
group_materi_hyper_volumetric_murnaghan <element_group>  K beta
group_materi_hyper_volumetric_ogden <element_group>  K beta
group_materi_hyper_volumetric_polynomial <element_group>  K_0 K_1 ...
group_materi_hyper_volumetric_simotaylor <element_group>  K
group_materi_hyper_stiffness <element_group>  -yes|-no
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `K_1, K_2, K_i` | Hyperelastic stiffness coefficients of the chosen model. |
| `alpha` | Exponent of the Besseling isochoric term. |
| `G` | Linear shear modulus of the Blatz-Ko model. |
| `beta` | Blatz-Ko `beta = 2*nu/(1-2*nu)`; or volumetric law exponent (`beta != 0`, `beta != 1`). |
| `K` | Bulk modulus of the volumetric law. |
| `group_materi_hyper_stiffness` | `-yes` (default) computes the consistent tangent `Chyper`; `-no` uses the elastic tangent only. |

## Example

```
element_group 1
group_materi_elasti_young   0  3.0
group_materi_elasti_poisson 0  0.46
group_materi_hyper_blatz_ko 0  1.0  14.0
group_materi_memory 0 -total
```

This is the `blatz1.dat` regression test: a single axisymmetric quad4 element
stretched to a stretch ratio of 2 must reach a true stress of `3.18`
(`target_value 1 3.18021 0.01`), which the implementation reproduces.
