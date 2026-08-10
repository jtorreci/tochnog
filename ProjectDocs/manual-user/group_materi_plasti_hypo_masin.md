# group_materi_plasti_hypo_masin (and _clay)

## Description

Clay hypoplasticity model of Masin (Masin 2014, *Géotechnique* 64(3):232-238)
for fine-grained soils. It is a rate-independent hypoplastic model based on
critical state soil mechanics with an explicit asymptotic state boundary
surface, extended with stiffness anisotropy.

The law is formulated in kPa; the rest of the input file must be consistent
with those units.

The model uses **8 history variables** (`materi_history_variables 8`):
`hisv[0..5]` = intergranular strain tensor (Voigt), `hisv[6]` = void ratio
`e`, `hisv[7]` = sensitivity `s`. The initial void ratio is given through the
`node_dof` record (the 7th history component).

Requires `materi_stress`, `materi_strain_total`, a velocity formulation, and
`group_materi_memory -updated_without_rotation` (or another updated
formulation; `-total` is not supported for hypoplasticity).

## Usage

```
materi_history_variables 8

group_materi_plasti_hypo_masin <element_group>
                        phi_c lambda* kappa* N r

group_materi_plasti_hypo_masin_clay <element_group>
                        phi_c lambda* kappa* N nu_pp

group_materi_plasti_hypo_masin_structure <element_group>  k A s_f
group_materi_plasti_hypo_masin_ocr <element_group>  OCR
control_materi_plasti_hypo_masin_ocr_apply <element_group> -yes|-no

group_materi_plasti_hypo_masin_clay_advanced_parameters <element_group>
                        alpha_G alpha_f ay oc
group_materi_plasti_hypo_masin_clay_avanced_direction <element_group>  diri
group_materi_plasti_hypo_masin_clay_ocr <element_group>  OCR
control_materi_plasti_hypo_masin_clay_ocr_apply <element_group> -yes|-no
group_materi_plasti_hypo_masin_clay_structure <element_group>  k A s_f

group_materi_plasti_hypo_strain_intergranular_masin_clay <element_group>
                        R A_g n_g m_rat beta_r chi [theta]
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `phi_c` | Critical state friction angle, in degrees. |
| `lambda*` | Slope of the normal compression line (dimensionless). |
| `kappa*` | Slope of the swelling/loading line (dimensionless, `< lambda*`). |
| `N` | Position of the normal compression line (void ratio at 1 kPa). |
| `r` | Parameter controlling the ASBS shape (basic law). |
| `nu_pp` | Poisson-like parameter of the clay law. |
| `alpha_G` | Anisotropy of stiffness (1 = isotropic, 2 = typical London clay). |
| `alpha_f` | Exponent for the peak-friction anisotropy (0 = auto). |
| `ay` | Shape parameter (default 0.30). |
| `oc` | Critical-state ratio parameter (default 2.0). |
| `diri` | Vertical direction: 0 in 1D, 1 in 2D, 2 in 3D (default). |
| `k` | Structure parameter. |
| `A` | Structure parameter, `0 <= A < 1`. |
| `s_f` | Final sensitivity, `>= 1`. |
| `OCR` | Overconsolidation ratio; initial void ratio is derived from it. |
| `R` | Intergranular strain radius (e.g. 5e-5). |
| `A_g` | Intergranular strain stiffness factor (e.g. 270). |
| `n_g` | Intergranular strain exponent (e.g. 1). |
| `m_rat` | Intergranular strain m_T/m_R ratio (e.g. 0.5). |
| `beta_r` | Intergranular strain beta (e.g. 0.08). |
| `chi` | Intergranular strain chi (e.g. 7). |

Defaults applied by the implementation when the optional records are absent:
`p_t = 0` (no cohesion shift), `alpha_G = 1` (isotropic),
`alpha_f = auto`, `ay = 0.30`, `oc = 2.0`, `s_f = 1` (no structure effect),
`A_g = 0` (intergranular strain off), vertical direction = z (3D).

## Example

```
materi_history_variables 8
...
node_dof -ra -from 1 -to 4 -ra
         0. 0.
         0. 0.
         -100. 0. 0. -100. 0. -100.  (kPa, initial stress)
         0. 0. 0. 0. 0. 0.
         0. 0. 0. 0. 0. 0. 0.7 0.   (hisv[6]=e0=0.7)
...
group_materi_memory  0  -updated_without_rotation
group_materi_plasti_hypo_masin 0
                        21.9
                        0.095
                        0.015
                        1.19
                        0.1
```

This reproduces the reference biaxial test
(`validation-suite/test-2014/hypomasin1.dat`) with London clay parameters:
`phi_c=21.9`, `lambda*=0.095`, `kappa*=0.015`, `N=1.19`, `r=0.1`, `e0=0.7`,
`p0=100 kPa`. The anisotropic variant (`hypomasin2.dat`, `alpha_G=2`) and the
intergranular-strain variant (`hypomasin3.dat`) are also regression-tested.

## Notes

- `group_materi_plasti_hypo_masin` (basic law, 5 params) and
  `group_materi_plasti_hypo_masin_clay` (anisotropic law, 5 params) share the
  same integration kernel; the clay variant additionally accepts the
  `_advanced_parameters`, `_avanced_direction`, `_clay_ocr` and
  `_clay_structure` records.
- The structure record is optional; without it the basic model is used.
- The OCR records require the matching `control_*_ocr_apply -yes`; they
  compute `e0 = exp(N - lambda* ln|OCR| - lambda* ln|p/pr|) - 1`.
- The intergranular-strain record activates the small-strain stiffness
  enhancement; `theta` is accepted for compatibility but has no direct slot
  in the kernel (chi is used).
- The visco extension (`Dr Iv`) is NOT yet exposed; see the developer manual.
