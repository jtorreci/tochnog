# group_materi_plasti_hypo_masin (and _clay)

## Description

Clay hypoplasticity models of Masin. The keyword family has two independent
routes with different kernels:

- **basic route** — `group_materi_plasti_hypo_masin` with the kernel
  `masin_basic.c`: the five-parameter Masin (2005) clay hypoplasticity model
  (Appendix C of Masin, D. (2005), *Géotechnique* 55(6):457-464).
- **_clay route** — `group_materi_plasti_hypo_masin_clay` with the kernel
  `masin.c`: the Masin (2014) HCEA formulation (*Géotechnique*
  64(3):232-238), a rate-independent hypoplastic model based on critical
  state soil mechanics with an explicit asymptotic state boundary surface,
  extended with stiffness anisotropy.

The rest of this document describes the `_clay` route unless stated
otherwise; the basic route has its own section below.

The law is formulated in kPa; the rest of the input file must be consistent
with those units.

The model uses **8 history variables** (`materi_history_variables 8`):
`hisv[0..5]` = intergranular strain tensor (Voigt), `hisv[6]` = void ratio
`e`, `hisv[7]` = sensitivity `s`. The initial void ratio is given through the
`node_dof` record (the 7th history component).

Requires `materi_stress`, `materi_strain_total`, a velocity formulation, and
`group_materi_memory -updated_without_rotation` (or another updated
formulation; `-total` is not supported for hypoplasticity).

## Basic route (Masin 2005)

The basic route implements the five-parameter Masin (2005) model from its
Appendix C rate equation in `masin_basic.c`. It is dispatched for
`group_materi_plasti_hypo_masin` when the group is neither `_clay` nor a
visco variant.

```
group_materi_plasti_hypo_masin <element_group>
                        phi_c lambda* kappa* N r
```

| Parameter | Meaning |
|-----------|---------|
| `phi_c` | Critical state friction angle, in degrees. |
| `lambda*` | Slope of the normal compression line (dimensionless). |
| `kappa*` | Slope of the swelling/loading line (dimensionless, `< lambda*`). |
| `N` | Position of the normal compression line (void ratio at 1 kPa). |
| `r` | Parameter controlling the asymptotic state boundary surface shape. |

The void ratio is stored in `hisv[6]`; the remaining history slots are not
used by the basic kernel. The dispatch still requires
`materi_history_variables >= 8`, `materi_stress` and `materi_strain_total`;
`group_materi_memory -total` is rejected.

**Limitations (partially implemented).** The structure evolution of Masin
(2005) Eq. (21) is NOT implemented: the
`group_materi_plasti_hypo_masin_structure k A s_f` record is parsed but has
no effect on the stress (the structure variable stays at its initial value).
The intergranular-strain extension is also not available for the basic
route. Both are pending work; see the developer manual.

**Example.** The corpus test `hypo7.dat` (drained axisymmetric triaxial,
`-updated_linear`) is the reference case:

```
group_materi_memory             0  -updated_linear
group_materi_plasti_hypo_masin  0
                        22.56    (phic deg)
                        0.11     (lambda*)
                        0.016    (kappa*)
                        1.375    (N)
                        0.4      (r)
```

With the reference-geometry strain increments for `-updated_linear` (R2),
the reported `eptyy` is `-0.300000000` and the final `sigyy` is `-231.789`
(Professional `-231.809`, target `±0.1`).

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

group_materi_plasti_hypo_masin_clay_visco <element_group>
                        Dr Iv
group_materi_plasti_hypo_masin_clay_visco_jm <element_group>
                        ocparam beta_deg ksi gama_deg Dref

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
| `ocparam` | Critical-state ratio parameter of the visco law (default 2.0). |
| `beta_deg` | Shear-band rotation angle, degrees (softening mechanism). |
| `ksi` | Visco exponent controlling the flow-direction weighting. |
| `gama_deg` | Flow-direction deviation angle, degrees (0 = auto). |
| `Dref` | Visco reference strain rate; activates the visco response (`> 0`). |
| `Dr` | Reference strain rate of the Niemunis visco law (e.g. 1e-6). |
| `Iv` | Viscosity index of the Niemunis visco law (e.g. 0.1). |

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

- The basic route (`group_materi_plasti_hypo_masin`) uses `masin_basic.c`
  (Masin 2005 Appendix C); the `_clay` route uses `masin.c` (HCEA). They do
  NOT share a kernel. The `_clay` variant additionally accepts the
  `_advanced_parameters`, `_avanced_direction`, `_clay_ocr` and
  `_clay_structure` records.
- The `_clay` structure record is optional; without it the clay model is
  used. For the basic route the structure record is accepted but currently
  ignored (see the limitations above).
- The OCR records require the matching `control_*_ocr_apply -yes`; they
  compute `e0 = exp(N - lambda* ln|OCR| - lambda* ln|p/pr|) - 1`.
- The intergranular-strain record activates the small-strain stiffness
  enhancement; `theta` is accepted for compatibility but has no direct slot
  in the kernel (chi is used).
- The visco extension comes in two formulations:
  - `_clay_visco Dr Iv` — the **Niemunis visco law** from the professional
    manual (creep rate `Dr*(1/OCR)^(1/Iv)`). The internal parameters
    (`ee0 pe0 lambda betaR`) are derived from the clay defaults. Requires
    `materi_history_variables 10`. Regression:
    `validation-suite/test-2014/hypomasin5.dat`.
  - `_clay_visco_jm ocparam beta_deg ksi gama_deg Dref` — the Jerman-Masin
    2020 formulation (rate scaling by `Dref`, shear-band rotation). The
    manual documents only `Dr Iv` for `_clay_visco`; the JM variant is kept
    under its own keyword. Regression: `hypomasin4.dat`.
  Both require `materi_history_variables 10`.
