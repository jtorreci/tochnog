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

The void ratio is stored in `hisv[6]` (hyhis0). The dispatch requires
`materi_history_variables >= 8`, `materi_stress` and `materi_strain_total`;
`group_materi_memory -total` is rejected.

### Structure (sensitivity) extension

The optional `group_materi_plasti_hypo_masin_structure k A s_f` record adds
the meta-stable structure model of Masin (2007) *Canadian Geotechnical
Journal* 44(5):530-540 (manual reference [12], "Masin law" theory Eqs.
12-21). The current sensitivity `s` is the history variable `hyhis4`
(initialised through `control_reset_dof ... -hyhis4`, or `node_dof`); it must
be `>= 1`. The three parameters are:

| Parameter | Meaning |
|-----------|---------|
| `k` | Rate of structure degradation (`>= 0`). |
| `A` | Relative weight of volumetric vs shear damage strain, `0 <= A < 1`. |
| `s_f` | Final (fully destructured) sensitivity, `>= 1` (typically 1). |

The structure enters the model in two ways:

1. **Barotropy / shape** — the sensitivity factor
   `Si = (s - k(s - s_f))/s` multiplies the barotropy factor `fs`, divides the
   `9 r` denominator of `c1`, and replaces `kappa*` by `kappa* Si` in the
   exponent `alpha` (manual Eqs. 12, 18, 19).
2. **Pyknotropy** — the Hvorslev equivalent pressure `p*_e` is replaced by
   `s p*_e`, i.e. `fd = (2p/(s p*_e))^alpha` (Masin 2007 Eq. 8). The
   manual Eq. (12) prints the reference (`s = 1`) form.

The sensitivity evolves with the accumulated damage strain (manual Eq. 21):

```
sdot = -(k/lambda*) (s - s_f) sqrt( (eps_v_dot)^2 + (A/(1-A)) (eps_s_dot)^2 )
```

with `eps_v_dot = tr D` and `eps_s_dot = sqrt(2/3)||dev D||`. A basic model
without structure effects is recovered with `s = s_f = 1`; the kernel also
treats `s < 1` (e.g. no structure record, `hyhis4 = 0`) as inactive, so the
five-parameter model is bit-identical.

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

The structure variant is exercised by `hypo8.dat` (Pisa clay, `k=0.4`,
`A=0.1`, `s_f=1`, `e0=1.03`, `s0=2`): the final `sigyy` is `-230.534`
(Professional `-230.588`, target `-230.58 ±0.1`) with `s` degrading from
`2.0` to `1.609`.

**Limitations.** The intergranular-strain extension is not available for the
basic route; see the developer manual.

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
  used. For the basic route the `group_materi_plasti_hypo_masin_structure`
  record activates the Masin (2007) structure model; without it (or with
  `s < 1`) the five-parameter model is used.
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

### History variables (`materi_plasti_hypo_history`)

With the `materi_plasti_hypo_history` initialisation the eight slots are
(manual Professional 4.23): `hyhis0` void ratio e, `hyhis1` substep size,
`hyhis2` mobilized friction angle (degrees), `hyhis3` stiffness size
sqrt(Mijkl Mijkl), `hyhis4` structure s, `hyhis5` OCR, `hyhis6` density
index (von Wolffersdorff only, 0 here), `hyhis7` intergranular strain
rho. Only `hyhis0` and `hyhis4` are state variables (read back each
step); the others are for post-processing. For the basic law `hyhis1`
is the time step / 32 (fixed Runge-Kutta substeps), `hyhis5` = p_e/p
without the structure factor and `hyhis7` = 0; for the clay laws they
come from the kernel (adaptive substep, OCR with s, rho). Test
`validation-suite/test-2026/hypo_masin_history_slots1.dat`.

## Differences with the Professional version

- **hypoplasticity history slots and stress path** (known difference, not
  decided; AUDIT K14). Since F4 T6a (K14a) the GNU fills the
  post-processing slots like the Professional 02-08-2026: `hyhis2`, `hyhis3`
  and `hyhis5` of `hypo7`/`hypo8` and of the clay law `hypo12`/`hypo13`
  agree to about 3e-4 (hypo7: 23.4037 / 22516.4 / 2.09388, Pro 23.4076 /
  22518.8 / 2.09449). Remaining: `hyhis1` is the GNU substep (basic law:
  time step / 32, clay law: the adaptive substep) where the Professional
  reports the time step; `hyhis3` with intergranular strain near rho = 1
  (`hypo12` 157632, Pro 74297: the GNU reports the m_R L tangent of the
  kernel); `hypo9` (basic law with `group_materi_plasti_hypo_strain_
  intergranular`, not used by the GNU basic kernel: `hyhis3`/`hyhis7`
  differ); stress path `hypo3` -0.1034 (Pro -0.1629) and `hypo10` disx
  4.3e-9 (Pro 9.84e-8) unchanged. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
