# group_materi_plasti_hypo_masin

## Description

Clay hypoplasticity model of Masin (Masin 2014, *Géotechnique* 64(3):232-238)
for fine-grained soils. It is a rate-independent hypoplastic model based on
critical state soil mechanics with an explicit asymptotic state boundary
surface, extended with stiffness anisotropy.

The law is formulated in kPa; the rest of the input file must be consistent
with those units.

The model uses **8 history variables** (`materi_history_variables 8`):
`hisv[0..5]` = intergranular strain tensor (Voigt, currently inert in the basic
model), `hisv[6]` = void ratio `e`, `hisv[7]` = sensitivity `s`. The initial
void ratio is given through the `node_dof` record (the 7th history component).

Requires `materi_stress`, `materi_strain_total`, a velocity formulation, and
`group_materi_memory -updated_without_rotation` (or another updated
formulation; `-total` is not supported for hypoplasticity).

## Usage

```
materi_history_variables 8

group_materi_plasti_hypo_masin <element_group>
                        phi_c lambda* kappa* N r

group_materi_plasti_hypo_masin_structure <element_group>  k A s_f
group_materi_plasti_hypo_masin_ocr <element_group>  OCR
control_materi_plasti_hypo_masin_ocr_apply <element_group> -yes|-no
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `phi_c` | Critical state friction angle, in degrees. |
| `lambda*` | Slope of the normal compression line (dimensionless). |
| `kappa*` | Slope of the swelling/loading line (dimensionless, `< lambda*`). |
| `N` | Position of the normal compression line (void ratio at 1 kPa). |
| `r` | Parameter controlling the asymptotic state boundary surface shape (equivalent to `nu_pp` of the clay law). |
| `k` | Structure parameter (structure record). |
| `A` | Structure parameter, `0 <= A < 1` (structure record). |
| `s_f` | Final sensitivity, `>= 1` (structure record). |
| `OCR` | Overconsolidation ratio; initial void ratio is derived from it. |

Defaults applied by the implementation when the optional records are absent:
`p_t = 0` (no cohesion shift), `alpha_G = 1` (isotropic), `s_f = 1` (no
structure effect), `A_g = 0` (intergranular strain off), vertical direction `z`.

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
`p0=100 kPa`.

## Notes

- The structure record is optional; without it the basic model is used
  (`s_f=1`, `k=A=0`).
- The OCR record requires `control_materi_plasti_hypo_masin_ocr_apply -yes`;
  it computes `e0 = exp(N - lambda* ln|OCR| - lambda* ln|p/pr|) - 1`.
- The intergranular-strain history slots are reserved but not active in this
  first integration (the basic Masin law).
