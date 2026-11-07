# group_materi_plasti_hypo_niemunis_visco

Also covers `group_materi_plasti_hypo_niemunis_visco_ocr` and
`control_materi_plasti_hypo_niemunis_visco_ocr_apply`.

## Description

Niemunis visco-hypoplasticity with intergranular strains (manual
Professional 6.720, 6.721, 6.148 and the theory section "Niemunis visco
law"; Niemunis, A. 2003, Extended hypoplasticity models for soils,
habilitation, section 4.2). It describes time dependent soil behaviour:
creep, relaxation and rate dependence, mainly for clays. The stress rate is

```
dT = M:D - L:Dvis
L  = fb [F^2 I + a^2 That That + b2 (I - 1/3 1 1)],  fb = -trT/((1 + a^2/3) kappa)
Dvis = Dr mhat (1/OCR)^(1/Iv),  OCR = pe/pe+
ln((1 + ee0)/(1 + e)) = lambda ln(pe/pe0)
```

with `That = T/trT`, `F` the Lode dependent von Wolffersdorff factor,
`a = sqrt3 (3 - sin phi)/(2 sqrt2 sin phi)`,
`b2 = (1 + a^2/3)(1 - 2 nu)/(1 + nu) - 1`, `mhat` the normalised flow rule
`m = -[F^2/a^2 (That + That*) + (That:That) That* - That (That*:That*)]`,
and `pe+` the Rendulic cap invariant with the flattening factor betaR
(`eta = q/(M p)`, `M = 6 F sin phi/(3 - sin phi)`). M is the stiffness
with intergranular strains (`group_materi_plasti_hypo_strain_intergranular`
R mR mT beta_r chi, as for the von Wolffersdorff law, without the N term).
The void ratio (history variable `hyhis0`) follows the volumetric strain.

## Syntax

```
group_materi_plasti_hypo_niemunis_visco index phi nu Dr Iv ee0 pe0 lambda betaR kappa
group_materi_plasti_hypo_niemunis_visco_ocr index OCR
control_materi_plasti_hypo_niemunis_visco_ocr_apply index -yes|-no
```

| Parameter | Meaning |
|---|---|
| `phi` | critical friction angle (degrees) |
| `nu` | Poisson ratio of the shear stiffness term b2 |
| `Dr` | reference creep rate |
| `Iv` | viscosity index |
| `ee0`, `pe0` | void ratio and pressure of a point on the reference isotach |
| `lambda` | compression index (Butterfield, ln(1+e) versus ln p) |
| `betaR` | flattening factor of the Rendulic cap (larger than 0) |
| `kappa` | swelling index (Butterfield, ln(1+e) versus ln p) |
| `OCR` | over-consolidation ratio at the start (`_ocr` record) |

Requires `materi_stress`, `materi_plasti_hypo_history` (or
`materi_history_variables`), `materi_strain_intergranular` and
`group_materi_plasti_hypo_strain_intergranular`. Give the initial void ratio
in `hyhis0` (`control_reset_dof`), or the OCR: with
`control_materi_plasti_hypo_niemunis_visco_ocr_apply -yes` (default `-no`)
the initial void ratio is set once, in the first timestep of that control
index, from `pe = OCR pe+` of the current stress:
`e0 = (1 + ee0) (pe/pe0)^-lambda - 1`.

History variables: `hyhis0` void ratio, `hyhis1` substep size, `hyhis2`
mobilised friction angle, `hyhis3` size of the stiffness M, `hyhis5` OCR.

Errors:

- `materi_strain_intergranular should be initialized.`
- `Error: phi and nu in group_materi_plasti_hypo_niemunis_visco give a negative b2 = (1+a^2/3)(1-2nu)/(1+nu)-1.`
  (the manual: the b2 equation only holds for a non-negative right-hand
  side; with phi 30 this needs nu <= 0.338).

## Example

From the Professional test `other/hypo6.dat` (relaxation, all boundaries
fixed):

```
group_materi_plasti_hypo_strain_intergranular 0 1.e-4 5.0 2.0 0.50 6.0 6.0
group_materi_plasti_hypo_niemunis_visco 0 30.0 0.2 100. 0.05 1.5 0.1 0.1 0.95 0.02
control_reset_dof             2  -hyhis0
control_reset_value_constant  2  1.48
```

Tests: `validation-suite/test-2026/group_materi_plasti_hypo_niemunis_visco1`-`4`,
`_ocr1`-`_ocr3`, `_error1`, `_error2`.

## Differences with the Professional version

- b2 < 0 is an error (manual); the Professional 02-08-2026 computes with the
  negative b2 (DIVERGENCES NIEMUNIS-VISCO-B2, class D, owner decision
  2026-11-07).
- The results agree with the Professional 02-08-2026 within about 5e-5
  relative (explicit integration details; for example hypo6 sigxx
  -0.0836179 against -0.0836218). The initial void ratio from the OCR is
  identical; in the first step of the index the Professional keeps that
  void ratio, the GNU already applies the strain of the step (1e-5 after
  1000 steps in `_ocr3`).
- `hyhis3` holds the size of M; the Professional writes another stiffness
  measure there (DIVERGENCES C120-NIEMUNIS-VISCO-HYHIS3, class B).
