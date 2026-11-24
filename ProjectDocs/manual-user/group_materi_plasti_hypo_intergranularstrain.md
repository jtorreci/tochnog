# group_materi_plasti_hypo_intergranularstrain (variable 5-6 parameters)

## Description

`group_materi_plasti_hypo_intergranularstrain` (manual Professional
6.715) sets the intergranular strain parameters of the hypoplastic law.
The record is **variable-length**: the Professional parser accepts a
short list and stops at the next keyword, so both of these forms are
valid:

```
group_materi_plasti_hypo_intergranularstrain <index>
  R  m_R  m_T  beta_r  chi  theta
group_materi_plasti_hypo_intergranularstrain <index>
  R  m_R  m_T  beta_r  chi
```

| Param | Meaning |
|-------|---------|
| R | radius of the intergranular strain surface |
| m_R | stiffness factor on reversal |
| m_T | stiffness factor on transverse loading |
| beta_r | interpolation exponent |
| chi | exponent of the rho^chi factor of the stiffness |
| theta | exponent of the rho^theta f_d N S_hat term (optional) |

When `theta` is omitted (the sfnet corpus form, 5 values) the GNU kernel
uses `chi` for the `rho^theta f_d N S_hat` term, which the manual
recommends for monotonic loading (`theta = chi`).

## Usage

Six values (Professional manual form):

```
group_materi_plasti_hypo_strain_intergranular 0
                       1.e-4 5.0 2.0 0.50 6.0 6.0
                       (R m_R m_T beta_r chi theta)
```

Five values (sfnet corpus form):

```
group_materi_plasti_hypo_intergranularstrain 0
                       1.e-4 5.0 2.0 0.50 6.0
                       (R m_R m_T beta_r chi)
```

(The Professional spelling of the record is
`group_materi_plasti_hypo_strain_intergranular`, accepted as an alias of
the GNU canonical name `group_materi_plasti_hypo_intergranularstrain`.)

## Notes

- Requires `materi_strain_intergranular` in the initialization part and
  `materi_plasti_hypo_history` (or `materi_history_variables`).
- `hypo2`/`hypo3`/`hypo4` of the sfnet corpus use the 5-value form and
  pass their targets; the validation suite uses the 6-value form.
- `R` must be positive.

## Differences with the Professional version

- **Strain reversal with intergranular strain (`hypo3`)** (deliberate
  divergence, class A, project decision 2026-10-12; DIVERGENCES K14b).
  Oedometric von Wolffersdorff + intergranular strain test with a strain
  reversal (test-2014 `hypo3`). Converged reference (independent Python
  integration of the same law, manual 2.2.4): sigxx -0.10466. GNU: -0.1034
  at dt 1e-3 (1.2 % off), -0.104235 at dt 1e-5. Professional 02-08-2026:
  -0.1629 at dt 1e-3 (56 % off), -0.10499 at dt 1e-5: a time-discretisation
  error at the reversal that vanishes only as dt -> 0. Guard tests
  `test-2026/hypo_wolff_igs_oedometer1/2`. The printed history slots
  (`hyhis*`) of `hypo3` still differ (AUDIT K14, open). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
