# group_materi_plasti_hypo_intergranularstrain (6 parameters)

## Description

`group_materi_plasti_hypo_intergranularstrain` (manual Professional
6.715) takes **six parameters**:

```
group_materi_plasti_hypo_intergranularstrain <index>
  R  m_R  m_T  beta_r  chi  theta
```

| Param | Meaning |
|-------|---------|
| R | radius of the intergranular strain surface |
| m_R | stiffness factor on reversal |
| m_T | stiffness factor on transverse loading |
| beta_r | interpolation exponent |
| chi | exponent of the rho^chi factor of the stiffness |
| theta | exponent of the rho^theta f_d N S_hat term |

The GNU previously accepted five values; the corpus tests (hypo2,
hypo10, triaxial_compression_drained_hypo_intergranular_isa) give all
six. The sixth value (`theta`) is the exponent of the
`rho^theta f_d N S_hat` term of the stiffness; for monotonic loading the
manual recommends `theta = chi`, and the GNU kernel evaluates that term
with `chi`, so both coincide for the corpus inputs.

## Usage

```
group_materi_plasti_hypo_strain_intergranular 0
                       1.e-4 5.0 2.0 0.50 6.0 6.0
                       (R m_R m_T beta_r chi theta)
```

(The Professional spelling of the record is
`group_materi_plasti_hypo_strain_intergranular`, accepted as an alias of
the GNU canonical name `group_materi_plasti_hypo_intergranularstrain`.)

## Notes

- Requires `materi_strain_intergranular` in the initialization part and
  `materi_plasti_hypo_history` (or `materi_history_variables`).
- hypo2/hypo3 use this record with 6 values and pass their targets.
