# group_materi_plasti_visco_exponential_name / _values / _limit

## Description

Exponential visco-plasticity data per plasticity model (manual Professional
6.744-6.746), the exponential twin of
`group_materi_plasti_visco_power_name` / `_value`. `_name` lists plasticity
models of the group; `_values` gives gamma and alpha of the exponential law
(as `group_materi_plasti_visco_exponential`) for each name, in the same
order.

When an integration point yields, the active plasticity model is looked up
in the names: a named model is evaluated elasto-viscoplastic with
`lambda = gamma * (-pressure) * exp(alpha*f)`, a model that is not named
stays elasto-plastic.

`group_materi_plasti_visco_exponential_limit` caps the exponent argument
`alpha*f` (default 3); it applies to the plain
`group_materi_plasti_visco_exponential` record and to the named data.

Precedence: `group_materi_plasti_visco_exponential`, then
`group_materi_plasti_visco_power`, then the power names, then the
exponential names; only the first one present is used.
`materi_plasti_visco_apply -no` / `control_materi_plasti_visco_apply -no`
also switch the named data off.

## Syntax

```
group_materi_plasti_visco_exponential_name index -name_0 -name_1 ...
group_materi_plasti_visco_exponential_values index gamma_0 alpha_0 gamma_1 alpha_1 ...
group_materi_plasti_visco_exponential_limit index limit
```

`group_materi_plasti_visco_exponential_value` (singular) is accepted as an
alias of `_values`. Each gamma and alpha must be larger than 0. Sub-surfaces
count for their record (`_druckprag` apex, `_matsuokanakai` apex, the
`_mohrcoul` corners).

## Example

```
group_materi_plasti_vonmises 0 1.
group_materi_plasti_visco_exponential_name 0 -group_materi_plasti_vonmises
group_materi_plasti_visco_exponential_values 0 0.5 0.1
group_materi_plasti_visco_exponential_limit 0 1.
```

The exponential law uses the pressure, so the material must be loaded in
compression for a positive rate. Tests:
`validation-suite/test-2026/group_materi_plasti_visco_exponential_name1`-`5`,
`_name_error1`, `_limit1`, `_limit2`.

## Differences with the Professional version

- The Professional 02-08-2026 reads `_name` as one integer and refuses
  model names ("Integer values are expected"), so the named exponential
  data cannot be used there. The GNU follows the manual (DIVERGENCES
  VISCO-EXPONENTIAL-NAME).
- `_values` must hold exactly two values per name (error otherwise).
- The exponential law without the limit active already differs from the
  Professional for some alpha (BACKLOG C117); the tests use alpha values
  where the plain record equals the Professional.
