# group_materi_plasti_visco_power_name / group_materi_plasti_visco_power_value

## Description

Visco-plasticity data per plasticity model (manual Professional 6.748,
6.749). `_name` lists plasticity models of the group
(`-group_materi_plasti_mohr_coul`, `-group_materi_plasti_tension`, ...);
`_value` gives eta and p of the power law (as
`group_materi_plasti_visco_power`) for each name, in the same order.

When an integration point yields, the active plasticity model is looked up
in the names: a named model is evaluated elasto-viscoplastic with its eta
and p, a model that is not named stays elasto-plastic. Useful when, for
example, only the shear failure should be rate dependent and the tension
cut-off not.

If `group_materi_plasti_visco_power` is also given, it wins and the names
are ignored (as the Professional). `materi_plasti_visco_apply -no` /
`control_materi_plasti_visco_apply -no` also switch the named data off.

## Syntax

```
group_materi_plasti_visco_power_name index -name_0 -name_1 ...
group_materi_plasti_visco_power_value index eta_0 p_0 eta_1 p_1 ...
```

Sub-surfaces count for their record (`_druckprag` apex, `_matsuokanakai`
apex, the `_mohrcoul` corners).

## Example

```
group_materi_plasti_tension 0 1.
group_materi_plasti_visco_power_name 0 -group_materi_plasti_tension
group_materi_plasti_visco_power_value 0 0.1 1.
```

Tests: `validation-suite/test-2026/group_materi_plasti_visco_power_name1`-`5`,
`_name_error1`.

## Differences with the Professional version

- `_value` must hold exactly two values per name (error otherwise); the
  Professional 02-08-2026 runs with a missing value.
- Groups that combine several plasticity models already differ from the
  Professional without visco data (BACKLOG C116), so the named data were
  verified on single-model groups.
