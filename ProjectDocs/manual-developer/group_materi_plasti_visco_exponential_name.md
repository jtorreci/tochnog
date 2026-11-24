# group_materi_plasti_visco_exponential_name / _values / _limit (developer)

## Where

- `database.cc`: `GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_NAME` (INTEGER
  list, DATA_ITEM_SIZE, not fixed; it was INTEGER length 1),
  `GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_VALUES` (DOUBLE list) and
  `GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL_LIMIT` (DOUBLE, 1), all
  `data_class` MATERI. The singular `..._value` alias is in `db_number`.
- `stress.cc` `set_stress()`: the names are read in the else-if chain after
  the power names (`visco_named = 2`; the power names set 1). The pressure
  is taken from `new_sig` there, as for the plain exponential record. At
  the first yield test (`plasti_iter==1`) `visco_named_select()` (the D2
  `visco_power_select()`, renamed) returns the pair of the active model;
  for `visco_named==2` it fills `gamma` / `alpha` and
  `plasti_visco_exponential[]`, which selects the exponential branch of
  the lambda update.
- The limit is read with `get_group_data()` into `visco_limit` (default
  `EPS_VISCO` 3.) and replaces the hardcoded `EPS_VISCO` in both lambda
  updates.

## Verification

1D bar in compression (E 2, vonmises 1, gamma 0.5). Named data = the plain
record (alpha 0.1: eppxx -0.192164, = Professional 02-08-2026 for the plain
record); unused name elasto-plastic; second name gamma 0.2 = plain 0.2;
plain record wins; visco_apply -no gates. Limit 1 with alpha 10: GNU eppxx
-0.453174, Professional -0.453047 (2.8e-4; before the change the GNU
ignored the limit, -3.34759). The Professional refuses the names, so the
named targets come from the GNU plain record.

## Pending

- BACKLOG C117: with the limit inactive (alpha 1 or 0.01 here) the plain
  exponential differs from the Professional; more iterations do not close
  it (the GNU fixed point overshoots to f<0).

## Differences with the Professional version

- DIVERGENCES VISCO-EXPONENTIAL-NAME (class D): the Professional 02-08-2026
  refuses model names in `_name`; the GNU follows the manual 6.745.
- `_values` needs two values per name (error).
- BACKLOG C117 / DIVERGENCES C117-VISCO-EXPONENTIAL: the plain exponential
  law with the limit inactive differs.
