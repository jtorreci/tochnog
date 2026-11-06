# group_materi_plasti_visco_power_name / _value (developer)

## Where

- `database.cc`: `GROUP_MATERI_PLASTI_VISCO_POWER_NAME` (INTEGER list,
  DATA_ITEM_SIZE, not fixed) and `GROUP_MATERI_PLASTI_VISCO_POWER_VALUE`
  (DOUBLE list); enums in `tochnog.h` / `tochnog-mod.h`.
- `stress.cc` `set_stress()`: read after `GROUP_MATERI_PLASTI_VISCO_POWER`
  (else-if, so that record wins). `visco_named` keeps `viscoplasti` 0 until
  the first yield test of a plastic loop pass (`plasti_iter==1`, after
  `GET_YIELD_RULE` set `plasti_type`); `visco_power_select()` maps the
  type to its record and returns eta / p of the matching name.
  `max_plasti_iter` switches between `max_iter_visco` (10, 2 with
  `materi_plasti_f_nonlocal`) and `max_iter_ep` (MAX_PLASTI_ITER); a
  user maximum overrides both, as before.

## Implementation details

- Names are compared with `labs()`, so `-group_materi_plasti_x` and its
  stored sign both match.
- Sub-surface types: DRUCKPRAG_APEX -> DRUCKPRAG, MATSUOKANAKAI_APEX ->
  MATSUOKANAKAI, MOHRCOUL_APEX / _01 / _12 / _20 -> MOHRCOUL.
- f_ref is 1 (Professional layout eta p).

## Verification

Professional 02-08-2026 probes on the 1D bar of `materi_plasti_visco_apply1`
(tension 1, E 2): named tension = visco_power (sigxx 2.666667), unused name
= elasto-plastic (1), second name with eta 0.2 (2.428571; GNU 2.428422 =
its own visco_power 0.2), visco_power wins (2.666667), visco_apply -no (1).
A full name list equals visco_power also on a tension + vonmises group
(GNU and Professional each), but that combination differs between them
already without visco (BACKLOG C116).

## Pending

- `group_materi_plasti_visco_exponential_name` / `_values` (6.745, 6.746):
  registered, not read (the Professional reads the name as an integer).
