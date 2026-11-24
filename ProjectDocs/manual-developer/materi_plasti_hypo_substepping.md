# materi_plasti_hypo_substepping

## Files and functions

- `database.cc`: `MATERI_PLASTI_HYPO_SUBSTEPPING` (F5 T3 C1, INTEGER,
  length 1, class MATERI, `no_index`), next to the existing
  `CONTROL_MATERI_PLASTI_HYPO_SUBSTEPPING` (registered since the GNU 2014
  but never read before F5 T3 C1). Enum appended in `tochnog.h` /
  `tochnog-mod.h` (F5 T3 C1 block).
- `hypoplas.cc`, `hypoplasticity()`: sets the thread-local
  `hypo_substepping_off` with `materi_apply_gate_off(
  MATERI_PLASTI_HYPO_SUBSTEPPING, CONTROL_MATERI_PLASTI_HYPO_SUBSTEPPING )`
  (`general.cc`: the control record of ICONTROL wins over the global one).
- `hypo.c`, `hypo_()`: the accuracy halving test (stress increment > 2 %
  of the stress) is skipped when `hypo_substepping_off` is set. The
  halving after a failed `sigma_()` stays.

## Implementation details

- The flags live in `hypo.c` as `__thread` C globals and are declared
  `extern "C" __thread` in `hypoplas.cc`: the element loop runs in
  parallel, so a plain global would race.
- Masin kernels (`masin_basic.c`, `masin_umat`) use a fixed substep count
  (`MASIN_BASIC_NSUB`) and ignore the switch.

## Verification

Professional 02-08-2026 probes on its `other/hypo1`: default -862.9291,
`-no` -863.7163, `control ... 20 -no` -863.7163, global `-no` + control
`-yes` -862.9291; GNU -862.9255 / -863.7163 / -863.7163 / -862.9255.
Tests `materi_plasti_hypo_substepping1`, `_2`.

## Pending

- The 2 % tolerance and the minimum substep `2e-6 dtime` are hardcoded in
  `hypo.c`.
- No switch in the Masin routes (DIVERGENCES HYPO-SUBSTEPPING-MASIN).
