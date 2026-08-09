# group_materi_maxwell_chain

## Files and functions

- `viscoela.cc` — `visco_elasticity()` (line 22): the whole linear Maxwell
  chain model.
  - Reads `dtime`, `group_materi_maxwell_chain` (pairs `E_m, t_m`), and
    `group_materi_elasti_poisson` (optional).
  - `n = db_len(GROUP_MATERI_MAXWELL_CHAIN, gr, VERSION_NORMAL)/2`.
  - For each chain: calls `visco_elastiticity_nonlinear()` (a stub, see below),
    builds the isotropic `C_matrix` with `em`, computes
    `inc_msig = C.inc_epe * (tm/dtime) - rotated_old_msig[m]` and
    `new_msig[m] = rotated_old_msig[m] + inc_msig*(1-exp(-dtime/tm))`.
  - `formulation==TOTAL`: adds `new_msig[m]` to `new_sig`; `INCREMENTAL`: adds
    only `inc_msig`.
  - Assembles `memmat` from the diagonal entries of `C_matrix` for the element
    tangent.
- Entry point called from `stress.cc:700` (the `visco_elasticity` branch of the
  `materi_stress` block).
- `database.cc:2456-2464` — registrations:
  - `GROUP_MATERI_MAXWELL_CHAIN`: `DOUBLE_PRECISION`,
    `data_length = materi_maxwell_stress*2` (pairs E,t).
  - `GROUP_MATERI_MAXWELL_CHAIN_NONLINEAR`: `DOUBLE_PRECISION`,
    `data_length = materi_maxwell_stress`.
- `check.cc:506-508` — both records require the unknown `materi_maxwell_stress`.
- Enums mirrored in `tochnog.h` and `tochnog-mod.h`.

## Implementation details

- The history variable is the per-chain Maxwell stress `msig`; its storage is
  sized by `materi_maxwell_stress` (see `materi.cc` indexes
  `rotated_old_msig`, `rotated_new_msig`, `old_msig`, `new_msig`, lines
  119-122). Rotation of the old stress is handled in `materi.cc:268-284` and
  359-364 (updated formulations).
- `em` is used both for the chain stress update AND for the tangent `C`, so the
  implementation is a generalized-Kelvin/Maxwell hybrid — each chain
  contributes the same `C` scaled by `(1-exp(-dtime/tm))` to the increment.
- `pois` is read with `GET_IF_EXISTS`; if absent it defaults to 0 in
  `C_matrix`.
- The relaxation update is exact for the linear model given a constant `inc_epe`
  within the step.

## Pending: group_materi_maxwell_chain_nonlinear

- Registered in the database (`database.cc:2462`) and validated in
  `check.cc:507`, but the model routine
  `visco_elastiticity_nonlinear()` (`visconon.cc:22`) is an EMPTY function
  body. The linear path calls it as a no-op.
- What remains to implement it: the nonlinear (e.g. rate-dependent / damage
  coupled) relaxation update in `visconon.cc`, plus wiring its data
  (`data_length = materi_maxwell_stress`) into the loop of `visco_elasticity()`.

## External dependencies

- Core `db()`/`get_group_data()` accessors, `C_matrix()` (miscel.cc), and the
  `materi_maxwell_stress` history infrastructure in `materi.cc`.
- Globals `formulation`, `dtime`, `MDIM`, `TOTAL`, `INCREMENTAL`.

## Hardcoded parameters / pending refactorings

- `data_length = materi_maxwell_stress*2` bakes the "pairs of E,t" layout into
  the database; the manual layout `E_0 t_0 E_1 t_1 ...` matches it.
- The nonlinear chain data is reserved but unused — this is the main gap in the
  viscoelastic family.
