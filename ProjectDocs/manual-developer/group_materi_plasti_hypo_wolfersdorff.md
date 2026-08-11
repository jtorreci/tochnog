# group_materi_plasti_hypo_wolfersdorff (and variants)

## Files and functions

- `hypo.c` — the hypoplasticity kernel (pure-C port of the original Fortran
  `hypo.f`, then f2c, then manually ported to C in P4-F). Contains:
  - `hypo_()` — entry point: parameter setup, pressure-dependent void-ratio
    initialisation, intergranular-strain substepping, calls `sigma_`.
  - `sigma_()` — the constitutive evaluation. `ihypotype` selects the law:
    0 = wolfersdorff, 1 = lowangles (uses `rval`=data[9], `powxi`=data[10]).
    `cohesion` is subtracted from the normal stresses (tcohesion), and the
    linear contribution is used when the pressure falls below `-3*c`.
  - Helpers: `inpro_`, `normvec_`, `power_`, `tra_`, `extract_`, `mul_`,
    `copy_`, `dev_`, `abdyadic_`, `unity4_`, `zero_`, `add_`, `minus_`.
- `hypoplas.cc` — dispatch:
  - `hypo_type[0]=0` for wolfersdorff, `=1` for lowangles.
  - `cohesion[0]` from `GROUP_MATERI_PLASTI_HYPO_COHESION`.
  - `use_epi[0]=1` + `epi_R/mr/mt/betar/chi` from
    `GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN` (requires
    `materi_strain_intergranular`).
  - `use_pres[0]=1` from `GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO
    -yes`.
- `database.cc` — registrations:
  - `GROUP_MATERI_PLASTI_HYPO_WOLFERSDORFF`: DOUBLE, length 8.
  - `GROUP_MATERI_PLASTI_HYPO_LOWANGLES`: DOUBLE, length 10.
  - `GROUP_MATERI_PLASTI_HYPO_COHESION`: DOUBLE, length 1.
  - `GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN`: DOUBLE, length 5.
  - `GROUP_MATERI_PLASTI_HYPO_PRESSUREDEPENDENTVOIDRATIO`: INTEGER, length 1.
- `tochnog.h` / `tochnog-mod.h` — enums (kept in sync).

## Implementation details

- `hypo_` requires `materi_history_variables >= 4` (void ratio in `his[0]`,
  intergranular strain in the remaining slots).
- The kernel uses the soil-mechanics convention (compression positive) — the
  caller (hypoplas.cc) converts the tochnog tensors accordingly.
- `lowangles` is the same law as wolfersdorff but with the exponent `xi` and
  `c1/c2` computed from `rval`/`powxi` (the "low angles" variant of Herle).
- `cohesion` is subtracted from the normal stresses BEFORE the law evaluation
  (tcohesion), and the linear part is used if the resulting pressure is below
  `-3*c` — this stabilises zones of small stress (free surfaces).
- `pressuredependentvoidratio -yes` activates the initialisation of the void
  ratio from the current mean stress at `time[1]==0` (hypo_ lines 120-139):
  `e = e_0 * exp(-(-trace/h_s)^n)` clamped to `[e_d, e_i]`.

## Validation

- `hypo1.dat` — wolfersdorff, `sigyy=-863` (validated vs the Herle program).
- `hypo2/3/4.dat` — wolfersdorff + intergranular strain.
- `hypo_cohesion.dat` — wolfersdorff + cohesion `c=5 kPa`, `sigyy=-3945.4`,
  `e=0.5848`. Added 2026-08-11 (P4-A1).
- `hypo_lowangles.dat` — lowangles with `rval=2, powxi=2`, `sigyy=-4263.7`.
  Added 2026-08-11 (P4-A1).
- `hypo_pdvr.dat` — pressure-dependent void ratio, `sigyy=-943.0`,
  `e=0.5622`. Added 2026-08-11 (P4-A1).

## External dependencies

- Core `db()`, `db_active_index()`, `pri()`, `array_add/subtract`.
- The kernel is self-contained C (no f2c, no libf2c) after P4-F.

## Hardcoded parameters / pending refactorings

- The `-3*c` linearisation threshold of cohesion is hardcoded in `sigma_`
  (matches the professional manual).
- `materi_history_variables >= 4` is enforced with an explicit error.
- The `lowangles` exponents `rval`/`powxi` default to 0/1 when the base
  wolfersdorff form is used.
