# group_materi_plasti_visco_exponential / _power / _always

## Files and functions

- `stress.cc` — `set_stress()`:
  - Lines 100-115: reads the records and sets `viscoplasti = 1`:
    - `GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL` (2 doubles: `gamma, alpha`),
      computes `pressure = (sig0+sig4+sig8)/3`; `db_error()` if `alpha<=0` or
      `gamma<=0`.
    - `GROUP_MATERI_PLASTI_VISCO_POWER` (3 doubles: `eta, p, f_ref`),
      asserts `f_ref!=0`.
  - Line 118: `GROUP_MATERI_PLASTI_VISCO_ALWAYS` (INTEGER, default `-NO`).
  - Lines 737-755: the plastic branch. When `viscoplasti`:
    ```
    lambda_previous = 0.; f_previous = f;
    if ( plasti_visco_exponential[0]>0. ) {
      tmp = alpha*f;
      if ( tmp>visco_limit ) tmp = visco_limit;
      lambda = gamma * (-pressure) * exp(tmp);
    } else {
      assert( plasti_visco_power[0]>0. );
      lambda = eta * pow(f/f_ref, p);
    }
    ```
  - The same `visco_limit` clamp is applied on the second pass
    (predictor/corrector update).
- `stress.cc` — `#define EPS_VISCO 3.` is the default of `visco_limit`,
  read from `group_materi_plasti_visco_exponential_limit` (F5).
- `check.cc:604-611` — both visco-plasti records require `materi_stress` and
  `materi_strain_plasti`.
- `database.cc:2713-2725` — registrations:
  - `GROUP_MATERI_PLASTI_VISCO_ALWAYS`: `INTEGER`, length 1.
  - `GROUP_MATERI_PLASTI_VISCO_EXPONENTIAL`: `DOUBLE_PRECISION`, length 2.
  - `GROUP_MATERI_PLASTI_VISCO_POWER`: `DOUBLE_PRECISION`, length 3.
- Enums mirrored in `tochnog.h` and `tochnog-mod.h`.

## Implementation details

- The viscoplastic multiplier `lambda` feeds the same `plasti_rule()`
  flow-rule machinery as rate-independent plasticity; the difference is only
  how `lambda` is computed (rate law vs. consistency condition). For
  viscoplasticity the plastic flow direction is NOT normalized
  (`if (!viscoplasti) array_normalize(...)`).
- `viscoplasti_always==-YES` makes the branch run even when `f<EPS_F`,
  effectively applying the rate law from the first iteration.
- When `materi_plasti_f_nonlocal` is active, `f` is replaced by the nonlocal
  yield function unknown `fn_indx` (line 768).
- The exponential `alpha*f` is capped at `visco_limit` (default 3) to avoid
  overflow of `exp()`.

## Missing / not implemented (professional parity)

- None left: `_limit`, `_name` / `_values` and the power names are
  implemented (F5 T3, see `group_materi_plasti_visco_exponential_name.md`
  and `group_materi_plasti_visco_power_name.md`).

## External dependencies

- Core `db()`/`get_group_data()`/`db_error()` accessors.
- `plasti_rule()` (plasticity machinery) and globals `viscoplasti`,
  `viscoplasti_always`, `pressure`, `gamma`, `alpha`, `eta`, `p`, `f_ref`,
  `materi_plasti_f_nonlocal`, `fn_indx`.

## Hardcoded parameters / pending refactorings

- BACKLOG C117: with the limit inactive the plain exponential law differs
  from the Professional 02-08-2026 for some alpha.
- The exponential law uses `-pressure` from the current stress; for purely
  tensile paths this may give `lambda=0` — a documented limitation.

## Time step (F4 T3, 2026-10-09)

`set_stress()` (stress.cc) integrates `lambda*dtime` with the time step
of the CURRENT step (`DTIME` VERSION_NEW). It read VERSION_NORMAL (the
previous step; 0 in the first step) before BACKLOG C33 was resolved.
Verified on the Professional 02-08-2026 (1D bar, tension 1, visco_power
0.1 1, dt 0.5 then 1: sigxx 2.666667, old GNU 2.818182); test
`test-2026/materi_plasti_visco_power_dt1`.
