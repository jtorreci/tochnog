# group_materi_plasti_visco_power

## Where implemented

- `database.cc` (db_initialize): record metadata. `data_length = 3`
  (storage width), `fixed_length = 0` — variable length so BOTH the
  Professional layout `eta p` (manual 6.748) and the legacy GNU layout
  `eta p f_ref` parse.
- `stress.cc` (`set_stress`, read block ~line 471): the record is read
  with a dedicated length variable (`visco_power_length`); `f_ref` is
  taken from the stored third value when present, otherwise set to 1.0
  (the Professional power law has no reference stress).

## Implementation details

The plastic multiplier of the viscoplastic branch is

```
lambda = eta * scalar_power(f/f_ref, p)          (power model)
lambda = gamma * (-pressure) * exp(alpha*f)      (exponential model)
```

and the plastic strain increment is `inc_epp = lambda * dtime *
plasti_dir` (`stress.cc`, "plastic strain part"). With `f_ref = 1` the
power law reproduces the manual 2.2.11 form `eps_dot^plas_kl =
eta*(f)^p * d f_flow/d sigma_kl`.

The two record branches are exclusive: `group_materi_plasti_visco_
exponential` (fixed 2-value) is tested first, then the power record.

## Verification

Corpus tests `visc_pl1` (eta=1.e-1, p=1: stress overshoots the tensile
strength and relaxes to the yield surface; Professional .dbs sigxx =
1.0243 vs GNU within the 0.1 target), `visc_pl2` (eta=1.e-6: creep
negligible over the analysed time span, elastic sigxx = 2.0 within
1.e-3) and `validation_12` (nonlocal tension with two material zones
and kap hardening; target kap = 0.0065 met). Legacy 3-value inputs
(examp19) keep the `f/f_ref` role.

## Pending

- `group_materi_plasti_visco_power_name` / `_value` (manual 6.749/6.750,
  per-plasticity-model visco parameters) are not implemented (records
  not registered).

## Hardcoded data

- `f_ref = 1.` when the Professional 2-value layout is used.
- `EPS_VISCO` caps the exponential argument (shared with the
  exponential law); `EPS_LAMBDA` is the plastic seed multiplier.
