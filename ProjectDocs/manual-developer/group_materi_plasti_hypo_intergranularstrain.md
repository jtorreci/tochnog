# group_materi_plasti_hypo_intergranularstrain (variable-length, 5 or 6 params)

## Implementation

- **Registration** (`database.cc`,
  `GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN`): `data_length = 6`
  (maximum) with `fixed_length = 0` (variable). A short 5-value list
  stops at the next keyword instead of swallowing it.
- **Kernel buffer** (`hypoplas.cc`): `LENGTH_INTERGRANULARSTRAIN` is 6.
- **Consumption** (`hypoplas.cc`, `hypo_...` group branch): the record
  is read with `GET` (not `GET_AND_CHECK`), which returns the actual
  stored length in `length_intergranularstrain`. The consumer maps
  `[0..4]` (`epi_R`, `epi_mr`, `epi_mt`, `epi_betar`, `epi_chi`) and
  ignores the optional `[5]` (`theta`): the kernel evaluates the
  `rho^theta f_d N S_hat` term with `chi`, which the manual recommends
  for monotonic loading (`theta = chi`).
- **Alias** (`database.cc`, `db_number`):
  `group_materi_plasti_hypo_strain_intergranular` (Professional spelling)
  -> `GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN`.

## History / gotcha

- The original GNU (jan-2014 sources) registered `data_length = 5`
  (fixed). Commit `22801c3` (2026-09-01) changed it to 6 per the
  Professional manual 6.715 and left `fixed_length` at its default of 1.
  The parser then read six values and swallowed the next keyword
  (`post_point`), failing hypo2, hypo3 and hypo4.
- The sfnet corpus gives **5** values (no `theta`); the validation suite
  gives **6**. A fixed registration cannot satisfy both, hence the
  variable-length form.
- The Professional parser is variable-length: the Professional binary
  (02-08-2026) gets past the data part of hypo2 (its error is in the
  initialisation part, not the intergranular record).

## Verification

- hypo2, hypo3, hypo4 of the sfnet corpus (5 values): rc=0 after the fix.
- Validation-suite hypo2 (6 values): rc=0 (needed by the build_safe
  suite, which would fail with a fixed `data_length = 5`).
- Corpus A/B: see the SEGUIMIENTO entry (parse-arity family).

## Differences with the Professional version

- **Strain reversal with intergranular strain (`hypo3`)** (deliberate
  divergence, class A, project decision 2026-10-12; DIVERGENCES K14b).
  Oedometric von Wolffersdorff + intergranular strain test with a strain
  reversal (test-2014 `hypo3`). Converged reference (independent Python
  integration of the same law, manual 2.2.4): sigxx -0.10466. GNU: -0.1034
  at dt 1e-3 (1.2 % off), -0.104235 at dt 1e-5. Professional 02-08-2026:
  -0.1629 at dt 1e-3 (56 % off), -0.10499 at dt 1e-5: a time-discretisation
  error at the reversal that vanishes only as dt -> 0. Guard tests
  `test-2026/hypo_wolff_igs_oedometer1/2` (wired in `scripts/build_safe.sh`;
  tolerances 2e-3 / 5e-4 absolute). No source change: the GNU substepping
  of the kernel `hypo.c` already converges to the reference. The printed history slots
  (`hyhis*`) of `hypo3` still differ (AUDIT K14, open). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
