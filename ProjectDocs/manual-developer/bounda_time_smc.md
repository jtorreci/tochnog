# bounda_time_smc family (SMC accelerograms)

Developer notes.

## Where

- `tochnog.h`/`tochnog-mod.h`: BOUNDA_TIME_SMC, BOUNDA_TIME_SMC_OFFSET,
  BOUNDA_TIME_SMC_UNITS enums.
- `database.cc`: registration (INTEGER switch / DOUBLE offset / DOUBLE x2
  units, class BOUNDA) next to the bounda_time_* records.
- `bounda.cc`: explicit runtime failure when a bounda_time_smc record is
  active (the reader is PENDING).

## What the reader must do (from the manual + Professional binary probing)

1. Open `<index>.smc` (the bounda_time_smc index) in the working
   directory.
2. Parse the SMC format: header lines, the "-32768" (missing) and
   1.7e38 (absent) placeholders are skipped; the remaining samples form the
   acceleration signal.
3. Time axis: sample k at time (offset + k*dt_file), converted by
   factor_time; acceleration samples converted by 1/factor_length
   (cm -> model length).
4. The signal feeds the -accx prescription (see materi_acceleration: the
   bound value is v_new = v_old + a*dt with a interpolated from the
   signal).

Probing notes: the corpus smc_1.smc is a synthetic constant-1.234 signal
(target identical to dynamic3, 61.7 +/- 1); smc_2.smc is a real Upland
accelerogram (target disx 0.08376 +/- 1e-5 with units 3600./100.).
Both tests bind the SMC record with index 10 -> they need a file literally
named `10.smc` at runtime (the corpus stores them as smc_1.smc/smc_2.smc),
which the corpus harness does not stage.
