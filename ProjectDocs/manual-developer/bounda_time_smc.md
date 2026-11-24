# bounda_time_smc family (SMC accelerograms) — developer

## Where

- `tochnog.h`/`tochnog-mod.h`: the BOUNDA_TIME_SMC, BOUNDA_TIME_SMC_OFFSET and
  BOUNDA_TIME_SMC_UNITS enums.
- `database.cc`: registration (INTEGER switch, DOUBLE offset, DOUBLE x2 units, class BOUNDA).
- `bounda.cc`:
  - **`bounda_time_smc_read()`** (static) reads `<index>.smc` into a time-value table in
    model units and returns the data type (0 acceleration, 1 velocity, 2 displacement).
  - **`bounda()`** has an SMC branch, taken when the switch is `-yes`, after bounda_sine and
    bounda_time_user. It copies the two samples around time_total into `bounda_time` (4
    values, ninc 2), so the ordinary time loop interpolates, and leaves the dofs free outside
    the table.
  - After `bounda_dof` is read, each acceleration, velocity or displacement dof is checked
    against the data type.
  - `bounda_baseline_load()` gets the whole SMC table (table = 1).

## Reader

1. 11 text lines. The data type is taken from the first one, by the words ACCELEROGRAM,
   VELOCITY and DISPLACEMENT.
2. 48 integers and 50 reals with `>>`. The integer 16 is the number of comment lines and
   the real 2 the samples per second.
3. Skip the rest of the last real line and the comment lines.
4. Read every remaining line into white-space words. A word longer than 10 characters is
   cut at the line columns 10, 20, ..., and each piece must be a whole number for `strtod`.
5. Scale and place the samples: t_k = offset + k/(rate*factor_time); value times
   factor_time^2/factor_length, factor_time/factor_length or 1/factor_length.

The file is read again at each call of `bounda()`, like `bounda_time_file`. The table is not
limited by DATA_ITEM_SIZE.

## Professional 02-08-2026 measurement (2026-11-23)

The probes were the examples other/smc_1 and smc_2, plus synthetic ramps (sample k = k) with
each rule varied on its own. The table gives the rule the Professional showed and the probe.

| Item | Professional | Probe |
|---|---|---|
| Sample time | k/rate, first sample at 0 | smc_2 history, ramps at 2 and 4 per second |
| Interpolation | Linear | smc_2 t = 1e-6 h: -4.93954 x 129600 |
| Units | Time / factor_time; value as in the table above | Ramps with units 2. 1. and 2. 4. |
| Offset | Model time, after the unit conversion | Offset 1. with units 2. 1. |
| Outside the signal | Dof free | A force on node 1 moves it before the offset and after the last sample |
| npts (integer 17) | Ignored, read to the end of the file | npts 3 with 10 samples |
| Comment lines | Integer 16 | 3 comment lines |
| Type | Words of the first line; mismatch with the dof is an error | First lines CORRECTED ACCELEROGRAM, VELOCITY, DISPLACEMENT, RESPONSE SPECTRA, FOO |
| Units record | Required (db error) | No `bounda_time_smc_units` |
| Missing file | `Error opening : 10.smc` | No file |

The GNU history equals the Professional one to 3e-12 in every probe, and the two examples
give the same final values: disx 61.7617 and 0.08376952307565.

## Differences with the Professional version: SMC-GLUED-LINE-END (deliberate)

With LF lines of 80 glued columns (`0.0000E+001.0000E+00...`), the Professional inserts a
zero sample after each line: at t = 4 in the ramp at 2 per second, accx is 0 where it should
be 8. The same file with CRLF lines reads correctly. The GNU reads both files alike. Guard
test `bounda_time_smc_glued2`, which the Professional fails.

The baseline fit with SMC differs from the Professional exactly as it does for a
`bounda_time` table. That divergence is BASELINE-CORRECTION-SIGNAL, which was already
recorded.
