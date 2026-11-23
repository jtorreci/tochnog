# bounda_time_smc family (SMC accelerograms)

Base excitation from a Strong Motion CD file (manual Professional 6.42, 6.43, 6.44).

## Syntax

```
bounda_dof               10  1 -accx
bounda_time_smc          10  -yes
bounda_time_smc_offset   10  1.          (optional, default 0)
bounda_time_smc_units    10  3600. 100.  (factor_time factor_length, required)
```

With `bounda_time_smc index -yes`, Tochnog reads the file `<index>.smc` in the working
directory, so index 10 reads `10.smc`. With `-no` the record is ignored.

## The file

The file follows the USGS SMC format (http://nsmp.wr.usgs.gov/smcfmt.html):

- 11 text lines. The first one names the data: `ACCELEROGRAM` (uncorrected or corrected),
  `VELOCITY` or `DISPLACEMENT`. Any other first line is an input error
  (`Error, dof format SMC file.`).
- 48 integers. The 16th is the number of comment lines.
- 50 reals. The 2nd is the number of samples per second.
- The comment lines, then the samples up to the end of the file. The number of samples in
  the integer header is not used.

The samples may be separated by white space (`1.234e0 1.234e0`) or glued in fields of 10
columns (`-3.1556E+0-5.6333E+0`). CRLF and LF line ends give the same samples.

## How the samples are applied

- Sample k lies at time `time_offset + k / (samples_per_second * factor_time)`.
- The SMC units are cm and seconds. A sample is multiplied by:

  | Data | Factor |
  |---|---|
  | acceleration | factor_time² / factor_length |
  | velocity | factor_time / factor_length |
  | displacement | 1 / factor_length |

  For hours and meters use `3600. 100.`.
- Between two samples the value is interpolated linearly.
- Before the first sample (time < time_offset) and after the last one, the dofs of the
  record are free.
- The dofs of `bounda_dof` must match the data of the file: an acceleration file needs
  `-accx`/`-accy`/`-accz` (with `materi_acceleration`), and likewise for velocity and
  displacement. Otherwise the run stops with
  `Error, the SMC file is for type: materi_acceleration` /
  `The BOUNDA_DOF is for: materi_velocity` / `That is not consistent.`
- `bounda_baseline_correction` works on the whole SMC signal, the same as on a
  `bounda_time` table; see [bounda_baseline_correction](bounda_baseline_correction.md).

Errors: a missing file stops with `Error opening : <index>.smc`; a missing
`bounda_time_smc_units` record with `Error detected for data item : bounda_time_smc_units`.

Example (the Professional example other/smc_2, the 1990 Upland accelerogram at 200 samples
per second): with units `3600. 100.`, a truss with node 1 bound by `-accx` reaches disx
0.0837695 at t = 6e-3 h, the same as the Professional 02-08-2026.

Tests: `bounda_time_smc_constant1`, `_upland1`, `_offset1`, `_free1`, `_velocity1`,
`_displacement1`, `_glued1`, `_glued2`, `_error1` to `_error4`.

## Differences with the Professional version

SMC-GLUED-LINE-END: when a file has LF lines of 80 glued columns (8 samples without spaces),
the Professional 02-08-2026 adds a zero sample after each line; with CRLF lines it does not.
The GNU version reads the same samples whatever the line end. Guard test
`bounda_time_smc_glued2`.
