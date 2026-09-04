# bounda_time_smc family (SMC accelerograms)

Base-acceleration input from Strong Motion CD files (manual Professional
6.42/6.43/6.44).

## Syntax

```
bounda_dof               10  1 -accx
bounda_time_smc          10  -yes
bounda_time_smc_offset   10  0.          (optional, default 0)
bounda_time_smc_units    10  1. 1.       (factor_time factor_length)
```

When `bounda_time_smc index -yes` is set, Tochnog reads the file
`<index>.smc` (the index of the bounda_time_smc record, so `10` reads
`10.smc`) next to the input file. The file follows the SMC format of
http://nsmp.wr.usgs.gov/smcfmt.html and contains base acceleration time
data. The SMC units are cm and seconds; `factor_time` converts the time
axis and `factor_length` the acceleration data to the units of the input
file (e.g. `3600. 100.` for hours and meters: cm/s2 -> m/s2, s -> h).

## Implementation status

PENDING — the three records are registered and parsed, but the SMC file
reader is not implemented. A run that activates `bounda_time_smc -yes`
aborts with an explicit "not implemented (PENDING)" error instead of
silently applying a zero acceleration.

The corpus tests smc_1/smc_2 additionally depend on the file `10.smc`
being present in the working directory (the SMC index is 10 in both tests);
the corpus harness stages only the `.dat` file.
