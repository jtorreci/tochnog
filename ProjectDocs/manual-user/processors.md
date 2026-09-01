# processors

## Description

Number of shared-memory CPUs/threads used by the solver element loop
(manual Professional 6.10xx, `processors nproc`; default 1). This is
the Professional name of the GNU record `options_processors`, which the
element loop already reads to size its per-thread accumulation arrays;
registering the alias makes Professional input files work unchanged.

If the implementation does not allow more processors the record is
ignored. Note that it sets the number of THREADS, not physical CPUs.

## Syntax

```
processors nproc
```

`nproc` is a positive integer. Example: `processors 1`.

## Example

```
processors 1
```

## Tests

Corpus tests that parse with this alias: `mpc4`, `interface11` (both
use `processors 1`). They remain RUNFAIL because of other missing
records (`mpc_linear_quadratic`, `mesh_interface_triangle_coordinate`)
— see SEGUIMIENTO-CONVERGENCIA.md.

## Notes

- `processors_used` (validation_14_mesh) is a different, OUTPUT record
  (the number of threads actually used, written by the solver) and is
  not covered by this alias.
- `processors_maximum` and `processors_partition` are separate
  Professional records, not implemented.
