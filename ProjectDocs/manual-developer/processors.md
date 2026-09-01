# processors — developer notes

Professional record `processors nproc` (number of threads, default 1).
Alias only — the record itself is the GNU `OPTIONS_PROCESSORS`
(INTEGER, `data_length = 1`, `no_index = 1`).

## Implementation

One line in the Professional-name translation chain of
`db_number()` (`database.cc`):

```c
else if ( !strcmp( str, "processors" ) )
  return OPTIONS_PROCESSORS;
```

Exact-match only: the output record `processors_used` must NOT be
translated (it is a separate Professional record not present in the
GNU; tests that use it stay blocked on it).

## How the record is consumed

- `elem.cc::element_loop()`: `db( OPTIONS_PROCESSORS, 0, &nthread, ... )`
  sizes the per-thread `node_nel`/`node_lhside`/`node_rhside`/
  `node_dof_tmp` buffers for the parallel element loop.
- `area.cc`: with `OPTIONS_PROCESSORS > 1` the
  `node_support_edge_normal_force` accumulation is single-threaded-only
  (warns once).

## Verification

- `mpc4` and `interface11` parse past `processors 1` (previously
  "I do not know : processors"); both still RUNFAIL on their other
  blockers (`mpc_linear_quadratic`, `mesh_interface_triangle_coordinate`).
- Corpus 118 -> 121 PASS (matrix2/temp2/matrix4), no regression.

## Gotchas

- The alias must be an exact `strcmp` — a prefix match on "processors"
  would wrongly translate `processors_used`, `processors_maximum` and
  `processors_partition`.
