# generated_element_dof_init

Zero initialization of the `element_dof` template of elements created by
the mesh generators (bug fix, F0 T4, 2026-10-05).

## Implementation

- File `generate.cc`, functions `generate_spring()`,
  `generate_beam_truss()`, `generate_interface()` and
  `generate_interface_triangle()`.
- Each function allocates a `tmp_element_dof` buffer of `mnolnuknwn`
  doubles with `get_new_dbl()` and writes it as the `ELEMENT_DOF` record
  of every generated element (`db( ELEMENT_DOF, ..., PUT )`).
- `get_new_dbl()` does not initialize memory. The buffer was written
  as-is, so the `element_dof` of generated elements carried heap garbage
  (denormals such as `6.9e-310` / `7.6e-321` in the printed `.dbs`).
- Fix: `array_set( tmp_element_dof, 0, mnolnuknwn );` right after the
  allocation in the four functions, mirroring `create_element()` in
  `create.cc` (which already zeroes its own `tmp_element_dof`).

## Verification

- Two runs of the new binary give byte-identical `.dbs` for `gen1`,
  `genbeam1`, `iface_gen*` (test-2014) and `conspr3`, `genera1`,
  `genera2`, `ho_othr1`, `spring2`, `spring3`, `spring5` (sfnet); no
  `e-3xx` values remain (before: 1 to 574 denormal lines per file).
- A/B against the previous binary: the only changed records are
  `element_dof` lines of generated elements (now all zero).
- Full corpora old vs new (clean builds of `b147ee4` and of the fix):
  test-2014 233/1 (the 1 is `print_mesh_dof.dat`, an output dump, not
  an input), test-2026 47/0, sfnet 187/10 — identical per-test status;
  every changed output is `element_dof`-only except wall-clock
  `time_calculation` lines and the pre-existing `dofid_no` `dof.21`
  garbage (BACKLOG C23).

## Notes

- The garbage sat in `element_dof` slots the calculation does not
  overwrite or read (it was still present in the final `.dbs`), so
  results (node dofs, targets) do not change. Only the printed record
  and the run-to-run determinism of the `.dbs` are affected.
- No new enums, records or hardcoded parameters.

## Pending

- None.
