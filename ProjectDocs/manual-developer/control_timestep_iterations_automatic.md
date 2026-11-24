# control_timestep_iterations_automatic (variable-length, 2 or 3 values)

## Implementation

- **Registration** (`database.cc`, `CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC`):
  `data_length = 3` (maximum) with `fixed_length = 0` (variable). The
  record is read until the next keyword token, so a short 2-value list
  stops instead of swallowing the next record name.
- **Consumption** (`top.cc`, `CONTROL_TIMESTEP` branch): the record is
  read with `GET`, which returns the actual stored length in `ldum`.
  The length-aware mapping is applied before `ldum` is reused by the
  `CONTROL_TIMESTEP_ITERATIONS_AUTOMATIC_STOP` `GET_IF_EXISTS`:

  ```c
  ratio_criterium = control_timestep_iterations_automatic[0];
  if ( ldum>=3 )
    maximum_timestep = control_timestep_iterations_automatic[2];
  else
    maximum_timestep = control_timestep_iterations_automatic[1];
  ```

  `minimal_timestep` (slot `[1]` of the 3-value form) is registered but
  not consumed: the GNU automatic scheme reduces/grows the step without a
  lower bound.

## History / gotcha

- The original GNU (jan-2014 sources) registered `data_length = 2`
  (fixed): `[0]=ratio_criterium`, `[1]=maximum_timestep`.
- Commit `22801c3` (2026-09-01) changed `data_length` to 3 per the
  Professional manual 6.386 but left `fixed_length` at its default of 1.
  The parser then read three values and swallowed the next keyword. The
  whole sfnet corpus gives **2** values, so condif8, delete1, truss7,
  disk1, ho_mech2 and taylor3 failed with "I don't know what to do with
  : <next keyword>".
- The Professional parser itself is variable-length: with
  `control_timestep_iterations_automatic 10 0.01 0.05` the Professional
  binary (02-08-2026) reports the *following* record
  (`control_options_solver`) as an unknown keyword, i.e. it consumed
  exactly two values.

## Verification

- condif8, delete1, truss7, disk1: rc=0 after the fix (previously the
  data-part error above).
- taylor3, ho_mech2 parse now; they still fail their targets for
  model/solver reasons unrelated to the record arity (taylor3 `kap`
  stays 0, ho_mech2 Bi-CG breakdown), see SEGUIMIENTO.
- Corpus A/B: see the SEGUIMIENTO entry (parse-arity family).
