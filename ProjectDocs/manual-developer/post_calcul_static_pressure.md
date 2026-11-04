# post_calcul_static_pressure

## Implementation (F5 T2 V8)

- **Record:** `database.cc`, DOUBLE, length 1, no index, class POST.
- **Consumer:** `calcul.cc`, static helper `calcul_phreatic_coord()`: it
  calls `groundflow_phreatic_coord()` and, when the record exists, sets the
  static pressure to its value and reports a level found. It is used by
  the `STATIC`, `DYNAMIC` and safety branches of `calculate_operat()`; the
  `TOTAL` branch keeps `groundflow_phreatic_coord()`. The record is not
  read by `groundflow_phreatic_coord()` itself, so boundary conditions,
  loads and the material stay unchanged (postprocessing only, as the
  manual says).

## Verification

`post_calcul_static_pressure1` (Professional `ground14` plus the record):
node 1 and node 41 `node_dof_calcul` equal the Professional 02-08-2026,
including node 41 above the level (total 0, dynamic 7, piping -0.928571).
