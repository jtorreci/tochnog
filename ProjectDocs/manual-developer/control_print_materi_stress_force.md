# control_print_materi_stress_force

## Implementación

- New file `print_materi_stress_force.cc`:
  `print_materi_stress_force( icontrol, method )`. Invoked from the
  control loop in `top.cc`, INSIDE the `control_print_frequency_allowed`
  gate (same as every other `control_print_*`), right after
  `control_print_beam_force_moment`:
  ```
  if ( frequency_allowed && db_active_index( CONTROL_PRINT_MATERI_STRESS_FORCE, icontrol, VERSION_NORMAL ) ) {
    db( CONTROL_PRINT_MATERI_STRESS_FORCE, icontrol, ival, ddum, ldum, VERSION_NORMAL, GET );
    print_materi_stress_force( icontrol, ival[0] );
  }
  ```
  `ival[0]` is the method (`-ALL` / `-PRIMARY`); anything else is
  `db_error`. The manual "index" is the RECORD index (icontrol): the
  example "materi_stress_force.100 if the index is 100" means the
  record `control_print_materi_stress_force 100 -all` -> file
  `materi_stress_force.100` (same convention as
  `control_print_beam_force_moment`).
- Registered in `database.cc`: CONTROL_PRINT_MATERI_STRESS_FORCE
  (INTEGER, length 1 = the method). Enum + PRIMARY (the -primary
  switch) in `tochnog.h` / `tochnog-mod.h` in sync.
  `print_materi_stress_force.o` and `calcul_force.o` added to the
  makefile.

## Diseño / decisiones

- **Data source**: `NODE_DOF_CALCUL` (VERSION_PRINT, pattern of
  print_vt.cc: `db_version_copy(VERSION_NORMAL, VERSION_PRINT)` +
  `renumbering(NO)` + `db_highest_index(NODE,...)` + read + delete).
  The -force block is located in `POST_CALCUL_UNKNOWN_OPERAT`: the
  first item with operat FORCE; its slot in NODE_DOF_CALCUL equals its
  item index (ONE slot per item, the flat layout that print_vt.cc
  reads). The print reads `nitems = 9 (2D) / 16 (3D)` consecutive
  slots (`post_calcul_materi_stress_force_items()` in calcul_force.cc).
  No -force block (or an incomplete one, e.g. cut by MCALCUL=20) -> NO
  file is written (documented decision, pattern of
  print_beam_force_moment).
- **File structure**: `materi_stress_force.<icontrol>`, append mode,
  one block per print call. The header comments (lines starting with
  `#`) explain the column structure, as the manual 6.328 requires
  ("the files themselves will contain comments explaining the detailed
  structure"); the column list is read from the GLOBAL
  `post_calcul_names` so the header always matches the actual items.
  First column = the node position in the compacted print version
  (0-based).
- **-all vs -primary**: `-primary` skips the averaged (non-primary)
  nodes through the `msf_node_is_averaged()` hook. LOT 1: no averaging
  exists, the hook returns 0 for every node, so both methods write the
  same lines (documented; the distinction becomes effective with the
  quad9/hex27 averaging of lots 2/3).
- **MCALCUL=20 limitation**: the 3D -force block uses 16 of the 20
  per-node slots of NODE_DOF_CALCUL; combining it with other multi-item
  post_calcul records (e.g. a 6-value -total) aborts with the
  pre-existing "MCALCUL too small" message of
  `parallel_calcul_node` (calcul.cc). 2D (9 items) leaves 11 slots.

## Detalles / gotchas

- The -force family is NODAL: `calculate()` (calcul.cc) rejects active
  POST_LINE_DOF/POST_POINT_DOF/POST_QUADRILATERAL_DOF records for a
  -force record BEFORE the per-node loop (clear error + exit), and the
  stub `post_calcul_materi_stress_force()` defends `inod<0` (the
  post-type branch marker) as well.
- `post_calcul` (and every `post_calcul_materi_stress_force_*` record)
  is `no_index=1`: the input syntax has NO leading index
  (`post_calcul -materi_stress -force`); a leading number becomes a
  record VALUE (length 3 -> the pre-existing `(length%2)` db_error).
- The config validation runs once per -force record in `calculate()`
  (BEFORE the per-node loop), so the errors are clear and cheap.

## Pendiente

- The numerical integration (element side normals, stress integration
  over the sides, reference-point orientation, quad9/hex27 averaging,
  the averaged-node flag for -primary) lands in lots 2/3. The print is
  the functional consumer: it reads whatever the calculation writes.
