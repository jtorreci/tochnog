# control_print_vtk

## Implementation

- `top.cc` (control loop): `control_print_vtk` other than `-no` calls
  `print_vtk( icontrol, switch )`; the end-of-run print
  (`print_vtk_calculation`, `miscel.cc` `exit_tn`) calls `print_vtk( -1 )`.
- `print_vt.cc` `print_vtk( icontrol, task=-YES )` (default argument in
  `tochnog.h` / `tochnog-mod.h`) builds the name with `snprintf`:
  `-SEPARATE_INDEX` -> `<data_file_base>_%03ld.vtk` (icontrol),
  `-SEPARATE_SEQUENTIAL` -> `<data_file_base>_%04ld.vtk` from the static
  counter `vtk_sequential` (shared by every index, from 0), otherwise
  `<data_file_base>.vtk`.

## Verification (C99, 2026-11-02)

Professional 02-08-2026 probes (`scratchpad/f5/c99`, two indices 0 and 12,
0.1 x 5 steps): `-yes` -> `ya.vtk` (t 0.5); `-separate_index` ->
`ysi_000.vtk` (t 0.2), `ysi_012.vtk` (t 0.5); `-separate_sequential` ->
`yss_0000.vtk` .. `yss_0004.vtk` (t 0.1 .. 0.5). GNU tests
`test-2026/control_print_vtk_names1`, `_separate_index1`,
`_separate_sequential1` (build_safe, isolated directories). The test-2014
`vtk_*1` inputs that compare two indices use `-separate_index`
(`vtk_coord1_030.vtk` ...).

## Hardcoded / pending

- Widths 3 (index) and 4 (counter) are the Professional ones; the manual
  6.339 has none (DIVERGENCES `VTK-TECPLOT-NAMES`, class A).
- `zip` / `control_zip` gzip `vtk*` (manual wording), which does not match
  `<base>.vtk`; not changed here.
