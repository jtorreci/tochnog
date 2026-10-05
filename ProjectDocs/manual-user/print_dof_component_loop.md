# print_dof_component_loop

## Description

Bug fix (2026-10-05) of the `dof.<index>` files written by
`control_print_dof`. The files now contain exactly the documented layout:
for every print, one line per node for each dof COMPONENT (vector dofs
per component `x`, `y`, `z`; the stress per Voigt component `xx yy zz xy
xz yz`), each component once.

Before the fix the vector and stress blocks were repeated (once per
stored value and derivative) and shifted into the next dof, and some
lines held random memory values such as `6.9e-310` that changed from run
to run (most visible with `derivatives` in the initialization part).
Files written by older binaries are therefore longer and not comparable
line by line.

## Syntax

No new input: the fix applies to every `control_print_dof` record,
including `control_print_dof_id` and `control_print_dof_smooth_dof`.

## Example

The validation-suite input `dofid_no` (2D, `derivatives`, velocity,
integrated velocity, displacement and stress on 4 nodes) writes `dof.21`
with 12 components x 4 nodes = 48 lines of `x y <value>` per print
(768 lines with garbage before the fix).

## Note

The FRD (`control_print_frd`) and gmsh (`control_print_gmsh`) outputs
still repeat the vector/stress result blocks; they are scheduled for the
same fix (BACKLOG C24).
