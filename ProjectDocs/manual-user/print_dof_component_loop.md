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

## FRD and gmsh files (2026-10-05)

The same fix applies to `control_print_frd` and `control_print_gmsh`
(BACKLOG C24). An FRD file now holds one result block per dof per print
(for example `VELO`, `DISP`, `STRESS`, `TOSTRAIN`, `NDTEMP`); in 2D the
third component of a vector block is 0. A gmsh file holds one
`node_`/`element_` field per dof component, named after the dof type:
`node_materi_velocity_0`, `node_materi_velocity_1`,
`node_materi_stress_xx`, ..., `node_condif_temperature`. Before the fix
both files repeated the blocks (68 FRD blocks instead of 5 in the
example below) with shifted, partly random values, and the gmsh names
used the component label with a wrong index (e.g. `node_velx_1`).

Example: `validation-suite/test-2026/print_frd_gmsh_dof_loop1.dat`.

## Differences with the Professional version

- `control_print_dof -separate_index`: the Professional writes one file
  per dof (BACKLOG C26).
- FRD: the Professional writes a `1PSTEP` header per block, an extra
  `ALL` vector component, the names `VELOCITY`, `S11`.., `T11`.., extra
  boundary/reaction blocks, and with `derivatives` one block per
  derivative. gmsh: the Professional writes vectors/tensors as one
  multi-component field (`node_materi_velocity`) and extra
  boundary/rhside fields. The values of the dofs agree (BACKLOG C27).
