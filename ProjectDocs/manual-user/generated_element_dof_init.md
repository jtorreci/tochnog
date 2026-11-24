# generated_element_dof_init

## Description

Elements created by the mesh generators
(`control_mesh_generate_contactspring`, `control_mesh_generate_truss` /
`_beam` / `_trussbeam`, `control_mesh_generate_interface` and
`control_mesh_interface_triangle`) now start with an `element_dof`
record filled with zeros.

Before this fix the initial `element_dof` of a generated element held
random memory content, visible in the printed database (`.dbs`) as tiny
numbers such as `6.9e-310`. They sat in positions the calculation does
not use, so the results (node dofs, targets) never depended on them,
but the printed record changed from run to run.

## Why it is useful

- The `.dbs` of models with generated elements is now reproducible
  (two runs give identical files), so regression comparisons work.
- No input change is needed.

## Syntax

No keyword. It applies whenever a mesh generator creates elements.

## Example

Run `validation-suite/test-2014/gen1.dat` (or any input with a
`control_mesh_generate_*` record) twice: the `element_dof` lines of the
generated elements in the printed
database show `0` instead of random values.
