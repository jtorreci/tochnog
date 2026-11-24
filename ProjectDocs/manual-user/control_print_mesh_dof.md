# control_print_mesh_dof

## Description

`control_print_mesh_dof` is the Professional name of the GNU keyword
`print_mesh_dof`: a one-shot dump of the node coordinates and the listed
dof values, written at the first evaluation of the calculation. The
dump is written to the file `print_mesh_dof.dat`, with one line per
node: node number, coordinates and the requested dof values (all dofs
when none are listed).

## Uso

```
control_print_mesh_dof -disy
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| (no index) | The record is not indexed (like the GNU `print_mesh_dof`). |
| `data_item_name ...` | Optional list of dof labels (e.g. `-disy`, `-temp`) whose values are added to the dump. When omitted, no dof values are printed. |

## Output

`print_mesh_dof.dat` — one line per node: `node x y [z] <dof values>`.
The dump is written once, at the first evaluation (before the first
step).

## Example

```
control_print_mesh_dof -disy
```

Produces `print_mesh_dof.dat` with the node numbers, coordinates and
the current `disy` value of every node.

## Validation

Test `meshdoff` (validation-suite/test-2014): smoke — the dump file is
generated with the Professional keyword and its first line is the first
node (`1 0 0` for node 1 at the origin in 2D).

## Differences with Professional

- The Professional syntax includes an index and a switch (`index
  switch`); the GNU record is not indexed and takes the dof list
  directly. The alias keeps the GNU layout (documented difference).
- The dof values of the dump are the values at the first evaluation
  (all zero in a fresh model before any load is applied).
