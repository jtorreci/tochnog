# control_input

## Description

`control_input index switch` (manual Professional 6.139): with `-yes`
the program reads an extra piece of input from the file `index.dat` (for
example `10.dat` for index 10) at the start of control index `index`.
Useful to add data produced during the run, for example a mesh written by
an external mesher (Professional examples `test/other/mesh_gid_1.dat` ...
`mesh_gid_8.dat` with `10.dat`, after `control_mesh_gid_batch 10`).

The file holds data-part records closed by `end_data`; the rest of the
file is ignored. Comments `( ... )`, defines, arithmetic and ranges
`-ra ... -ra` are not allowed. Measured on the Professional 02-08-2026 and
done the same way in the GNU:

- the file is read at the start of the index (an index 0 file affects the
  first timestep);
- a record replaces an existing record with the same index; new records
  are added (also nodes and elements: the mesh is extended, new nodes
  start with zero dofs);
- each record is checked as in the data part;
- a missing file is skipped without a message; `-no` does not read it;
- a comment, an unknown word, a defined word, a range or a missing
  `end_data` stop the program with `Error for control_input detected.` /
  `Error detected near: <word>`.

## Usage

```
control_input 1 -yes
```

with `1.dat`:

```
group_materi_elasti_young 0 2000
end_data
```

Tests (run by `scripts/build_safe.sh` in an isolated directory that holds
the `1.dat`): `validation-suite/test-2026/control_input1` (Young's modulus
replaced at index 1: sigyy -50 = Professional), `control_input2` (a node
pair and an element added), `control_input_error1` (comment in the file).

## Differences with the Professional version

- A new element starts its integration-point state from the node values
  (as every element a GNU mesh operation adds); the Professional starts it
  otherwise. In `control_input2` the nodes of the new element only agree
  (sigyy 20); the shared nodes average to -10 in the GNU and -5 in the
  Professional (DIVERGENCES `CTRL-INPUT-NEWELEM`, class C; BACKLOG C96).
- The manual says the piece of input "needs to be closed by two end_data
  statements"; both programs stop at the first `end_data`.
