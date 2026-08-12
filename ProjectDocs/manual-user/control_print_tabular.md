# control_print_tabular

## Description

`control_print_tabular` exports the nodal results to a tabular format
(CSV and/or SQLite) for programmatic post-processing with pandas, Python,
R or any spreadsheet tool. It is the recommended alternative to
`control_print_plotmtv`/`control_print_gid` when the results are to be
analyzed numerically instead of visualized.

The derived stress magnitudes (von Mises, Tresca, principal stresses)
produced by this keyword are also written by `control_print_vtk` as
`SCALARS vmises/tresca/sig1/sig2/sig3` in the POINT_DATA section when a
stress tensor dof (`materi_stress`) is present.

Each `control_print_tabular` record produces two files next to the output
base name:

- `<base><icontrol>.csv` — one header row (column names) and one row per
  node per time step, in append mode (the file grows across time steps).
- `<base><icontrol>.sqlite` — an SQLite database (requires the build to be
  compiled with `SQLITE_USE=1`; see the developer manual). The CSV is
  always written.

## Uso

Place it in the data part, using the same `icontrol` index as the
`control_timestep`/`control_print` records you want to export:

```
control_print_tabular 0 -yes
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `0`       | Index of the control record. Must match an active `control_timestep` record. |
| `-yes`    | Enable export. Omit or use `-no` to disable. |

## Output

### CSV

The first row is the header. Columns:

- `node` — node number.
- `t` — current time.
- one column per nodal dof: `velx_0, velx_1, ...` for vector dofs,
  `sigxx_00, sigxx_01, ...` for matrix (stress) dofs.
- if a stress tensor dof (`materi_stress`) is present: `vmises, tresca,
  sig1, sig2, sig3` — von Mises equivalent stress, Tresca
  (2*tau_max = sig1-sig3), and the principal stresses in descending
  algebraic order (compression-negative in soil mechanics sign
  convention).

Each call appends one row per node at the current time. To get a clean
single time step, remove the file before running.

### SQLite

Tables (long format):

- `coords(node, x, y, z)` — nodal coordinates (0.0 for unused dimensions).
- `primary_data(node, dof, t, value)` — every exported dof, one row per
  (node, dof, t).
- `derived(node, t, vmises, tresca, sig1, sig2, sig3)` — derived stress
  magnitudes, one row per node per time step.
- `user_data(node, variable, t, value)` — user variables computed by
  `tools/postprocess.py user`.
- `meta(key, value)` — metadata (ndim, sign convention, file base).

`node` is the original node number as written in the input mesh.

## Example

Export the results of control step 20 of a biaxial test to CSV and SQLite:

```
control_print_tabular                    20  -yes
control_timestep                         20  0.001 0.04
control_print                            20  -time_current
```

Produces `<base>20.csv` and `<base>20.sqlite`.

Post-process in Python:

```python
import pandas as pd
df = pd.read_csv("biax20.csv")          # CSV
con = sqlite3.connect("biax20.sqlite")  # SQLite
df = pd.read_sql("SELECT * FROM primary_data WHERE dof='sigyy_11'", con)
```
