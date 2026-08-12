# postprocess.py

## Description

`tools/postprocess.py` is a Python command-line tool for post-processing
the SQLite database produced by `control_print_tabular`. It lets you
inspect, summarize, plot and extract profiles from the results without
opening the CSV manually.

Requires Python 3 with `pandas`; `numpy` and `matplotlib` are needed only
by the `plot`/`line` subcommands.

## Subcommands

```
tools/postprocess.py info   <file>.sqlite
tools/postprocess.py stats  <file>.sqlite [--dof DOF ...] [--csv out.csv]
tools/postprocess.py plot   <file>.sqlite --dof DOF [--node N ...] [--out png]
tools/postprocess.py line   <file>.sqlite --dof DOF --x1 --y1 --x2 --y2
                           [--t TIME] [--tol T] [--csv out.csv] [--plot png]
tools/postprocess.py user   <file>.sqlite --name VAR --expr "PYTHON-EXPR"
```

### info

Prints the database summary: metadata, available dofs, time steps and nodes.

### stats

Descriptive statistics (count, mean, std, min, percentiles, max) of the
selected dofs. By default all dofs are reported. `--dof` may be repeated
and accepts both primary dofs (e.g. `sigyy_11`) and derived magnitudes
(e.g. `vmises`, `tresca`, `sig1..sig3`). `--csv` also writes the table.

### plot

Time series of one dof. By default all nodes are plotted; use `--node`
(repeatable) to restrict. The plot is saved to `--out` (default
`<dof>.png`).

### line

Profile of one dof along a straight geometry line. The line is defined by
its two end points `--x1 --y1 --x2 --y2` (only the x,y projection is
used). Nodes within a perpendicular distance `--tol` (default 0.1) of the
segment are selected; `s` is the distance along the line. `--t` selects a
time step (default: the last one). `--csv` saves the profile, `--plot`
also renders it.

### user

Computes a user variable from a Python expression over the exported dofs
and stores it in the `user_data` table (long format
`(node, variable, t, value)`). The expression may use any column name as a
variable, e.g. `0.5*(velx_0**2+vely_0**2)`. Existing rows for the same
variable are replaced.

## Example

Export the results with tochnog:

```
control_print_tabular 20 -yes
control_timestep      20 0.001 0.04
```

Then, from Python:

```
tools/postprocess.py info biax20.sqlite
tools/postprocess.py stats biax20.sqlite --dof sigyy_11 --dof vmises
tools/postprocess.py plot biax20.sqlite --dof vmises --node 0 --out vmises.png
tools/postprocess.py line biax20.sqlite --dof sigyy_11 \
    --x1 0 --y1 0 --x2 1 --y2 0 --t 0.04 --plot sig.png
tools/postprocess.py user biax20.sqlite --name kin \
    --expr "0.5*(velx_0**2+vely_0**2)"
```
