# control_print_dof_line

## Description

`control_print_dof_line` prints the values of the `node_dof` records and
the `node_dof_calcul` records (manual Professional 6.273) **along a line
in space** to files. The line is a POLYLINE: the start point of the first
segment is `x_0 y_0 z_0` and its end point is `x_1 y_1 z_1`; the second
segment starts at `x_1 y_1 z_1` and ends at `x_2 y_2 z_2`, etc. In 1D only
the x-coordinates are needed.

The printed files contain lines like `x y z <dof>` (in 1D only `x`),
where `dof` is the value of the dof, e.g. `temp`. One file is written per
dof label.

The companion records of this record are:

- [`control_print_dof_line_coordinates`](control_print_dof_line_coordinates.md) — the polyline vertices (required).
- [`control_print_dof_line_n`](control_print_dof_line_n.md) — number of points along the line (optional).
- [`control_print_dof_line_element_group`](control_print_dof_line_element_group.md) — restrict the interpolation to element groups (optional).
- [`control_print_dof_line_eps_iso`](control_print_dof_line_eps_iso.md) — point-in-element tolerance (optional).
- [`control_print_dof_line_method`](control_print_dof_line_method.md) — coordinate frame (optional).
- [`control_print_dof_line_move`](control_print_dof_line_move.md) — follow material particles (optional).
- [`control_print_dof_line_time`](control_print_dof_line_time.md) — write the time as a gnuplot comment (optional).

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_timestep` record:

```
control_print_dof_line            30  -separate_index
control_print_dof_line_coordinates 30  0.5 0. 0.5 1.
control_print_dof_line_n          30  3
control_timestep                  30  1.0 1.0
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `30`      | Index of the control record. Must match an active `control_timestep`. |
| switch    | `-separate_index`: filenames like `dof.index`, e.g. `temp.10`, `velx.10` (the "dof" is the dof LABEL). `-separate_sequential`: filenames sequentially numbered, `temp.0`, `temp.1`, ... one number per call. `-yes` behaves like `-separate_index`. |

## Output

- One file per dof label (e.g. `disy.30`) and per `node_dof_calcul` item,
  containing one `x y z <dof>` line per printed point. Each call appends
  the current state, so the file grows as a time series (one block of
  `n` lines per time step).
- Points that are not accepted as part of any element (outside the mesh
  or outside the filtered element groups) are omitted from the files.
- With `control_print_dof_line_time ... -yes` the first line of each
  block is `# time <time_current>` (gnuplot comment format).

## Example

A vertical line in a 2D model with `n = 3`:

```
control_print_dof_line            30  -separate_index
control_print_dof_line_coordinates 30  0.5 0. 0.5 1.
control_print_dof_line_n          30  3
control_timestep                  30  1.0 1.0
```

Produces `disy.30` (and one file per other dof label) with, for each time
step, one block of `x y <disy>` lines for the 3 points `(0.5,0)`,
`(0.5,0.5)` and `(0.5,1)`.

## Differences with the Professional version

- The Professional default `node_type` of `_method` is
  `-node_start_refined`; when no `node_start_refined` record exists the
  stored node coordinates are used (geometrically linear analyses).
- The default of `_n` is 5 (the legacy GNU `post_line_n` default; the
  manual does not state one).
- `-separate_sequential` numbering is a static counter shared per
  process (reset per run).
