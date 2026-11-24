# control_print_dof_point

## Description

`control_print_dof_point` prints the values of the `node_dof` records and
the `node_dof_calcul` records (manual Professional 6.281) **in a point in
space** to files. The point is given by the companion record
[`control_print_dof_point_coordinates`](control_print_dof_point_coordinates.md).

The printed files contain lines like `x y z <dof>` (in 1D only `x`),
where `dof` is the value of the dof, e.g. `temp`. One file is written per
dof label. In practice the file holds ONE line per time step per dof: a
time series of the value at the point.

The companion records of this record are:

- [`control_print_dof_point_coordinates`](control_print_dof_point_coordinates.md) — the coordinates of the point (required).
- [`control_print_dof_point_time`](control_print_dof_point_time.md) — write the time as a gnuplot comment (optional).

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_timestep` record:

```
control_print_dof_point            44  -separate_index
control_print_dof_point_coordinates 44  0.5 0.5
control_timestep                   44  0.1 0.2
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `44`      | Index of the control record. Must match an active `control_timestep`. |
| switch    | `-separate_index`: filenames like `dof.index`, e.g. `temp.10`, `velx.10` (the "dof" is the dof LABEL). `-separate_sequential`: filenames sequentially numbered, `temp.0`, `temp.1`, ... one number per call. `-yes` behaves like `-separate_index`. |

## Output

- One file per dof label (e.g. `disy.44`) and per `node_dof_calcul` item,
  containing one `x y z <dof>` line per call (a time series of the value
  at the point).
- If the point is not accepted as part of any element (outside the mesh)
  no data line is written.
- With `control_print_dof_point_time ... -yes` the first line of each
  block is `# time <time_current>` (gnuplot comment format).

## Example

```
control_print_dof_point            44  -separate_index
control_print_dof_point_coordinates 44  0.5 0.5
control_timestep                   44  0.1 0.2
```

Produces `disy.44` with one `x y <disy>` line per time step at
`(0.5, 0.5)`.

## Differences with the Professional version

- `-separate_sequential` numbering is a static counter shared per
  process (reset per run).
- A point outside the mesh (not accepted within the point-in-element
  tolerance) produces no data line in that call (the Professional
  behaviour is not documented; the point-in-element tolerance of
  `control_print_dof_line_eps_iso` does not exist for the point family).
