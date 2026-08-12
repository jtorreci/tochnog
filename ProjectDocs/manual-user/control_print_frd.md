# control_print_frd

## Description

`control_print_frd` writes the mesh and the nodal results in the
**CalculiX result format** (`.frd` file). The `.frd` files can be plotted
by CGX (the CalculiX postprocessor), FreeCAD and prepomax.

This is the CalculiX counterpart of `control_print_gid` /
`control_print_vtk` / `control_print_gmsh`: use it when the results are to
be viewed or post-processed in the CalculiX ecosystem or FreeCAD.

Only results for 2D and 3D isoparametric elements are written. For
structural elements (trusses, beams, ...) nothing is printed.

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_timestep` record whose results you want to export:

```
control_print_frd 20 -yes
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `20`      | Index of the control record. Must match an active `control_timestep`. |
| switch    | `-yes`: single file `<base>.frd`, mesh written only the first time, results appended per time step. `-separate_index`: file `<base><icontrol>.frd`. `-separate_sequential`: files `<base>0.frd`, `<base>1.frd`, ... one per call. |

## Options

- `control_print_frd_freecad index switch` — write result names suited
  for FreeCAD: `DISP`, `STRESS`, `TOSTRAIN`, `NDTEMP`. Default `-yes`
  behaviour is embedded in `control_print_frd` (the standard names are
  always used).
- `control_print_frd_prepomax index switch` — write result names suited
  for prepomax.

Result names written:

- `DISP` — for `materi_displacement` or `materi_velocity_integrated`
- `VELO` — for `materi_velocity`
- `STRESS` — for `materi_stress`
- `TOSTRAIN` — for `materi_strain_total`
- `NDTEMP` — for `condif_temperature`
- other dofs use the Tochnog names truncated to 8 characters

## Output

File structure (CalculiX FRD format):

- Model header `1C` + user record `1U...`
- Nodal coordinates block `2C` (`-1` node lines, `-3` terminator)
- Element definition block `3C` (`-1` element, `-2` connectivity)
- Per time step: `1PSTEP` + `100CL` headers, then one result block per
  dof with `-4` (dataset) / `-5` (component) records and `-1` nodal
  values, terminated by `-3`.
- `9999` end marker is written by CGX/CalculiX; Tochnog writes `-3`.

## Example

```
control_print_frd                       20  -yes
control_timestep                        20  0.001 0.04
control_print                           20  -time_current
```

Produces `<base>.frd`. Open it in CGX (`cgx`), FreeCAD or read with
CalculiX tools.
