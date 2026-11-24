# control_print_tecplot

## Description

`control_print_tecplot index switch` (manual Professional 6.338) writes
the nodal results in tecplot format at the time steps of
`control_timestep index`. The switch selects the file:

| switch | file (input `excavation.dat`) |
|---|---|
| `-yes` | `excavation.plt`: one zone per print, appended; the first print of a run empties a file left by an earlier run |
| `-separate_index` | `excavation_<index>.plt` with the index in 4 digits (`excavation_0000.plt`, `excavation_0012.plt`), rewritten by every print of that index |
| `-separate_sequential` | one file per print, counter in 4 digits: `excavation_0000.plt`, `excavation_0001.plt`, ... |
| `-mesh` | as `-yes`, with the mesh written again |
| `-no` | nothing |

Each new file starts with its title, variables and the mesh; later zones
of `excavation.plt` share the connectivity of the first one.
`print_tecplot_calculation -yes` appends the final state to
`excavation.plt`.

## Usage

```
control_timestep      0  0.1 0.5
control_print_tecplot 0  -separate_sequential
```

Tests `validation-suite/test-2026/control_print_tecplot_names1`,
`control_print_tecplot_separate_index1`,
`control_print_tecplot_separate_sequential1`.

## Differences with the Professional version

- Names = Professional 02-08-2026 (owner decision 2026-11-02, DIVERGENCES
  `VTK-TECPLOT-NAMES`, class A against the manual, which gives no name
  format for tecplot). Before F5 the GNU appended every print to
  `tecplot.plt`, across runs, and ignored the `-separate_*` switches.
- The zone layout differs (one `ZONE T="t=<time>"` per print against the
  Professional `nodal_data` / `element_average_data` zones with strand
  ids): BACKLOG C27.
