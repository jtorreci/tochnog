# print_frd_freecad

## Description

`print_frd_freecad switch` (manual Professional 6.979) is the default of
`control_print_frd_freecad index switch` (6.289) for every
`control_print_frd` file; a control record overrules it for its index.
Default `-no`. With `-yes` the frd file is meant for FreeCAD, which only
reads files with displacements, stresses and total strains: when
`materi_displacement` (or `materi_velocity_integrated`), `materi_stress`
and `materi_strain_total` are not all initialised, no frd file is written
(manual 6.288). Otherwise the file is the normal frd file (the result names
DISP, STRESS, TOSTRAIN, NDTEMP are always used, see
[control_print_frd](control_print_frd.md)).

## Usage

```
control_print_frd 0 -yes
print_frd_freecad -yes
```

| Parameter | Meaning |
|---|---|
| switch | `-yes`: FreeCAD file (needs disp + stress + total strain). `-no` (default). |

Tests `validation-suite/test-2026/print_frd_freecad1/2`,
`control_print_frd_freecad1`.

## Differences with the Professional version

None measured: the Professional 02-08-2026 writes no file without
`materi_strain_total` and an unchanged file with it. Not implemented: the
`PEEQ` name for `post_calcul -groundflow_pressure -total_pressure` (the
GNU frd writer prints no post_calcul fields; BACKLOG C27).
