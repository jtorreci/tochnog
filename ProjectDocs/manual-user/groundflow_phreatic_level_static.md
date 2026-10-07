# groundflow_phreatic_level_static

Static groundwater from a phreatic level, without solving the groundwater
flow (manual Professional 6.574). The pore pressures below the
`groundflow_phreatic_level` are set to the hydrostatic (static) pressure of
that level. You do not need the `groundflow_pressure` unknown, its initial
values or the hydraulic heads.

This is the quick way to get the initial effective stresses of a soil
under a water table, or of a soil under free water (the level may lie
above the mesh).

```
groundflow_phreatic_level_static -yes
```

## What it does

For the element groups whose `group_type` holds `-groundflow`:

- **Pore pressure on the soil skeleton.** At every integration point the
  total (pore) pressure is `p = rho_w * g * (z_level - z)`, where `rho_w`
  is `groundflow_density` and `g` is the vertical `force_gravity`. Above
  the level the pressure is the atmospheric pressure (0 by default). The
  pressure acts on the skeleton like the pressure of a groundflow
  analysis. The stress dofs hold the effective stress.
  `groundflow_pressure_factor` scales it, as in a groundflow analysis.
- **Wet or dry density.** With `group_materi_density_groundflow
  rho_wet rho_dry`, an element takes `rho_wet` when the average of the
  static total pressures of its nodes is negative (below the level), and
  `rho_dry` otherwise.
- **Output.** `post_calcul -groundflow_pressure -total_pressure`
  (`-static_pressure`, `-dynamic_pressure`) gives the static pressures.

Groups without `-groundflow` get no pore pressure and keep the wet
density everywhere (a wall that is impermeable, for example).

The switch has no effect when `groundflow_pressure` is initialized: the
normal groundflow analysis then applies.

## Parameters

| value | meaning |
|---|---|
| `-yes` | static pore pressures from the phreatic level |
| `-no` | default: no static pore pressures |

## Example

A 1 x 100 soil column, water table at y = -50, dry density 1.5 above the
table and wet density 2.0 below it:

```
groundflow_density 1.
force_gravity 0. -10.
groundflow_phreatic_level   -50.
groundflow_phreatic_level_static   -yes

group_type 0  -materi -groundflow
group_materi_density_groundflow 0  2.0 1.5
group_materi_elasti_young 0  1.0e6
group_materi_memory 0  -total_linear

post_calcul -groundflow_pressure -total_pressure
```

The effective vertical stress at the bottom is
`-1.5*10*50 - (2.0-1.0)*10*50 = -1250`. The pore pressure there is -500.
Full inputs: `validation-suite/test-2026/phreatic_static1.dat` and
`phreatic_static_submerged1.dat` (level 10 m above the top of the column,
with `force_edge_water` on the top).

## Differences with the Professional version

The Professional 02-08-2026 applies the static pore pressure to the soil
skeleton **only when** the input also holds
`post_calcul -groundflow_pressure -total_pressure`. Without that
post-processing record it keeps the wet/dry density choice but applies no
pore pressure (in the example above: -1745 instead of -1250). The GNU
always applies the pore pressure, as the manual describes. With the
`post_calcul` record, both give the same results. Test `excavate1` of the
Professional has no such record, so its target differs in the GNU.

- **`post_calcul -groundflow_pressure -total_pressure` with a pressure dof**
  (GNU extension / input dialect; AUDIT 1.1; BACKLOG C10). GNU: Not
  required. Professional 02-08-2026: Required when a pressure dof exists
  (and couples the static pore pressure, C10). Status: Input dialect. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
