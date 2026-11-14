# groundflow_phreatic_level

Defines the groundwater level (phreatic line) of a groundflow analysis
(manual Professional 2.4.1): the level where the pore water meets the free
air. The record holds either one value (a horizontal level) or a table of
`x y` (2D) / `x y z` (3D) points describing a non-horizontal phreatic line.

```
groundflow_phreatic_level <water_level>          (horizontal)
groundflow_phreatic_level <x0> <y0> <x1> <y1>    (2D table)
```

## Effect on the calculation

- The static groundwater pressure used by the post-processing split
  (`post_calcul -groundflow_pressure -static_pressure`) is
  `p_static = rho*g*(level - coord_vertical)` for the nodes covered by the
  level; the total pressure is `p_total = pres_dof + p_static` and the
  dynamic pressure `p_dynamic = p_total - p_static = pres_dof`.
- Free surface / dry zone: nodes at or above the phreatic line are dry and
  carry no water pressure. The GNU imposes `p_total = 0` there by bounding
  the pore pressure dof to 0 (`pres_dof = 0`; the static part is clamped to
  the atmospheric pressure above the level). This confines the saturated
  flow domain below the level and reproduces the Professional behavior
  (its hydraulic head on/above the phreatic line is `h = rho*g*level`, i.e.
  `p_dynamic = 0`). Explicit pressure boundary conditions
  (`bounda_dof -pres` / `-topres`) win over this default.
- `groundflow_phreatic_bounda -yes` / `-no` (manual 6.570): switches this
  free-surface condition on (`-yes`, the default, also when the record is
  absent, as in the Professional) or off (`-no`: only the explicit
  pressure boundary conditions hold). The GNU name
  `groundflow_phreaticlevel_bounda` also takes the older `-method1` (bound
  only the nodes on the level) and `-method2` (on and above the level),
  which prescribe `pres_dof = rho*g*level` there.
- Use `groundflow_phreatic_level_multiple` when different parts of the
  domain have different levels.

## Differences with the Professional version

- With `groundflow_phreatic_bounda -no`, the Professional reports a total
  pressure of 0 at every node above the level; the GNU reports
  `pres_dof + rho*g*(level - z)` capped by `groundflow_total_pressure_limit`
  (for example -0.5 instead of 0 in `groundflow_phreatic_bounda2`). The
  heads (nodal dofs) are the same (DIVERGENCES PHREATIC-BOUNDA-NO-POST).

## Example

```
groundflow_density 1.
force_gravity 0. -10.
groundflow_phreatic_level 0.25
bounda_dof 10 -lower_edge -topres
bounda_time 10 -10.
```

Soil column with the water table at 0.25 m and the total pore pressure -10 at
the bottom: the saturated zone (below 0.25) carries the prescribed flow
state, the zone above the level stays dry (p_total = 0). Verified against
the Professional on ground14/15/16 of the corpus.
