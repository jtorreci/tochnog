# post_calcul_safety_method

Determines how the hydraulic safety factors of `post_calcul -materi_stress
-safety_piping` / `-safety_lifting` are computed (manual Professional 6.919):
`-vertical` (default, the vertical normal stress; one value), `-prival` (the
three principal stresses; three values `..._prival_0..2`) or `-global` (the
three global normal stresses; three values `..._global_x/y/z`).

## Safety factors

`post_calcul -materi_stress -safety_lifting` computes, for the stress selected
by the method, `(sigma_i + p_total)/p_total`, and
`post_calcul -materi_stress -safety_piping` computes `(sigma_i + p_dynamic)/
p_dynamic`, with the tochnog sign conventions (compression negative, water
pressures negative):

- `p_total` is the total pore pressure of the groundflow split
  (`post_calcul -groundflow_pressure -total_pressure`); with a phreatic level
  or a static height it equals `pres_dof + p_static`, otherwise
  `pres_dof - rho*g*z`.
- `p_dynamic = p_total - p_static`, the excess over the hydrostatic state of
  the level (`p_static = rho*g*(level - coord_vertical)`).

`sigma_i` is the vertical normal stress (sigma_yy in 2D, sigma_zz in 3D,
sigma_xx in 1D), one of the three principal stresses listed from the MOST
compressive (`-prival`: `prival_0` is the minimum), or the global normal
stress sigma_xx/sigma_yy/sigma_zz (`-global`: slots `_x`, `_y`, `_z`).

When the pressure in the denominator is zero the factor is set to 0
(`post_calcul_safety_default` eps very small, value 0). The optional
`post_calcul_safety_maximum` caps the factor.

## Usage

```
post_calcul_safety_method -vertical | -prival | -global
post_calcul -materi_stress -safety_piping -materi_stress -safety_lifting
post_calcul_safety_maximum <value>
```

The generated items are named (post_calcul_label of the .dbs):
`safety_piping`/`safety_lifting` (vertical, one value),
`safety_piping_prival_0..2`/`safety_lifting_prival_0..2`,
`safety_piping_global_x/y/z`/`safety_lifting_global_x/y/z`, and can be used
as `target_item` labels of `node_dof_calcul`/`post_point_dof_calcul`.

## Example

```
groundflow_phreatic_level 0.25
groundflow_density 1.
force_gravity 0. -10.
bounda_dof 10 -lower_edge -topres
bounda_time 10 -10.
...
post_calcul_safety_method -prival
post_calcul -groundflow_pressure -total_pressure -materi_stress -safety_piping
            -materi_stress -safety_lifting
post_point 10 0.05 0.0
...
target_item 10 -post_point_dof_calcul 10 -safety_piping_prival_0
target_value 10 2.3333 1.e-3
target_item 20 -post_point_dof_calcul 10 -safety_lifting_prival_0
target_value 20 2. 1.e-3
```

Verified against the Professional .dbs on ground15 (prival) and ground16
(global): at the bottom point of the piping/lifting benchmark the factors are
2.3333 (piping) and 2.0 (lifting) with `sigma_i = p_total = -10` and
`p_dynamic = -7.5`.
