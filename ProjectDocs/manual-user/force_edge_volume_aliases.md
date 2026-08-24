# force_edge / force_volume (Professional names) + control_materi gates

## Professional keyword names via prefix translation

The Professional manual names `force_edge_*` and `force_volume_*` are
accepted directly: the input parser resolves them inside `db_number` to
their GNU equivalents `force_element_edge_*` / `force_element_volume_*`
(same physics, same record layout). This works everywhere — keyword
detection AND the end-of-values detection of variable-length records
(both call db_number).

## New restriction variants

For the three edge families (`force_edge`, `force_edge_normal`,
`force_edge_water`) and `force_volume`:

- `_element index element_0 ...` — restrict to the listed elements.
- `_element_group index group_0 ...` — restrict to elements of the
  listed element groups.
- `_element_side index element side ...` — pairs element/side (edge
  families).
- `_node index node_0 ...` — only the listed global nodes of the edge.
- `_element_node index element node_0 ...` — element + local node
  numbers.
- `_node_factor index element f_0 f_1 ...` (edge/normal) — per-local-node
  multiplication factors (first value = element number).
- `force_edge_water_factor` — spatial polynomial factor for the water
  edge load (new; the GNU had none).

## control_materi_*_apply per-timestep gates

- `control_materi_viscosity_apply index -no` — ignore any viscosity in
  the input file for these timesteps.
- `control_materi_plasti_tension_apply index -no` — ignore any
  tension-plasticity data (the Mohr-Coulomb direct cutoff stays).
- `control_materi_plasti_visco_apply index -no` — ignore any
  visco-plasticity data.
- `control_materi_damage_apply` / `control_materi_failure_apply index
  -no` — ignore any damage data (the damage stays at its old value).
- `control_materi_updated_apply index -no` — any -updated material
  memory is set to -updated_linear for these timesteps. (The -yes
  direction — defaulting UNSPECIFIED memory to -updated — is not wired:
  the GNU cannot distinguish unspecified from explicitly-set after
  input; documented difference.)

Registered and parsed, behaviour NOT wired (documented partials):
`control_materi_dynamic`, `control_materi_elasti_k0`,
`control_materi_plasti_hardsoil_gammap_initial`,
`_hypo_niemunis_visco_ocr_apply`, `_hypo_pressure_dependent_void_ratio`,
`_hypo_substepping`, `control_materi_undrained_apply`.

## Examples

```
force_edge_normal_geometry 0  -geometry_line 1
force_edge_normal 0  5.0
```
(test `fedge_alias`, Professional syntax only; sigxx = 5.0)

```
force_volume 0  1.0 0.
force_volume_element_group 0  1
```
(test `fvol_elem`: body force only in group 1; each fixed support takes
half -> sigxx 0.5, the unforced column 0)

```
control_materi_plasti_tension_apply 1  -no
```
(test `cmat_gate`: the tension cutoff is ignored during timestep 1 ->
linear sigxy -88.4 instead of the capped ~1)
