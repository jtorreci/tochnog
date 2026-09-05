# group_materi_undrained_capacity

Capacity `C` of the UNDRAINED groundwater analysis of a soil/water model
(manual Professional 6.760 + theory 2.2.7 "Undrained groundwater analysis").

```
group_materi_undrained_capacity <element_group> <C>
```

When the element group carries the capacity `C` (and the undrained analysis
is not switched off), the total groundwater pressure change of the element
follows from the groundflow storage equation **without permeability**:

```
C * p_dot = div(v_material) = eps_vol_dot
```

which is solved on the ELEMENT level, so the groundwater pressure does NOT
need to be a dof of the system matrix. This is the tool of choice for
UNDRAINED (u-p) soil analyses when you do not want the large, ill-conditioned
coupled soil-groundwater matrix: the displacement/velocity equations are
solved alone and the pore pressure response is computed element-wise from the
volumetric strain rate.

## What the capacity does

Per step and per integration point the excessive undrained pressure is
accumulated as `p_u += d(eps_vol)/C` (eps_vol compression negative). The
pressure acts isotropically on the skeleton: the TOTAL stress used in the
momentum equilibrium is

```
sigma_total = sigma_effective + (p_fixed + p_u) * I
```

where `p_fixed` is the groundflow total pressure of the hydraulic heads when
a groundflow pressure field is present (manual 2.2.7: "the fixed total
pressure from the hydraulic pressure heads plus the excessive undrained
pressure of the remainder of the calculation as the full total pressure").
The EFFECTIVE stress result (`-sigxx`, ... dofs) keeps the drained
constitutive value, so the undrained response is the drained one plus the
volumetric stiffness `1/C` (for an isotropic linear soil the undrained bulk
modulus is `K_undrained = K_drained + 1/C`).

The pressure of each element/step is stored in
`element_intpnt_materi_undrained_pressure` (one value per integration point;
`..._average` holds the element mean).

## Control records

- `control_materi_undrained_apply <index> -yes/-no` (6.153): switches the
  undrained analysis on/off for the timestep with the same index. Default
  `-yes`.
- The undrained analysis is active when the element group carries
  `group_materi_undrained_capacity`. A group with the capacity but a
  `-groundflow` group type keeps the usual coupled flow when the flow dofs
  are solved; the capacity adds the element-level excessive pressure.

## Typical strategy (manual 2.2.7)

Set the fixed hydraulic pressure heads once, then solve the material steps
with the undrained machinery:

```
groundflow_phreatic_level <level>          (optional)
control_reset_dof  ... -pres              (set the initial heads)
bounda_dof   40  -all -tpres              (fix the total pressure)
bounda_time  40  <p_total>
group_type   0  -materi -groundflow
group_materi_elasti_young 0 <E>
group_materi_undrained_capacity 0 <C>
...
control_timestep 30 ...
control_materi_undrained_apply 30 -yes
```

## Example

undrained2 of the corpus (1D bar, E = 3, C = 1, unit force): the drained
response `sigma = E*eps` and the excessive pressure `p = eps/C` share the
force, so `eps = F/(E + 1/C) = -0.25`, `sigxx = -0.75` and the undrained
pressure is `-0.25`.

```
element 1  -bar2 1 2
group_materi_elasti_young 0  3.
group_materi_undrained_capacity 0  1.0
bounda_dof   10  1 -velx
bounda_force 20  2 -velx
bounda_time  20 -1.0
control_timestep 10  1. 1.0
```

## Related items

- `element_intpnt_materi_undrained_pressure` (6.441): output record with the
  undrained total pressure per element and integration point (target_item /
  control_print item).
- `control_materi_undrained_apply` (6.153): per-step switch.
- `-topres`/`-tpres` in `bounda_dof` (2.4.1): prescribe the fixed total
  pressure (see bounda_dof.md).
- `post_calcul -materi_stress -young_apparent/-poisson_apparent` (6.903):
  apparent E and Poisson ratio of the last time step (used by ground17 of
  the corpus to check the effective response of the undrained skeleton).
