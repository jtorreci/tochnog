# group_materi_plasti_mohr_coul_hardening_softening

## Description

`group_materi_plasti_mohr_coul_hardening_softening` is the hardening /
softening version of the classic Mohr-Coulomb law
(`group_materi_plasti_mohr_coul`). The yield surface is the SAME classic
Mohr-Coulomb surface

```
f = 0.5*(sig_max - sig_min) + 0.5*(sig_max + sig_min)*sin(phi) - c*cos(phi)
```

(where `sig_max`/`sig_min` are the largest/smallest principal stresses,
tension positive), but the cohesion and the friction angles vary LINEARLY
with the accumulated equivalent plastic shear strain `kappa`:

- at `kappa = 0` the law uses `phi_0 c_0 phi_flow_0`;
- at `kappa = kappa_shear_crit` the law uses `phi_1 c_1 phi_flow_1`;
- for `kappa > kappa_shear_crit` the final values are kept constant.

`kappa` is accumulated by the standard plasticity driver as
`kappa = integral sqrt(0.5 * depsp:depsp)` over the plastic strain increment
(the `materi_plasti_kappa` dof). Softening is obtained with
`c_1 < c_0` (or `phi_1 < phi_0`); hardening with `c_1 > c_0`.

## Uso

Place it in the data part, in the element group definition, together with
`materi_plasti_kappa` in the initialization part:

```
materi_plasti_kappa
end_initia
...
group_type 0  -materi
group_materi_elasti_young 0  1000.0
group_materi_elasti_poisson 0  0.0
group_materi_plasti_mohr_coul_hardening_softening 0
                        0.    ( phi_0 = 0 deg )
                        80.   ( c_0 )
                        0.    ( phi_flow_0 = 0, no dilatancy )
                        0.    ( phi_1 = 0 deg )
                        20.   ( c_1: softening )
                        0.    ( phi_flow_1 = 0, no dilatancy )
                        0.5   ( kappa_shear_crit )
```

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_plasti_mohr_coul_hardening_softening` | `phi_0 c_0 phi_flow_0 phi_1 c_1 phi_flow_1 kappa_shear_crit` | Initial (index `_0`) and final (index `_1`) yield/flow friction angles (rad), cohesions, and `kappa` at which the final values are reached. |

- `phi` and `c` enter the YIELD surface; `phi_flow` enters the FLOW
  direction (associated flow when `phi_flow = phi`, non-associated
  otherwise).
- The interpolation is linear in `kappa/kappa_shear_crit`, clamped to
  `[0, 1]`.
- The accumulated `kappa` is the `materi_plasti_kappa` node dof; the
  standard plastic driver (cutting plane, `plasti_kappa`, boundary
  reduction) performs the integration.
- For softening runs, set `control_timestep_iterations` (e.g. 8) so the
  plastic return can chase the downward-moving surface; with the default
  single global iteration the stress stays above the current surface
  (weak kappa coupling). Use several small time steps so `kappa` grows
  gradually instead of jumping past `kappa_shear_crit` in one step.

## Validation

Uniaxial tension rig with `phi = 0` (Tresca-like, smooth path): the tensile
strength is `sig_t = 2*c(kappa)`. With `c_0 = 80`, `c_1 = 20`,
`kappa_shear_crit = 0.5` and a 20-step rig (`dt = 0.02`, total strain 0.4),
the measured `kappa` is 0.318 -> ratio 0.636 -> `c = 41.8` -> analytic
`sig_t = 83.6`; the measured `sigxx = 81.9` (98%). Test `mmchs_soft` in
`scripts/build_safe.sh`.
