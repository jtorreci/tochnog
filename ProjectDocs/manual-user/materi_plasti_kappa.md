# materi_plasti_kappa

## What it does

`materi_plasti_kappa` (manual Professional 4.24) is an initialization
option that adds the scalar dof `kap` to the `node_dof` records: the size
of the plastic strain,

```
d(kappa)/dt = sqrt( 0.5 * d(eps_p)/dt : d(eps_p)/dt )
```

(manual theory section "Isotropic hardening and softening"). It is the
isotropic hardening parameter of the plasticity laws: a
`dependency_item ... -kap` / `dependency_diagram` makes a yield parameter
(for example `group_materi_plasti_vonmises`) a function of kappa, and the
Hardening-Soil and Mohr-Coulomb hardening-softening laws read it
directly.

## Why it is useful

It gives isotropic hardening/softening and lets the user plot or target
the accumulated plastic strain.

Since 2026-10-11 (F4 T6f) kappa is kept per INTEGRATION POINT, as the
plastic strain it integrates: each point hardens with its own history.
Before, the GNU version read back the nodal (averaged) kappa interpolated
to the point, so the hardening lagged and was smeared; `validation_3`
(von Mises plate with a hole) had its displacements up to 1.25 % (of the
field maximum) off the Professional (now 0.06 %) and its nodal kappa
11.6 % off (now 2.7 %).
The nodal `kap` written to the output is the mean of the adjacent
elements' values, with the same recovery operator as the strains.

## Input

In the initialization part:

```
materi_velocity
materi_displacement
materi_stress
materi_strain_plasti
materi_plasti_kappa
end_initia
```

| Record | Meaning |
|--------|---------|
| `materi_plasti_kappa` | Adds the `kap` dof (size of the plastic strain) to the node_dof records; per integration point it is the hardening parameter. |

## Example

`validation-suite/test-2026/kappa_ip_state1.dat`: two quad4 (membrane),
von Mises 24.3 hardening to 248.3 at kap 1:

```
group_materi_plasti_vonmises 0 24.3
dependency_item 1 -group_materi_plasti_vonmises 0 -kap 2
dependency_diagram 1 0. 1. 24.3 248.3
```

The top is pushed down 0.01; the bottom reactions equal the Professional
02-08-2026 within 1e-3 (node 1 vely -11.912).

## Differences with the Professional version

- The nodal `kap` of quad4 meshes differs a few per cent from the
  Professional (`validation_3` node 1 0.01600 vs 0.01644) while the
  solution agrees within 0.1 %: the quad4 nodal recovery operator differs
  (BACKLOG C42/C64), as for the nodal stresses.
- The other plasticity state variables (`kapsh`, `pc` of cap1, `rho`,
  `dam`) are still read back from the nodes (BACKLOG C68).
