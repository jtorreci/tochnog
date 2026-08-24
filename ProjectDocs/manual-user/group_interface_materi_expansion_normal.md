# group_interface_materi_expansion_normal

## Description

Thermal strain expansion of the interface in its thickness direction
per unit temperature (manual Professional 6.630). The temperature is
the average of the temperatures of both sides at the interface.

Two effects:

- The MECHANICAL normal strain seen by the gap / tension /
  Mohr-Coulomb criteria is the accumulated strain minus the total
  thermal expansion `alpha * T_avg`.
- The thermal expansion INCREMENT of each step acts as a pseudo-load on
  the normal force (the same incremental pattern as the material
  thermal strains in stress.cc): a temperature rise reduces the
  compressive force / generates tension at constant opening.

Only meaningful with `condif_temperature` (and
`group_interface_materi_memory` `-updated_linear`/`-total_linear`, both
supported by the GNU interface elements; the manual also asks for
`materi_strain_elasti`).

## Usage

```
group_interface_materi_expansion_normal <index> <expansion_coefficient_normal>
```

## Example

```
group_interface 10  -yes
group_interface_materi_elasti_stiffness 10  1000. 1000. 1000.
group_interface_materi_expansion_normal 10  1.e-4
```

Interface heated to T=10 between two fixed blocks: it wants to expand
eps = 1e-4*10 = 1e-3 and, fully constrained, develops the normal force
kn*eps = 1.0 per node pair — the blocks carry the corresponding tension
(test `iface_expansion`: sigxx ~ 2*kn*alpha*T; without the record the
stress is exactly 0).
