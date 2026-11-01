# group_materi_plasti_mohr_coul

## Description

Classic Mohr-Coulomb plasticity (manual Professional 6.725). The yield
function is

    f = 0.5 (sig1 - sig3) + 0.5 (sig1 + sig3) sin(phi) - c cos(phi)

with sig1 the largest and sig3 the smallest principal stress (tension
positive). The flow rule uses the same expression with `phi_flow` instead
of `phi` (non-associative flow when `phi_flow` differs from `phi`). It is
useful for soils and rocks where the strength grows with the confining
pressure.

## Usage

```
group_materi_plasti_mohr_coul index phi c phi_flow
```

| parameter | meaning |
|---|---|
| `index` | element group |
| `phi` | friction angle, radians |
| `c` | cohesion |
| `phi_flow` | dilatancy angle of the flow rule, radians |

The uniaxial tensile strength (no tension cut-off) is
`2 c cos(phi) / (1 + sin(phi))`. Hardening/softening of `c` and `phi` with
kappa: see `group_materi_plasti_mohr_coul_hardening_softening`.

## Example

Uniaxial tension of one quad4, plane strain, `poisson` 0
(`validation-suite/test-2014/mmc_tension.dat`): sigxx reaches
`2 * 40 * cos(30 deg) / (1 + sin(30 deg)) = 46.188`.

```
group_materi_elasti_young  0  1000
group_materi_elasti_poisson  0  0
group_materi_plasti_mohr_coul  0  0.5236 40 0
```

## Differences with the Professional version

- **Vertex stress state in uniaxial tension (MC-VERTEX-ZZ, deliberate,
  project decision 2026-10-12).** In `mmc_tension` (plane strain, `poisson`
  0, no tension cut-off) the stress reaches the Mohr-Coulomb vertex with
  sigyy = sigzz = 0. Any plastic eppzz > 0 would move the stress off that
  vertex, so the continuum solution has sigzz = 0 and eppzz = 0. The GNU
  gives exactly that for every time step (sigzz ~1e-9). The Professional
  02-08-2026 leaves a time-step dependent artefact: sigzz 0.974 and eppzz
  -9.7e-4 at dt 0.05, 0.367 at dt 0.047, 1.44 at dt 0.0047, and 0 at dt 0.1.
  Both programs give the same sigxx 46.19. Guard tests:
  `test-2026/mmc_tension_vertex1` (dt 0.05) and `mmc_tension_vertex2`
  (dt 0.047). The Professional fails both. See
  `ProjectDocs/DIVERGENCES-PROFESSIONAL.md` (MC-VERTEX-ZZ).
