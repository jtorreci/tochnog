# group_materi_membrane

## Description

Plane stress (manual Professional 6.685). With `-yes` the zz stress
becomes zero in 2D, and the yy and zz stresses become zero in 1D. Without
the record the 2D model is plane strain. Useful for thin plates loaded in
their own plane.

The zz strain `eptzz` is a real result: the total zz strain of the final
state, that is the elastic part (`-nu (sigxx + sigyy) / E` in 2D) plus,
for a plastic material, the accumulated plastic part `eppzz`.

## Usage

```
group_materi_membrane index switch
```

| parameter | meaning |
|---|---|
| `index` | element group |
| `switch` | `-yes` plane stress; `-no` (default) plane strain |

Cannot be combined with `group_materi_elasti_compressibility`,
`group_materi_hyperelasticity`, `group_materi_viscosity` or
`group_axisymmetric` (1D).

## Example

Single quad4, von Mises with hardening, plane stress, stretched to
eptxx 0.02 (`validation-suite/test-2026/membrane_eptzz_total1.dat`):

```
group_materi_elasti_young   0  7.e3
group_materi_elasti_poisson 0  0.2
group_materi_memory 0  -total_linear
group_materi_membrane 0 -yes
group_materi_plasti_vonmises 0  24.3
```

Analytic: sigxx 27.42, eptzz = -0.2 * 27.42 / 7000 - 0.0080414 =
-0.0088249. GNU: -0.00863 at dt 0.01, -0.00880 at dt 0.001 (the plastic
strain split between yy and zz converges with the time step).

## Differences with the Professional version

- **eptzz of a plastic membrane (MEMBRANE-EPTZZ, deliberate, project
  decision 2026-10-12).** The GNU stores the total zz strain (elastic +
  accumulated plastic). The Professional 02-08-2026 stores only the zz
  correction of the last time step, which depends on the time step: in the
  example above -9.6e-5 at dt 0.01 and -7.9e-4 at dt 0.1, while the analytic
  value is -0.0088. Guard tests `test-2026/membrane_eptzz_total1` (dt 0.01)
  and `membrane_eptzz_total2` (dt 0.1); the Professional fails both. See
  `ProjectDocs/DIVERGENCES-PROFESSIONAL.md` (MEMBRANE-EPTZZ).
- Known GNU limitation (not a Professional difference): with plasticity the
  split of the plastic strain between yy and zz is first order in the time
  step (eptzz -0.00721 at dt 0.1 in the example). Use smaller steps when
  eptzz matters. BACKLOG C69.
