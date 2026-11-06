# group_materi_plasti_mpc

## Description

Wall friction reduction for granular material, as
`group_materi_plasti_bounda` (manual Professional 6.687), but the wall is
found through mpc records (6.733). With `group_materi_plasti_mpc index -yes`
an element of group `index` is on a wall when one of its nodes is in an
`mpc_node_number` record. The friction of the element is then reduced by
`group_materi_plasti_mpc_factor` (6.734), default 2/3.

An element that is also on a wall through `group_materi_plasti_bounda`
uses the bounda factor.

## Syntax

```
group_materi_plasti_mpc  index  switch
group_materi_plasti_mpc_factor  index  factor
```

| Parameter | Meaning |
|---|---|
| `index` | Element group of the granular material. |
| `switch` | `-yes` to find walls through mpc records; `-no` (default) ignores them. |
| `factor` | Reduction factor; default 2/3. |

## Example

```
mpc_node_number 0 6 -velx 5 -velx
group_materi_plasti_mpc 0 -yes
group_materi_plasti_mpc_factor 0 0.5
```

Tests: `validation-suite/test-2026/group_materi_plasti_mpc1` to `_mpc4`
(default factor, factor 0.5, `-no`, bounda precedence).

## Differences with the Professional version

- Which elements are on the wall, which factor applies and the reduced
  parameters (phi, phi_flow, c; C114) agree with the Professional
  02-08-2026. The stresses differ as much as the GNU Mohr-Coulomb law
  differs from the Professional without a wall (2e-4 on the test mesh, up
  to about 10 % on other probes).
- The wall detection through `group_materi_plasti_bounda` is broader in
  the GNU (BACKLOG C115); the mpc detection is not affected.
- The mpc records that a `-veln` bounda generates internally also count as
  mpc records (not measured in the Professional).
