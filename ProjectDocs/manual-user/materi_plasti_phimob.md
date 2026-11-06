# materi_plasti_phimob

## Description

Initialization record (manual Professional 4.26). Adds the mobilized
friction angle phi_mob, in degrees, to the `node_dof` records (label
`-phimob`). phi_mob is the angle for which the Mohr-Coulomb yield function

    f = 0.5 (s2 - s0) + 0.5 (s2 + s0) sin(phi_mob) - c cos(phi_mob)

becomes zero, with s0 / s2 the smallest / largest principal stress
(tension positive) and c the cohesion of the group. At yield with small
cohesion phi_mob reaches the friction angle; with substantial cohesion it
can exceed it. Negative values mean the stress state is far from yield.

Available for `group_materi_plasti_mohr_coul` and
`group_materi_plasti_mohr_coul_direct`; other groups get 0. Requires
`materi_stress`.

## Syntax

```
materi_stress
materi_plasti_phimob
end_initia
```

## Example

```
target_item 0 -node_dof 4 -phimob
```

Tests: `validation-suite/test-2026/materi_plasti_phimob1` (at yield: phi
0.5 rad = 28.6479), `_2` (elastic, c 1: 5.39561), `_3`
(`mohr_coul_direct`), `_4` (Drucker-Prager: 0).

## Differences with the Professional version

- The manual says "in radians" and "in degrees"; the Professional
  02-08-2026 writes degrees, as the GNU.
- When f has two roots in [-90, 90] degrees the GNU takes the one closest
  to 0, as the Professional on the measured cases (5.3956 instead of 90);
  90 when f has no root.
- The GNU `post_calcul` operator `-phimob` (not accepted by the
  Professional) keeps its own definition.
