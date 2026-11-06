# materi_plasti_hypo_substepping

## Description

The von Wolffersdorff hypoplasticity integrates each timestep in substeps.
A substep is halved while its stress increment exceeds 2 % of the current
stress. `materi_plasti_hypo_substepping -no` (manual Professional 6.803)
switches this accuracy substepping off, which is faster but less
accurate. `control_materi_plasti_hypo_substepping index switch` (6.149)
does the same for the timesteps of one `control_timestep` index and wins
over the global record. The default is `-yes`.

A substep in which the hypoplastic stress update fails is still halved
with `-no`.

## Syntax

```
materi_plasti_hypo_substepping  switch
control_materi_plasti_hypo_substepping  index  switch
```

| Parameter | Meaning |
|---|---|
| `switch` | `-yes` (default) or `-no`. |
| `index` | `control_timestep` index (control record only). |

## Example

```
materi_plasti_hypo_substepping -no
control_materi_plasti_hypo_substepping 20 -yes
```

Tests: `validation-suite/test-2026/materi_plasti_hypo_substepping1`
(biaxial test of the Professional `other/hypo1`: sigyy -863.716 instead of
-862.926) and `materi_plasti_hypo_substepping2` (the control record wins).

## Differences with the Professional version

- The manual header shows `materi_plasti_hypo_substepping index switch`;
  the Professional 02-08-2026 rejects the index and accepts
  `materi_plasti_hypo_substepping switch`, as the GNU.
- The Masin routes (`group_materi_plasti_hypo_masin*`) integrate with a
  fixed number of Runge-Kutta substeps in the GNU, so the switch has no
  effect there. In the Professional it moves the `other/hypo7` sigyy from
  -231.8090 to -231.8126 (DIVERGENCES HYPO-SUBSTEPPING-MASIN).
