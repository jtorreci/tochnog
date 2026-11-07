# group_materi_plasti_hypo_void_ratio_linear

## Description

How the void ratio of the von Wolffersdorff (and low angles) hypoplastic
law changes with the deformation (manual Professional 6.722).

- `-yes` (default): linear approximation from the trace of the strain
  increment, `e_new = e_old + (1 + e_old) tr(d_eps)`. Hypoplastic
  parameters are normally calibrated with this approximation.
- `-no`: exact volume change from the determinant of the deformation,
  `e_new = (1 + e_old) det(I + d_eps) - 1`.

## Syntax

```
group_materi_plasti_hypo_void_ratio_linear index -yes|-no
```

## Example

```
group_materi_plasti_hypo_wolffersdorff 0 33. 1000. 0.25 0.95 0.55 1.05 0.25 1.
group_materi_plasti_hypo_void_ratio_linear 0 -no
```

Tests: `validation-suite/test-2026/group_materi_plasti_hypo_void_ratio_linear1`, `2`
(Professional `other/hypo11` without the Niemunis record).

## Differences with the Professional version

- None measured: `other/hypo11` sigyy -1.0198992 / -1.0199336 (GNU) against
  -1.019954 / -1.019924 (Professional 02-08-2026) for `-yes` / `-no`.
- Before F5 T3 D3 the GNU read this name as an alias of
  `group_materi_plasti_hypo_pressuredependentvoidratio`, so `-yes` switched
  on the pressure dependent initial void ratio correction (BACKLOG C118).
- The Masin routes keep their own void ratio update.
