# group_materi_plasti_hypo_minimum_void_ratio

## Description

Sets a minimum for the void ratio calculated by the hypoplasticity laws
(manual Professional 6.715). A void ratio that would drop below the
minimum is reset to the minimum. This helps in difficult nonlinear
calculations where the void ratio becomes negative in local zones and
stops the calculation. Check the final results.

It applies to the von Wolffersdorff route (`group_materi_plasti_hypo_wolfersdorff`,
each substep) and to the Masin routes (`group_materi_plasti_hypo_masin*`,
after the integration of the step).

## Syntax

```
group_materi_plasti_hypo_minimum_void_ratio  index  minimum_void_ratio
```

| Parameter | Meaning |
|---|---|
| `index` | Element group. |
| `minimum_void_ratio` | Lowest void ratio allowed. |

## Example

```
group_materi_plasti_hypo_minimum_void_ratio 0 0.599
```

Tests: `validation-suite/test-2026/group_materi_plasti_hypo_minimum_void_ratio1`
(von Wolffersdorff, Professional `other/hypo1` with minimum 0.599) and
`group_materi_plasti_hypo_minimum_void_ratio2` (Masin, `other/hypo7` with
minimum 1.1207).

## Differences with the Professional version

None measured: sigyy and the void ratio agree to 2e-5 relative on both
tests. The SANISAND route (`group_materi_plasti_hypo_sanisand`) ignores
the minimum in the GNU (not measured in the Professional).
