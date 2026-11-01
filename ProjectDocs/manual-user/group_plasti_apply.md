# group_plasti_apply

## Description

`group_plasti_apply index switch` (manual Professional 6.762): with `-no`,
any plasticity data of element group `index` is neglected; the other groups
keep their plasticity. Default `-yes`.

## Usage

```
group_plasti_apply index switch
```

| Parameter | Meaning |
|-----------|---------|
| `index` | the element group (see `element_group`) |
| `switch` | `-yes` (default) or `-no` (the group's plasticity is neglected) |

What is neglected is the same list as for `plasti_apply -no` (yield-surface
laws, direct cut-offs, hypoplasticity, incremental models, visco-plasticity),
for that group only. A group `-yes` does not overrule a global
`plasti_apply -no` or `control_plasti_apply -no`.

## Example

```
group_materi_plasti_tension 0 1.
group_materi_plasti_tension 1 1.
group_plasti_apply 1 -no
```

Self-checking inputs in `validation-suite/test-2026/`:
`group_plasti_apply1.dat` (two bars: group 0 yields at 1, group 1 elastic
3) and `group_plasti_apply2.dat` (group -yes with `plasti_apply -no`:
elastic). Both equal the Professional 02-08-2026.

## Differences with the Professional version

None measured. Truss, spring, interface and beam plasticity are not gated
(as `plasti_apply`).
