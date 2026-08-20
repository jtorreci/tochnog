# group_groundflow_consolidation_apply

## Description

Element-group level switch that enables or disables the consolidation term (the
material divergence part of the groundflow equation) for the elements of the
group.

- `-no`: consolidation is not applied for the elements of the group.
- `-yes`: consolidation is active.

The record is indexed by `element_group`. The global
`groundflow_consolidation_apply` switch applies everywhere else, and
`control_groundflow_consolidation_apply` overrides per timestep.

## Usage

```
group_groundflow_consolidation_apply <element_group> <switch>
```

## Parameters

| Parameter       | Meaning                                                |
|-----------------|--------------------------------------------------------|
| `element_group` | Element group to which the switch applies.             |
| `switch`        | `-yes` or `-no`.                                       |

## Example

```
group_type 0  -groundflow -materi
group_groundflow_consolidation_apply 0  -no
```

The material divergence term is skipped for the elements of group 0, e.g. when
that group represents a region where the skeleton is rigid (no consolidation).
