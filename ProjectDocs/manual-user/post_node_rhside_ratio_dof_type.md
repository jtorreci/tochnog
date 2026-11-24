# post_node_rhside_ratio_dof_type

`post_node_rhside_ratio_dof_type dof_type_0 dof_type_1 ...` (manual 6.951)

Restricts the dofs that enter
[post_node_rhside_ratio](post_node_rhside_ratio_method.md) to the dof types you list. For example, the
following line uses only the temperature:

```
post_node_rhside_ratio_dof_type -condif_temperature
```

Default when the record is absent:

- With `materi_velocity` initialised, only the `materi_velocity` dofs are used.
- Otherwise, every principal dof is used.

The GNU name of the same record is `post_node_rhside_ratio_unknowntypes`. See
[names_2026_aliases](names_2026_aliases.md).

Tests: `post_node_rhside_ratio_dof_type1` (default), `post_node_rhside_ratio_dof_type2`
(`-condif_temperature`).
