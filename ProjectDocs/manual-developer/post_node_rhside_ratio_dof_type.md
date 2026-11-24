# post_node_rhside_ratio_dof_type — developer

The record is an alias of `POST_NODE_RHSIDE_RATIO_UNKNOWNTYPES` (`database.cc`, see
[names_2026_aliases](names_2026_aliases.md)).

`post_node_rhside_fixed_free()` in `post.cc` reads it. A principal dof whose `DOF_TYPE` is not in the
list does not enter the ratio. The signed sums `post_node_rhside_fixed` and `post_node_rhside_free`
are not filtered.

When the record is absent and `materi_velocity` is initialised, the list becomes
`{-MATERI_VELOCITY}`. This is the default of manual 6.951, measured on the Professional 02-08-2026:

| Input | Ratio |
|---|---|
| No record | 2.56e-15 |
| `-materi_velocity` | 2.56e-15 |
| `-condif_temperature` | 0.980198 |

The input is an elastic velocity field plus one iteration with a temperature-dependent
conductivity. Before this change, the GNU version used every principal dof without the record, which
gave 0.980198 for this input.
