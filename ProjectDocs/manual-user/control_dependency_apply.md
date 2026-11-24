# control_dependency_apply

Per-timestep switch of the dependency machinery (manual Professional 6.126):
`-yes` includes the `dependency_item`/`dependency_diagram` dependencies for the
control steps with the same index, `-no` excludes them. Precedence: the control
record overrides the global dependency_apply for its index.

Consumed in `get_group_data()` (group.cc): when disabled, the dependency
diagram lookup is skipped and the group data is read directly. ground11 of the
corpus disables the dependency during the first control block.

PENDING: the dependency monitor labels are the DOF labels of the model; the
Professional groundflow tests monitor the post item `-to_pres` (total pressure),
which is not a dof of the GNU groundflow formulation.
