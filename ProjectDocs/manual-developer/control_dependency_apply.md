# control_dependency_apply — developer

Enum + registration (INTEGER, 1 value, class CONTROL; NO data_required — the
same-index combination check would otherwise require a dependency_item at the
control index, which ground11 does not have: dependency_item is at index 10,
the gate at index 100). Gate evaluated inside get_group_data() (group.cc);
see dependency_apply. PENDING: monitoring a non-dof item (-to_pres) — the
Professional groundflow tests monitor the total pressure post item, which is
not a DOF label of the GNU model.
