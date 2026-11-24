# dependency_apply — developer

Enum + registration (INTEGER, 1 value, no_index 1, class DEPENDENCY).
Consumed in get_group_data() (group.cc) together with
control_dependency_apply: the current icontrol is read and the control record
overrides the global; -NO short-circuits the whole dependency lookup (the
fallback db(idat,...) read runs directly). Default -yes when neither record
exists.
