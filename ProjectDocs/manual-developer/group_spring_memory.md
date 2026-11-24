# group_spring_memory — developer

Enum GROUP_SPRING_MEMORY appended at the end of the item enum (tochnog.h /
tochnog-mod.h) + registration in db_initialize (database.cc): INTEGER, 1
value, class SPRING, data_required GROUP_TYPE. Parse-only: spring.cc does not
consume the memory model. Since C88 (F4 T6f, 2026-10-19) the spring
elongation is the relative displacement projected on a fixed direction
(`group_spring_direction`, else the initial axis), which is what the
Professional 02-08-2026 gives for every memory value (the 90-degree rotation
probe carries no force with `-total_linear`, `-updated_linear` and
`-updated`). Before C88 the comment here said the GNU spring was
"incremental on the current configuration": the length was the node
distance, which is what C88 removed. See
[group_spring_direction](group_spring_direction.md).
