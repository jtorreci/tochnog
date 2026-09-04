# group_spring_memory — developer

Enum GROUP_SPRING_MEMORY appended at the end of the item enum (tochnog.h /
tochnog-mod.h) + registration in db_initialize (database.cc): INTEGER, 1
value, class SPRING, data_required GROUP_TYPE. Parse-only: spring.cc does not
consume the memory model (the GNU spring is incremental on the current
configuration = -updated_linear; -total_linear coincides in the 1D corpus
models). Full -total_linear/-updated rotation models PENDING.
