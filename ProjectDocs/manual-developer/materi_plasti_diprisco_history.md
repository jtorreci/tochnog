# materi_plasti_diprisco_history

## Implementación

- **Initia**: the `materi_history_variables` branch in `input.cc` now
  accepts BOTH strings (`materi_history_variables` and
  `materi_plasti_diprisco_history`, manual Professional 4.18): same
  parser branch, same `materi_history_variables` value, same shared
  `hisv_indx` dof (`dof_type = -MATERI_HISTORY_VARIABLES`, basenames
  `hisv0 .. hisvN-1`). Using the alias sets the flag
  `materi_plasti_diprisco_history` (initia.cc).
- **Enum**: `MATERI_PLASTI_DIPRISCO_HISTORY` in `tochnog.h` /
  `tochnog-mod.h`; name registered in `database.cc`.
- **check.cc**: the alias self-check requires `materi_stress` (the
  history dof is updated inside `set_stress`). The
  `GROUP_MATERI_PLASTI_DIPRISCO` / `_RT` group checks now accept
  EITHER initia name (`check_unknown_atleastone`), so a Professional
  input with the alias passes the checks.
- The di Prisco law (plasti.cc) is untouched: it reads
  `old_hisv`/`new_hisv` sized by `materi_history_variables`.

## Física

Manual (4.18): "The history variable di Prisco plasticity models are
added to the node_dof records. For group_materi_plasti_diprisco you
need to set number_of_history_variables to 11. For
group_materi_plasti_diprisco_density you need to set
number_of_history_variables to 12." The 11 variables: chi tensor
(9 values), beta (yield surface form factor), rc (preconsolidation mean
pressure); the 12th (density) is only used by the NOT-implemented
diprisco_density model.

## Gotchas

- The camclay group checks still require the exact string
  `materi_history_variables` (the alias is diprisco-specific, the
  manual gives no cross-model alias).
