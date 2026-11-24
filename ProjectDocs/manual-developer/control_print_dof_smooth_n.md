# control_print_dof_smooth_n

## Implementación

- Stored as a CONTROL INTEGER record (one value). Read with
  `GET_IF_EXISTS` inside `print_dof_smooth_apply()` (`print_hi.cc`):
  the number of smoothing passes of
  [`control_print_dof_smooth_dof`](control_print_dof_smooth_dof.md).
- Registered in `database.cc` with
  `data_required = CONTROL_PRINT_DOF_SMOOTH_DOF` (so it can only be
  used together with the smoothing record, same index).

## Diseño / decisiones

- Default = 10 passes when the record is not given (manual 6.271: "if
  you don't specify this optional number of smoothings it is done 10
  times").
- `number_of_smoothings < 1` -> `db_error(CONTROL_PRINT_DOF_SMOOTH_N)`.

## Detalles

- Each pass averages the values of the previous pass over the node
  neighbourhood (see
  [`control_print_dof_smooth_dof`](control_print_dof_smooth_dof.md) for
  the formula). Repeated passes converge to the mean value of the field
  (verified: the 1D chain `0..4` converges to `2` after the default 10
  passes, `|v-2| <= 0.0625`).

## Pendiente

- None.
