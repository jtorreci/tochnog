# groundflow_pressure_gradient

## Where

- Initia parsing: `input.cc` (`groundflow_pressure_gradient` branch):
  declares the vector dof `pres_grad` (ndim components) with the dof
  type `GROUNDFLOW_PRESSURE_GRADIENT`.
- Enum: `GROUNDFLOW_PRESSURE_GRADIENT` in `tochnog.h` +
  `tochnog-mod.h` (in sync).
- Names: `database.cc` name entry + `db_initialize()` basename branch
  (`pres_gradx`, `pres_grady`, `pres_gradz` via `pres_grad_indx`).

## Pending

The dof is registered but not yet FILLED: the groundflow element
(`groundfl.cc`) computes the pressure gradient at the integration
points (`grad_new_unknowns[...pres_indx]`); writing the recovered
gradient into `pres_grad_indx` per node (element-average recovery) is
not done. Nothing in the current corpus targets the value of the
gradient dof, so the registration is parse-only for now.
