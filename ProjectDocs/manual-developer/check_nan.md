# check_nan

## Archivos y funciones

- `check.cc` → `check_nan_results( long int check_nan_switch )`
  (`check.cc:1099`) — iterates all active `NODE_DOF` entries and tests
  `node_dof[iuknwn] != node_dof[iuknwn]`; on a match it prints
  `Error: NAN detected in node_dof of node <inod> dof <iuknwn>.` plus
  `The solution may have diverged.` to `std::cout` and calls
  `exit(TN_EXIT_STATUS)`.
- `top.cc:805-806` — called at the end of `step_close()` (defined at
  `top.cc:636`): `if ( check_nan==-YES ) check_nan_results( check_nan );`.
- `initia.cc` (`initia.cc:94`) — global `long int check_nan=-NO;`.
- `top.cc:128` — reads the keyword once at startup with
  `db( CHECK_NAN, 0, &check_nan, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS )`.
- `database.cc:290-293` — keyword registration:
  `strcpy(name[CHECK_NAN], "check_nan")`, `type = INTEGER`,
  `data_length = 1`, `no_index[CHECK_NAN] = 1`.
- Enum `CHECK_NAN` in `tochnog.h:167` / `tochnog-mod.h:160` (must stay in sync).

## Detalles de implementación

- Early return: if the switch is `-NO`, or if there are no nodes
  (`db_max_index( NODE, max_node, VERSION_NORMAL, GET )` returns `max_node < 0`),
  the function does nothing.
- NaN detection uses the self-inequality idiom `node_dof[iuknwn] !=
  node_dof[iuknwn]`, which is `true` only for NaN in IEEE-754. It does NOT rely
  on `<cmath>` / `std::isnan`.
- Iterates dofs from `0` to `nuknwn-1` (global degrees of freedom count) for
  every active `NODE_DOF` entry. A non-active node is skipped via
  `db_active_index( NODE_DOF, inod, VERSION_NORMAL )`.
- The check runs inside `step_close()` after the previous-RHS handling, so it
  inspects the converged step results.

## Dependencias externas

Uses the database layer (`db_dbl()`, `db_active_index()`, `db_max_index()`) and
the global `nuknwn`. No external libraries and no `<cmath>`.

## Parámetros hardcodeados / refactorizaciones pendientes

- The `!=` self-comparison is deliberate (NaN-safe) but cryptic; a named helper
  such as `is_nan(double)` with a comment would document the intent.
- Error output goes to `std::cout`; routing it through `pri()` / `tn.log`
  would be consistent with `exit_tn()` error reporting.
- Only `node_dof` is scanned. Element results (`element_dof`) and auxiliary
  arrays are not covered; extending the scan would catch NaNs earlier.
- `max_node` is recomputed on every step via `db_max_index`; caching it across
  steps would avoid the repeated scan when the mesh is static.
