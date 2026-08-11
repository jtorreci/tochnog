# control_print_tabular

## Implementación

- **CSV + SQLite**: `print_tabular()` in `print_tb.cc`, invoked from the
  control loop in `top.cc` (same pattern as `control_print_vtk`):
  `db_active_index(CONTROL_PRINT_TABULAR, icontrol)` && `ival[0]==-YES`.
- **Keyword registration**: `database.cc`, data_class CONTROL, type INTEGER,
  data_length 1, `ival[0]` holds the `-yes`/`-no` flag.
- **Derived magnitudes**: `calc_derived()` in `derived.cc` — a
  `template<int D> Tensor` helper, `matrix_jacobi` for the 3x3
  eigenvalues, von Mises from the invariants, Tresca as `sig1-sig3`.
- **SQLite access**: `sqlite.cc` / `sqlite.h` — a C++ RAII `SqliteDB`
  wrapper, guarded by `SQLITE_USE` (set in `tn_sqlite.h`). Enabled at
  build time via `makefile` (`SQLITE_OBJ`, `-lsqlite3`).
- **Dof detection**: same approach as `print_vtk.cc` — `dof_label` +
  `dof_scal_vec_mat` (`-SCALAR`/`-VECTOR`/`-MATRIX`).

## Diseño / decisiones

- CSV is always written (no external dependency). SQLite only when the
  build has `SQLITE_USE=1`; otherwise a warning is printed.
- The SQLite schema is **long format** `primary_data(node,dof,t,value)`
  to support arbitrary dofs without a fixed-column schema.
- CSV is written in **append mode**; the header is emitted only on the
  first call (detected by whether the file already exists). This makes the
  CSV grow across time steps like a time series.

## Detalles

- `db_version_copy(VERSION_NORMAL, VERSION_PRINT)` +
  `renumbering(...)` before reading `NODE_DOF` (same as print_vtk).
- `stress_indx(kdim,ldim)` maps Voigt components; the 6-component `sig[]`
  passed to `calc_derived` is `[xx,yy,zz,xy,xz,yz]`.
- `matrix_jacobi` does **not** sort its eigenvalues; `calc_derived`
  sorts them descending (sig1 >= sig2 >= sig3) before reporting.
- Hardcoded: `data_length=1` for the keyword; derived columns hardcoded to
  the first stress dof found (`sig_indx`).

## Pendiente

- SQLite `meta` table is currently empty (no metadata rows written).
- Derived magnitudes only for the first `-MATRIX` dof; multi-material
  problems with several stress dofs would need a loop.
- No option yet to select specific dofs to export (exports all).
