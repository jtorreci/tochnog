# inertia_apply

## Files and functions

- `general.cc` — in `general()` (lines 63-68). Read with the global record
  and the per-timestep override:
  ```c
  db( OPTIONS_INERTIA, 0, &options_inertia, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  db( CONTROL_OPTIONS_INERTIA, icontrol, &options_inertia, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  ```
  The lumped transient (inertia) term of each principal dof is added only
  when `options_inertia==-YES` (line 293):
  ```c
  if ( options_inertia==-YES || !principal_unknown ) { ... }
  ```
  with the per-dof inertia weight: material density for `-MATERI_VELOCITY`,
  groundflow capacity C for `-GROUNDFLOW_PRESSURE`, `dens·capacity` for
  `-CONDIF_TEMPERATURE`. Non-principal dofs (history variables, strains,
  displacement integration, ...) always keep their transient term.
- `truss.cc` — in `truss()` (lines 31-32, 86-88): the same record gates the
  truss mass term.
- `database.cc` — keyword registration: `options_inertia` (INTEGER,
  `data_length = 1`, `no_index = 1`) and `control_options_inertia` (CONTROL
  class). The Professional names `inertia_apply` / `control_inertia_apply`
  are aliases resolved in `db_number()`.
- `check.cc` — requires unknowns to be initialized.

## Implementation details

- `GET_IF_EXISTS`: absent records keep the local initializer `-NO` (since the
  u-p consolidation sprint, commit `0560dcf`; the GNU legacy initializer was
  `-YES`, aligned with the Professional default `-no`).
- The GNU legacy default `-yes` made every analysis with a material density
  pseudo-dynamic (mass matrix active without `inertia_apply`). Measured on
  ground14 of the corpus (safety piping/lifting, artesian column): with the
  legacy default the momentum oscillated around the drained equilibrium
  (sigma'(0) = -9.53 at t = 1 s vs the drained -10) and the corpus targets
  failed; with the static default the drained state is reached in the first
  step (-10.0000001423), reproducing the Professional binary to 8 digits.
- The Professional record is multi-valued (one switch per principal dof);
  the GNU registration stores a single value that applies to every principal
  dof (the Pro manual allows the single-value form explicitly).

## External dependencies

- Core `db()` accessor; the per-dof inertia weights come from
  `get_materi_density` / `groundflow_data` / `group_condif_capacity`.

## Hardcoded parameters / pending refactorings

- Multi-value `inertia_apply` (one switch per principal dof, Professional
  form) is not supported: the GNU reads only the first value and applies it
  to every principal dof. Extending `OPTIONS_INERTIA` to a variable-length
  record would need the per-dof switch bookkeeping of `dof_type`.
