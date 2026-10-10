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
- Reconstruction gotcha (validation-suite `groundflow_interface`, 2026-09-15):
  a transient groundflow model that declares `group_groundflow_capacity` /
  `group_interface_groundflow_capacity` but relies on the legacy GNU default
  silently loses its storage term after `0560dcf`, and the free nodes do not
  fill from the interface flux. Adding `inertia_apply -yes` restores the
  Professional `.dbs` pressure field (node5 pres `1.99999894661` vs the
  archived `1.9999990002`). The reconstructed test had been labelled
  "faithful" only because the `interface_element()` abort kept it from
  running.
- Reconstruction gotcha (the 15 solver-divergence validation-suite tests,
  2026-09-16): `iface_3d`, `iface_3d_slip`, `iface_3d_stress`, `iface_gen`,
  `iface_gen_geom`, `iface_gen_geom_off`, `iface_gen_method`,
  `iface_gen_method_gen`, `iface_gen_quad9`, `iface_mc`, `iface_mc_mem`,
  `iface_tangref`, `contact_block`, `contact_ctrl_apply`, `contact_heatgen`
  have a material density but no `inertia_apply` record and no complete
  Dirichlet set. After `0560dcf` their assembled matrix is genuinely singular
  (rigid-body modes -> `d^T A d ~ 0` -> Bi-CG breakdown -> LU retry diverges),
  which had been misattributed to the honest Bi-CG criteria. Adding
  `inertia_apply -yes` restores the `rho/dtime` diagonal and makes 13/15
  converge to their stored targets. Two residuals remain and are **not**
  solved by the record: `iface_tangref`
  (`element_interface_force_tang2` became a per-IP array in `bf74f28`; the
  target is tied to the old scalar) and `iface_gen_method_gen` (the
  reconstruction drop-list dropped the input `element_geometry` records).
  Details: `CORPUS-STALE-TARGETS.md` §"Inertia-apply reconstruction".

## External dependencies

- Core `db()` accessor; the per-dof inertia weights come from
  `get_materi_density` / `groundflow_data` / `group_condif_capacity`.

## Hardcoded parameters / pending refactorings

- Multi-value `inertia_apply` (one switch per principal dof, Professional
  form) is not supported: the GNU reads only the first value and applies it
  to every principal dof. Extending `OPTIONS_INERTIA` to a variable-length
  record would need the per-dof switch bookkeeping of `dof_type`.

- The `tet10` lumped inertia/capacity of the principal dofs uses the HRZ
  fractions 1/36 (corner) and 4/27 (mid-edge) instead of `h`
  (`general.cc` `general()`, hardcoded; C31, F4 T4). A generic
  consistent-mass-diagonal HRZ for any element would replace the
  hardcoded fractions.

## Differences with the Professional version

- **Lumped mass / capacity of the `tet10`** (deliberate divergence,
  project owner 2026-10-09; DIVERGENCES TET10-HRZ). GNU: positive HRZ
  lumping, corner V/36 and mid-edge 4V/27 (the diagonal of the
  consistent mass scaled to the element total). Professional 02-08-2026:
  row-sum lumping of its 5-point rule, corner -V/20 and mid-edge V/5;
  the negative corner masses make dynamics unstable (sudden gravity on a
  48-tet10 column: the Professional solver fails, a GNU row-sum build
  diverges to 5e8). Static results are unaffected. Code: `general.cc`
  `general()`. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
