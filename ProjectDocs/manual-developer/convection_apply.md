# convection_apply / control_convection_apply / convection_stabilization (developer)

## Where

- `database.cc` — `db_number()` (the alias translation chain around
  line 7947, next to the `inertia_apply` / `control_inertia_apply`
  aliases):
  - `"convection_apply"` -> `OPTIONS_CONVECTION` (manual 6.395)
  - `"control_convection_apply"` -> `CONTROL_OPTIONS_CONVECTION`
    (manual 6.113, indexed per-timestep override)
  - `"convection_stabilization"` -> `OPTIONS_STABILIZATION`
    (manual 6.396)
- `general.cc` — `general()` already consumes the targets: it reads
  the switch through `convection_apply_switch()` (default `-NO`, K03),
  computes the convective velocity from the materi
  velocity dofs / `condif_flow` / groundflow velocity and assembles the
  convective + artificial-diffusion terms; `OPTIONS_STABILIZATION`
  switches the peclet-based upwind (`-NO` disables it, `-MAXIMAL` also
  clamps condif temperatures).
- `top.cc` — reads/validates `OPTIONS_STABILIZATION` (`-DYNAMIC` is not
  implemented and errors).

## Implementation details

The GNU data items and their consumers already existed (options_* names);
the corpus Professional files use the `convection_apply` /
`convection_stabilization` names, so only the three name translations
were missing. Both data items are INTEGER no-index switches whose values
(`-yes`/`-no`/`-maximal`) resolve through the generic value reader
(input.cc: `db_number(&str[1])`).

Defaults: `OPTIONS_CONVECTION` defaults to `-NO` (K03, 2026-10-06; it
was `-YES`). See the K03 section below. `OPTIONS_STABILIZATION` defaults to `-STATIC`
(top.cc), which in general.cc behaves as the minimal-peclet stabilization
(`-yes` of the Professional); only `-NO` and `-MAXIMAL` are
special-cased.

## K03: default -no and scope (2026-10-06)

- `general.cc` `convection_apply_switch()`: single source of truth for
  every reader (`general()`, `beam.cc`, `truss.cc`, `spring.cc`,
  `conspr.cc`; prototype in `tochnog.h` / `tochnog-mod.h`). Returns
  `OPTIONS_CONVECTION` overruled by `CONTROL_OPTIONS_CONVECTION` at the
  current `ICONTROL`; default `-NO` (manual 6.395). The local defaults of
  the readers were `-YES` (hardcoded in five places).
- `general()`: the `lagrange` flag (`options_convection==-NO ||
  materi_displacement || (materi_velocity && -follow_material)`) now
  gates only the material-velocity part of the convective velocity. The
  `condif_flow` and groundflow-velocity part of the temperature
  convection uses `lagrange_flow` (same condition WITHOUT the switch).
  Measured on the Professional 02-08-2026: sfnet condif1
  (`group_condif_flow`) and condif5 (groundflow velocity) give temp 1
  with `convection_apply` `-no`, `-yes` and absent; condif2 (material
  velocity) gives 0.5 absent/`-no` and 1 with `-yes`.
- Side effect: `stress_recovery_absolute` (nodal stress recovery) uses
  the same switch, so inputs without the record on a `-fixed_in_space`
  mesh with velocity unknowns switch to the absolute recovery and no
  longer convect the stress/strain dofs (validation-suite: ~57 tests
  change outputs, mostly round-off, a few by 1-13 %; where the
  Professional runs them - cmat_gate, mesh_act_grav, mdirect_comp - the
  new values are closer to it).
- The beam/truss/spring/contact-spring checks ("set OPTIONS_CONVECTION
  to -NO") now only fire with an explicit `-yes` on a `-fixed_in_space`
  mesh.
- Tests: `test-2026/convection_apply_default1.dat` (RED with the old
  binary: temp 1, target 0.5 fails; GREEN 0.5 = Professional) and
  `convection_apply_yes1.dat` (1 = Professional), run by `build_safe.sh`.

## Verification

Verified against the Professional binary (user-supplied 25-10-2023,
.dbs): condif2 gives node-2 temperature 1.000000 (GNU 0.992028, within
the 2e-2 tolerance); tube1 gives post-point vely 1000.49 (GNU 1066.59,
within 100). The .dbs of the Professional stores the records under the
same input names (`convection_apply -yes`, `convection_stabilization
-yes`).

## Pending

- `validation_1` (eulerian NS, quad9, p+h refinement): the GNU transient
  of the velocity self-advection diverges from the Professional from
  step 2 (final vx at the post point +0.0723 vs -0.0391; the no-convection
  variant of both codes matches). The convective assembly of the mixed
  velocity/stress eulerian formulation in `general()` needs a dedicated
  A/B sprint.
- `validation_8`: blocked before the physics by the `-veln` bounda_dof
  pseudo-dof (normal-velocity MPC generation, manual 6.22) - not
  implemented in the GNU nor in the GNU-2014 lineage.
- `node_convection_apply` (6.880) stays registered-without-behaviour.
