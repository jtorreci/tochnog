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

- `validation_1` (eulerian NS, p+h refinement): since F4 T6f K21 WU2 /
  C72 (2026-10-12/13) the GNU is steady (-0.057, GAMM range) and the
  Professional drifts (class B C72-V1, see below). The GNU prints "Bi-CG
  iteration broke down" about 600 times on the t = 0..3 run (BACKLOG
  C74). Convected stress dofs are recovered absolutely since `c1a1ae9`
  (`general()` skips their convective terms; see
  [nodal_stress_recovery](nodal_stress_recovery.md)).
- `validation_8`: blocked before the physics by the `-veln` bounda_dof
  pseudo-dof (normal-velocity MPC generation, manual 6.22) - not
  implemented in the GNU nor in the GNU-2014 lineage.
- `node_convection_apply` (6.880) stays registered-without-behaviour.

## Differences with the Professional version

- **explicit `convection_apply -no` next to `materi_displacement`** (GNU
  extension / input dialect; AUDIT 1.1). GNU: Accepted (convection ignored
  with `materi_displacement`). Professional 02-08-2026: Refuses the explicit
  record (`-no` is its default). Status: Input dialect. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`validation_1` steady state** (class B, DIVERGENCES C72-V1,
  Professional-suspect, owner 2026-10-13). GNU velx -0.0570386 /
  -0.0570934 / -0.0571954 at t = 1 / 2 / 3; Professional -0.0391 /
  -0.0179 / +0.0034, dt-converged (dt 0.02 / 0.01 / 0.005: 0.00313 /
  0.00344 / 0.00358 at t = 3), identical with `global_element_dof_apply
  -no` (1e-14), unchanged with `convection_stabilization -no`,
  `timestep_predict_velocity -yes` and density 0. Both codes accumulate
  the compressibility pressure at the integration points linearly in
  time (max |p| 1.17 / 1.27 / 1.75 GNU, 1.17 / 1.28 / 1.77 Pro at t =
  0.3 / 1 / 3), so that is not the cause. Guard:
  `test-2026/validation1_steady1` (+ `stokes_group2` at t = 1).
