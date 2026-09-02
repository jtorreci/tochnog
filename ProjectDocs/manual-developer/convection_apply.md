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
  `OPTIONS_CONVECTION` / `CONTROL_OPTIONS_CONVECTION` per element
  (default `-YES`), computes the convective velocity from the materi
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

Defaults: `OPTIONS_CONVECTION` is `-YES` in the GNU element integration
(an eulerian `-fixed_in_space` analysis convects unless the switch says
`-no`; with a follow-material mesh the convection is suppressed anyway by
the `lagrange` flag). `OPTIONS_STABILIZATION` defaults to `-STATIC`
(top.cc), which in general.cc behaves as the minimal-peclet stabilization
(`-yes` of the Professional); only `-NO` and `-MAXIMAL` are
special-cased.

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
