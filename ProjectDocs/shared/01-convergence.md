# The convergence project: GNU → Tochnog Professional

## Goal

Tochnog Professional (the proprietary line, Dennis Roddeman) documents
~1090 keywords in its manual. The GNU open-source line (the last sfnet
release, January 2014) covers ~304 of them. Our project converges the GNU
line towards the Professional feature set: implementing the missing
keywords, with analytical validation and full documentation.

## Methodology

Every feature is implemented in a "work unit" with:

1. **The code** (keyword registration, dispatch, physics) in a `feat(...)`
   commit, with enums kept in sync between `tochnog.h` and
   `tochnog-mod.h`, and a clean build.
2. **Analytical verification**: a dedicated test where the target values
   come from hand calculations (statics, closed-form solutions), not from
   the code's own output. When the Professional manual provides an example,
   we validate against it; when it does not, we build a dedicated test.
3. **Documentation**: a user manual and a developer manual per feature,
   in English, committed with the feature.
4. **A verification log** entry (`SEGUIMIENTO-CONVERGENCIA.md`), recording
   the commit, the date, and how the feature was verified.

The test suite (`scripts/build_safe.sh`) currently runs **201 tests** with
file-level checks, all green on a clean build.

## Scope of work (August 2026)

Implemented, among many others:

- **Material models**: Mohr-Coulomb (classical + hardening/softening),
  Drucker-Prager, cap models (cap1, cap2), Hardening-Soil (elastic +
  plastic + initial history), hypoplasticity Masin family, Sanisand,
  Young's power law, K0, shear factor, stress-pressure history, per-model
  plastic strain initialisers.
- **Control families**: `control_data_*` (activate, arithmetic, copy),
  `control_distribute_*` (statistical distributions with correlation),
  `control_reset_*` (18 keywords), `control_mesh_*` (generate, extrude,
  mirror, rotate, interface conversion), `control_repeat_*` (save +
  calculate), `control_materi_*` gates.
- **Boundary conditions and loads**: `bounda_*` families (unknown, alternate,
  normal, water, radial/cylindrical, time factor), `force_element_edge/volume`
  with restriction variants, `groundflow_flux_edge_normal*` (11 keywords),
  `condif_heat_*` (34 keywords complete), `groundflow_*` families
  (van Genuchten, phreatic level multiple, seepage, total pressure limit,
  consolidation switches).
- **Interface elements (Carril A)**: a complete interface infrastructure —
  `group_interface_*` (13 keywords), 2D/3D interface elements with
  Mohr-Coulomb and tension laws, thermal/groundflow coupling, mesh
  generation with quadratic-face subdivision, and the
  `control_print_interface_stress*` post-processing family.
- **Post-processing**: the full `control_print_*` family — 37 new keywords
  including `control_print_dof_line/point` (field interpolation along lines
  and at points), `control_print_node*` (with angular/sort/geometry/zero
  variants), VTK extensions, print frequency control, and
  `post_calcul -materi_stress -force` (section forces and moments — see the
  solver finding).

The `control_print` checklist alone went from 21/85 to 63/85 with the
remaining items being dependencies (documented) or deliberate discards
(e.g. the proprietary GiD family, replaced by Gmsh/FRD/VTK/CSV exports).

## Validation against Tochnog Professional

We obtained the Tochnog Professional binary (version 02-08-2026, from the
author's public drive) and built a repeatable comparison harness
(`scripts/compare_professional.sh`). The same models run on both binaries;
section forces, stresses and displacements are compared. Results:

- On well-conditioned configurations, the GNU line now matches the
  Professional (e.g. the tunnel section forces **digit-for-digit**).
- The pathological 1-element-in-thickness bending configuration — which
  exposed a deep solver defect in the GNU line — gives exact statics in the
  Professional (the defect does not survive there).

This harness is the acceptance test for every solver-related change.

## Current state

- 209-test suite, all green (201 + the 8 qsri3d SRI-hex8 tests).
- `control_print` family: 63/85 (implementable subset complete).
- The deep solver defect of the open-source line: diagnosed, fixed, and
  validated against the Professional (see `02-solver-finding.md`).
- The SRI (selective reduced integration) extension: DONE for the quad4
  (2D) and the hex8 (3D, 2026-08-29 — the loaded cantilever recovers
  0.897× of the Euler-Bernoulli deflection vs 0.221× locked). Measured
  caveat: the shear-only SRI hex8 retains zero-energy twist/warping
  modes (the classic limitation of the shear-only selective
  integration of the 8-node brick; the 2D quad4 SRI is stable).
- Open fronts: the coarse-mesh multi-step stress equilibrium (the
  recovered σ of the staggered scheme is not in equilibrium — a solver
  issue, not post-processing); the remaining Professional families
  (safety, support, etc.).
