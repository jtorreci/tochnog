# group_materi_memory -total_linear

## Semantics (manual Professional 2.2.1)

`-total_linear` is the total Lagrange formulation with LINEAR engineering
strains and NO rigid-body rotation:

- the total deformation gradient `F` refers to time 0;
- the linear engineering strain is `0.5(F + F^T) - I`;
- the DIFFERENCE of these strains between two time points is the
  incremental strain;
- the stress at the new time point is the previous stress PLUS the extra
  stress from the incremental strain (an INCREMENTAL update, not a
  recomputation from the total strain).

This is what makes a stress set by `control_reset_dof` survive into the
next solve: the reset only changes the stored stress dof, not the
displacement field, so the incremental strain is zero and the new stress
equals the reset value.

## Implementation

Two loci, both required (item 3 is a later, separate fix of the
compressibility increment):

1. **stress.cc, `set_stress()`** — the formulation selection used to force
   `formulation = TOTAL` for `-TOTAL`/`-TOTAL_PIOLA`/`-TOTAL_LINEAR`. The
   `-TOTAL_LINEAR` case was removed from the forcing condition, so the
   `INCREMENTAL` default (declared at the top of the function) is kept.
   The `check_unknown` calls for the three memory types are preserved in
   a separate block. With `formulation==INCREMENTAL` the elastic update is
   `new_sig = C*inc_epe + rotated_old_sig` (stress.cc:1406-1411); for
   `-TOTAL_LINEAR` `rotated_old_sig == old_sig` (no rotation,
   materi.cc:423-427), which is exactly the manual's "previous stress plus
   extra stress".

2. **materi.cc, `set_deften_etc()`** — the `-TOTAL_LINEAR` strain branch
   computed `new_ept` (total) and `inc_ept = new_ept - old_ept`. `old_ept`
   is loaded from the solved dof only when `materi_strain_total` is
   present (materi.cc:255-256); without that keyword `old_ept` is
   zero-initialized, so `inc_ept` collapsed to the TOTAL strain. Feeding
   that into the (now incremental) stress update double-counts the
   previous stress. Fix: when `!materi_strain_total`, derive `old_ept`
   from the OLD deformation gradient
   (`0.5*(old_deften + old_deften^T) - I`). `old_deften` is always built
   from the previous converged displacements (materi.cc:1462-1473), so
   `inc_ept` is a true increment independent of the strain dof. The
   subsequent `array_add(old_ept, inc_ept, new_ept)` (materi.cc:1600) then
   reconstructs the same total as before.

3. **stress.cc, `set_stress()`, compressibility (2026-09-29)** — the
   linear-memory branch (`-TOTAL_LINEAR`/`-UPDATED_LINEAR`/`-UPDATED_AREA`)
   now takes the incremental volumetric stress from `inc_epe`
   (`tr(inc_epe)/compressibility`) instead of `tr(new_epe - old_epe)`.
   `old_epe` is 0 unless `materi_strain_elasti` is initialized, so the old
   form added the TOTAL elastic strain every step. **Deliberate divergence
   from the Professional**: the Professional 02-08-2026 shows the same
   defect for `-total_linear`. Probe = Pro `interac1` (2 bars in series,
   bar 2 compressibility 1, analytic -0.005) with `-total_linear` and
   `materi_displacement`, loading to t = 0.02:

   | dt | Pro sigxx | GNU before | GNU after |
   |---|---|---|---|
   | 1e-2 | -0.0075 | -0.0075 | -0.005 |
   | 1e-3 | -0.0099002 | -0.0099002 | -0.005 |
   | 1e-4 | -0.009999 | -0.009999 | -0.005 |

   The Professional (and the pre-fix GNU, to all digits) is time-step
   dependent and tends to -0.01 (bar 2 rigid, its strain ~0), which is the
   signature of the accumulated total strain. The Professional's own
   `-updated_linear` gives the analytic -0.005 on the same model. No corpus
   test combines `-total_linear` with compressibility. See BACKLOG C9.

## Why not just switch to INCREMENTAL

Switching the formulation alone (without the `old_ept` fix) leaves
`inc_ept` as the TOTAL strain when the strain dof is absent, so the
incremental update adds the full total stress on top of the stored stress
and the result is too stiff. Measured on `total3` (1D bar, `-total_linear`,
target `disx = 0.5`): formulation-only switch gives `disx = 0.1`
(RUNFAIL); with both loci it is `0.5` (PASS).

## Verification (2026-09-13)

Base `abf9c29`, branch `dev/total-linear-formulation`, clean build.

- `total3` (only sfnet-corpus test with `-total_linear`): PASS,
  `disx = 0.5` (no regression).
- Reset-only probe (`delete3` with the delete records removed, reset
  `-sigyy 1`): GNU baseline `sigyy = 0`, fixed GNU `sigyy = 1`,
  Professional 02-08-2026 `sigyy = 1`. The reset now survives.
- sfnet 197-test corpus: 160 PASS / 37 RUNFAIL / 0 PARSE, byte-identical
  per-test classification to the baseline (only `total3` uses
  `-total_linear`). Internal suite 16/16.
- `delete3` (`-total_linear` + reset + delete factor 0.4, target
  `sigyy = 0.6`): with the test's own 5 timestep iterations
  `sigyy = 0.59862`; with 30 iterations `0.59999996`. Professional
  converges to `0.6` in its 5 iterations.
- `mpc5` (`-total_linear` + reset + delete factor 0.5, target
  `sigyy = 0.5`): GNU default `control_timestep_iterations` is 1, so
  `sigyy = 1.0`; with 20 iterations `0.50008`. The Professional's default
  is 3 (`number_iterations 3` in its `.dbs`) and reaches `0.5`.

The residual at the tests' own iteration counts is CONVERGENCE SPEED, not
the formulation: the fixed stress update converges to the Professional
value. Orthogonal mechanisms documented elsewhere contribute to the slow
fixed point (`mpc_linear_quadratic` value-constrained ties; the
`element_delete_factor` momentum scaling which is diagonal-only vs the
Professional's full-matrix scaling). See
`ProjectDocs/SEGUIMIENTO-CONVERGENCIA.md`.

## Pending

- `control_mesh_delete_geometry_factor` momentum scaling: the GNU stores
  the element rhside plus only the matrix DIAGONAL (`element_matrix_delete`,
  elem.cc:1044-1067) and on the next assembly overrides the diagonal and
  the rhside with `factor * stored` and skips the element constitutive
  computation (elem.cc:352-374). The Professional stores only
  `element_rhside_delete` (no `element_matrix_delete`) and scales the full
  matrix, which converges faster. The exact Professional scaling is not
  uniquely determined from the `.dbs` alone.
