# group_beam_force_moment_plasti — developer

## Files and functions

- `database.cc` — `GROUP_BEAM_FORCE_MOMENT_PLASTI` (DOUBLE_PRECISION, 12,
  class BEAM); enum after `GROUP_BEAM_SHEAR` in `tochnog.h` /
  `tochnog-mod.h`.
- `beam.cc`:
  - `beam_3d()`: a group with the record goes to `beam_local()` also in
    2D.
  - `beam_local()`: with the record, `beam_plastic_hinges()` replaces the
    elastic update f_l += K du_l and turns K into the tangent; the global
    matrix Kg = T^t K T is built after it.
  - `beam_plastic_hinges()`: end moments per bending plane from the
    converged `element_beam_force_moment` (VERSION_NORMAL) plus the
    elastic increment of the relative end rotations
    theta = r -/+ (u2 - u1)/L; return mapping per plane; torsion clipped;
    rebuilds the 12 end forces (shear = (m1 + m2)/L, fx = 0) and
    K = B^t C B.
  - `beam_plastic_plane()`: closest point projection in the energy norm
    of k = EI/L [4 2; 2 4] on the polygon |m1| <= Y1, |m2| <= Y2,
    |m1 + m2| <= L min(Vy1, Vy2) (no active limit, one, or a corner of
    two); exact backward Euler for associative perfect plasticity.
  - `beam_plastic_active()`, `beam_plastic_admissible()`: limit tests with
    relative tolerance `BEAM_PLASTIC_TOL` (1e-10).

## Implementation details

- No extra history: for perfect plasticity the end forces are the state.
- Tangent: a limit is condensed (k - k a a^t k / a^t k a, or 0 at a
  corner) only when the converged moments of the previous step already
  reached it. The Tochnog predictor (previous velocity) overshoots in the
  first iteration of the step where a hinge forms; condensing those
  transient hinges gave a false mechanism (Lapack singular in
  `group_beam_force_moment_plasti1`). Condensing the limits active in the
  previous iterate (VERSION_NEW) was tried and makes Newton cycle with
  period 4; rejected.
- A node whose element ends are all plastic has a zero rotation row; the
  solver regularization (`so.cc` |diag| < 1e-8 -> 1e-8) keeps its
  increment 0 and the node equilibrium m - m = 0 holds
  (`group_beam_force_moment_plasti4`).
- Professional 02-08-2026 (probes, cantilevers E 1 Izz 4): pure clipping
  of the components of `element_beam_force_moment` with elastic
  stiffness; tip displacement 0.1, limit mz1 0.6: mz1 -0.6 with fy 1.2
  (not in equilibrium); force control past the limit gives m1 -3.2
  instead of -2 and diverges (rotations -9.6). Owner decision 2026-11-22:
  a real plastic hinge model (DIVERGENCES BEAM-PLASTIC-HINGE).
