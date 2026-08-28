# materi_stress staggered u-σ scheme — fixed point consistency (lot C/D)

**Status**: implemented 2026-08-28 (branch `documentation-improvement`).
Companion docs: `ProjectDocs/DIAG-SOLVE-MIXTO.md` (§12: the exact
mechanism and the fix), `ProjectDocs/VALIDACION-PROFESIONAL.md` (§9:
the post-fix arness), `ProjectDocs/SEGUIMIENTO-CONVERGENCIA.md`.

## What this document covers

The `materi_stress` "mixed" formulation of the GNU is NOT a monolithic
saddle-point system: the stress dofs never enter the global matrix
(`dof_principal` is not assigned to them, `input.cc:439-444`; `so.cc:298`
excludes them), so the linear solve is a velocity-only SPD system and the
stress field is advanced by a lumped diagonal update after each solve.
The scheme is STAGGERED:

```
v_{k+1} = v_k + (dt*K_uu)^-1 * (P - B^T*sigma_constit(v_k))
sigma_{k+1} = sigma_constit(v_k) = sigma_old + C_full*eps(v_k*dt)
```

This lot makes the FIXED POINT of that scheme equal to the solution of
the ELEMENT formulation (the displacement-stiffness solution `K_elem*u =
P`), instead of the full-constitutive (locked) solution `K_full*u = P`
that the old scheme converged to regardless of the element.

## The mechanism (measured, code lines)

1. **Momentum matrix** (elem.cc → materi.cc): `dt*K_uu` with
   `K_uu = B^T*D_elem*B` — `matrix_atba(new_b, ddsdde_total, stiffness,
   work, MSTRAIN, nnol*ndim)` (materi.cc:519) and
   `tmp = volume*dtime*stiffness[...]` (materi.cc:651). For the SRI
   quad4 the tangent is split `D = D_norm + D_shear` and the shear is
   integrated with 1 Gauss point at the centroid (materi.cc:654-662).
2. **Momentum RHS** (materi.cc:508 + 617): `P - V*B^T*sigma` where
   sigma = `new_sig`, the CONSTITUTIVE stress of the current velocity
   iterate: `new_sig = sigma_old + C:inc_epe` (stress.cc:1182-1183)
   with `inc_ept = B*v*dt` (set_deften_etc, materi.cc:1131).
3. **Stress update** (dof.cc:155, `parallel_new_dof_diagonal`):
   `sigma_node += RHS_sigma/LHS_sigma` with
   `RHS_sigma = V*h*(new_sig - sigma_old)/dt` (materi.cc) minus the
   inertia `V*h*(sigma_iterate - sigma_old)/dt` (general.cc:253-255,
   always active for non-principal dofs) and `LHS_sigma = V*h/dt`
   (general.cc:256). The inertia term cancels the accumulated iterate,
   so **sigma_new = new_sig** (a REPLACEMENT, not an accumulation).
4. **Velocity solve** (so.cc:740): `node_dof_new += solve_b` — the
   velocity ACCUMULATES across the equilibrium iterations.

### The old fixed point

At the fixed point `dv = 0`, so the momentum matrix drops out:

```
P = B^T*sigma_constit(v*) = B^T*sigma_old + dt*B^T*C_full*B*v*
v* = (dt*K_full)^-1 * (P - B^T*sigma_old)        (K_full = B^T*C_full*B)
sigma* = sigma_old + C_full*eps(u*)
```

The fixed point is governed by the FULL constitutive — the SRI matrix
only shapes the transient (the old "0.3125x" SRI value at 2 iterations),
and at equilibrium the full-rule shear re-enters through `-B^T*sigma`,
cancelling the SRI (measured 0.2316x ~ the plain 0.2315x).

### The stress recovery defect

The sigma dofs are recovered by the h-weighted average
`sigma_node = sum_gp V*h*sigma_gp / sum_gp V*h`. With the
node-containing quadratures (the default 2x2 Lobatto corners of the
quad4, the quad9/hex8 Lobatto rules) h is the Kronecker delta and the
recovery is exact. With the 2x2 GAUSS rule (interior points at
+-1/sqrt(3), switched by the SRI quad4 — polynom.cc:421) the h-weighted
average dilutes the corner values (measured 0.577x for the bilinear):
the section forces (post_calcul -materi_stress -force) read low nodal
stresses even when the displacement field is correct.

## The fix (D)

### D-c — element-consistent momentum feedback (materi.cc, `sri_quad4`)

The momentum right-hand side must carry the ELEMENT internal force
`B^T*sigma_old + dt*K_elem*v` instead of the full-constitutive stress.
For the SRI quad4:

- the current-iterate shear increment of the feedback stress is zeroed:
  `sigvec[stress_indx(0,1)] -= 2*sri_g*inc_ept[stress_indx(0,1)]`
  (the old shear prestress `sigma_old_xy` is kept);
- the reduced 1-point shear internal force
  `-dt*K_shear*v` (the same matrix-vector product the SRI shear block
  adds to the momentum matrix, scaled by 1/npoint) is added to the
  momentum RHS in the velocity block.

With the element-consistent feedback the velocity map collapses:
`v(k) = v* = (dt*K_elem)^-1*(P - B^T*sigma_old)` for every k >= 1 —
ONE-pass convergence, no drift, and the fixed point is the element
solution:

```
AFTER:  v* = (dt*K_elem)^-1 * (P - B^T*sigma_old)   (the element)
        sigma* = sigma_old + C_full*eps(u*)         (the physical)
```

For the plain quad4/quad9/hex8 (`K_elem = K_full`) the feedback is
unchanged — the fixed point and the transients are byte-identical to
the old scheme (verified on gforce7q4_ref/gffq4/gforce7_ref/qsri_beam2d).

### D-b — consistent sigma recovery (materi.cc + general.cc)

The sigma dof update uses the bilinear Lagrange EXTRAPOLATION of the
Gauss-point values to the nodes (the "same B at the node") for the
NORMAL stress components, which are superconvergent at the Gauss points
(measured: the SRI gp sigma_xx = 96% of the analytic beam stress; the
extrapolated nodal sigma_xx = 93.9%). The SHEAR components keep the
h-weighting: the Q4 shear is NOT superconvergent (the interpolation
error dominates; the centroid-biased average is the better estimate).

Implementation: `sri_stress_recovery_weight()` (miscel.cc) computes the
2D tensor product of the 1D Lagrange extrapolation weights on the
2-point Gauss grid (row-major-from-bottom convention, the same as pol()
and the SRI centroid B). The weight is applied in the sigma-RHS
(materi.cc, normal components only) and in the inertia term of
general.cc (MATERI_STRESS dofs, normal components only — the components
are identified by `kcomp = iuknwn - stres_indx` in {0,3,5}). For every
node-containing quadrature the weight reduces to h, so ONLY the SRI
quad4 (2x2 Gauss) changes. The momentum feedback does not read these
dofs (it uses the fresh constitutive stress), so the recovery change is
purely an OUTPUT improvement.

Note: `materi()` and `general()` gained the `ipoint` parameter (the
integration-point index of the current call) and `materi()` also the
`new_dof[]` array (the current velocity iterate, needed for the
`K_shear*v` product). Signatures updated in tochnog.h/tochnog-mod.h and
the calls in elem.cc.

## Verification

- SRI quad4 cantilever (qsri_beam2d_sri): clamp moment 0.0750 =
  0.9375*P*L (the classic Hughes SRI reference; the msf moment about the
  element centroid = P*(L-0.5) exactly), stable at 32 equilibrium
  iterations (no drift). The normal stress field matches the analytic
  bending stress to 94%.
- Plain quad4 / quad9 / hex8: byte-identical (the locked element
  solution for the plain — the lock is ELEMENT physics, opt-in SRI).
- Arness vs Professional (scripts/compare_professional.sh):
  gforce10/gforce13 converge with the axial N EXACT (1.0000x);
  gforce7q4 converges; gforce7_ref keeps N/M >= 0.99x.
- Suite: 199/199 + file checks green.

## Known limitations (NOT fixed by this lot)

- The plain quad4/hex8 shear LOCK (element physics; the opt-in SRI
  keyword fixes the quad4; hex8 SRI is future work). The scheme now
  converges to the element's OWN solution, locked or not.
- The section SHEAR pollution: the raw sigma_xy of the bilinear Q4 is
  dominated by the interpolation error (the FE shear strain of the
  beam field); the nodal recovery cannot fix it. The Professional's own
  nodal sigma_xy is equally polluted (force7q4 clamp sigma_xy = -90.9
  vs the beam tau(y=0) = 0); its exact section statics do not come from
  the raw nodal stress.
- The monolithic mixed solve (option C) is not needed: the staggered
  scheme with the element-consistent feedback converges to the element
  solution in one pass.

## Files touched

- `materi.cc` — element-consistent feedback (D-c); sigma recovery weight
  (D-b); `ipoint`/`new_dof` parameters.
- `general.cc` — inertia weight for the MATERI_STRESS dofs (D-b);
  `npoint`/`ipoint` parameters.
- `miscel.cc` — `sri_stress_recovery_weight()`.
- `tochnog.h` / `tochnog-mod.h` — signatures.
- `elem.cc` — call sites with the new parameters.
- `scripts/build_safe.sh` — qsri_beam2d_sri check updated to the fixed
  0.9375x moment.
- `scripts/compare_professional.sh` — magnitude (s) indices + correct
  per-model nodes; the arness input `.dat` BC fixes (gforce7q4: velx on
  the left edge {1,4}; gforce10/13: bottom face {1,2,3,4}) are local
  (the test .dat files are gitignored).
