# direct_plastic_strain — developer notes (K13, 2026-10-06)

## Where

- `stress.cc` `set_stress()`: inside the plastic-iteration loop, before
  the direct cut-offs, the stress is saved in `direct_sig_before` and
  `direct_f` is reset to `NO_YIELD_F`; after the loop (block "K13",
  before the material stiffness fill) `direct_plastic_strain()` solves
  `C : deps = sig_before - sig_after` and adds `deps` to `inc_epp`
  (subtracted from `inc_epe` and `new_epe`); `new_f` takes `min(direct_f,
  0)` when no incremental law set it (`NO_YIELD_F`) or when larger.
- `stress.cc` `direct_plastic_strain()`: 6x6 tensor-form system (rows
  sig_11 sig_22 sig_33 sig_12 sig_13 sig_23, shear columns `C_ijkl +
  C_ijlk`), Gauss-Jordan with partial pivoting; C = `Cmem + Cuser + Chyper
  + Chypo` (+ `1/compressibility` on ii,kk), the same elastic stiffness
  as the plane-mode phi_flow return. Singular C -> no plastic strain
  (old behaviour).
- `materi_direct_full_mc()`, `materi_direct_cutoff()`,
  `materi_compression_cutoff()`: new `double &f_direct` output, the
  largest direct-law yield value on the FINAL stress of that function.
  `materi_direct_cutoff()` was restructured (no early returns) so the
  tension and MC-plane values are evaluated after both steps.
- `materi.cc` stores `inc_epp` / `new_f` exactly as for the incremental
  laws: the nodal `eppxx` / `f` are the recovered element values. Follow-up
  (same day): the `f` block writes 0 instead of `NO_YIELD_F` (-1e10) for
  an element without any yield law (Professional `tension_direct1` +
  `materi_plasti_f`: f 0 at every node; GNU before -1e10 / -3.3e9 at the
  shared nodes); test `test-2026/plasti_f_no_law1` RED/GREEN.

## Verification (Professional 02-08-2026, single quad4 probes, E 1, nu 0.25)

| probe | GNU new | Professional |
|---|---|---|
| tension full, eps_xx 2e-3, sigy 1e-3 | eppxx 1.4e-3, eppyy/zz -3.5e-4, epexx 6e-4, f 0 | same |
| tension full elastic (eps 2e-4) | f -7.6e-4 | same |
| tension on plane n=(1,0) / (0,1) | epp as above / f -8e-4 | same |
| compression full, cut / elastic | eppxx -1.4e-3 / f -7.6e-4 | same |
| MC full, cut (sigxy 2.4e-3 -> 0.767e-3) | eppzz 3.29191382e-5, eppxy 2.04114892e-3, f 0 | eppzz same, eppxy 4.08229785e-3 |
| MC full / plane elastic | f -3.0158e-5 / -7.29647e-4 | same |
| tension + MC (c 1 or 1e-3) elastic | f -6.70557e-4 / -6.65124e-4 | same |
| tension + compression elastic | f -6.70557e-4 | same |
| von Mises (incremental, reference) | unchanged | same |

`visc_pl4` (Pro test): eppxx 0.5 (Pro 0.499976, its visco stress is
1.00005), f 0 (Pro 0). Tests `test-2026/direct_plastic_strain_tension1`,
`_mc1`, `direct_yield_value1`: RED on the old binary (rc 1), GREEN new;
the Professional passes them except the `_mc1` eppxy target (shear
convention).

## Hardcoded / pending

- Shear components: tensor convention (Professional: twice) — deliberate
  divergence, see the user page.
- C33 (BACKLOG): `set_stress()` reads `DTIME` from `VERSION_NORMAL` (the
  previous step: the first step of a `_visco` direct law does not cut,
  later steps use the old dt) and the principal tension cap
  (`direct_cap_principal_values`) ignores `_visco`; not changed here
  (stress path, own A/B).
- With formulation TOTAL (`-total` memory without `materi_strain_total`)
  the stored plastic strain now enters `new_epe` of the next step, so a
  `_visco` direct cut persists instead of being re-derived from the
  total strain every step.
