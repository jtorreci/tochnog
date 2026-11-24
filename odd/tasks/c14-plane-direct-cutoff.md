# C14 — plane-mode direct cut-off (mohr_coul_direct / tension_direct with a normal)

## Objective
Make the plane-mode direct cut-off follow the Professional manual (6.725–6.740) and
the Professional binary 02-08-2026, except the phi_flow factor, where the GNU uses
the consistent elastic return (M/G) as a documented divergence.

## Problem (evidence: scratchpad c14 probes, BACKLOG C14 item 4)
- `materi.cc:508-522` merges the MC and tension normals into one `direct_normal`.
- `stress.cc` `materi_direct_cutoff`: no default sigy = 0 when `mohr_coul_direct`
  is present, no principal-stress tension cap when the tension law has no own
  normal, phi_flow ignored.
- `database.cc:5078`: `tension_direct_normal` requires `tension_direct`
  (the Professional accepts it alone).

## Rule to implement
1. Tension first. sigy = `tension_direct` value, else 0 when `mohr_coul_direct`
   is present. Own tension normal (`tension_direct_normal` / `_automatic`) →
   cap sig_n <= sigy on that plane; otherwise cap the principal stresses at sigy.
2. Then MC on the MC plane with the updated stress:
   k = (M/G)·tan(phi_flow);
   tau_new = (c − sn·tanφ + k·tanφ·tau) / (1 + k·tanφ), floor max(·,0);
   sn_new = sn − k·(tau − tau_new). Only sn and tau change.
   Professional uses k = 2·tan(phi_flow) (constant; equals M/G only at nu = 0)
   → documented divergence (BACKLOG).

## Decisions
- 2026-10-02 user: phi_flow factor M/G (consistent elastic return), not the
  Professional constant 2; document the divergence.

## Tasks
- [x] T1 Separate MC and tension normals; tension step (default sigy, plane vs principal cap); drop data_required (tension_direct_normal only; Pro rejects _automatic alone, probe c23). Route: delegated (multi-file). Commit `d17461e`.
- [x] T2 phi_flow return with k = (M/G)·tanψ, M/G from the element elastic stiffness (stress.cc direct_flow_factor; C = Cmem+Cuser+Chyper+Chypo+1/compressibility built in set_stress), fallback 2. Route: delegated (same writer). Commit `77f9282`.
- [x] T3 Closed-form regression tests validation-suite/test-2026/direct_cutoff_plane_*.dat (9, on-disk only, gitignored); corpus A/B. No commit (untracked corpus).
- [x] T4 Docs (user + developer), SEGUIMIENTO, BACKLOG (C14 items 1 and 4 resolved; C15 phi_flow divergence; C16 large-strain elastic; C17 stale materi_direct_mc target). Commit `3bde0bb`.

## Checks
Clean build; new tests RED on old binary / GREEN on new; validation-suite 2014
(233/1), test-2026, sfnet (187/10) unchanged except explained changes.

## Progress
- Branch `dev/c14-plane-direct-cutoff` created from `9264b8d`.
- Extra Pro probes (c20-c25): TN alone without any law = no effect; MC without normal + tension normal = full principal MC (no principal cap); separate planes confirmed (c24); floor beyond apex sn_new = sn - k tau (c22, c25).
- Tests: old binary 8/9 RED (tension_normal1 passes on both), T1 binary phi_flow1/2/floor1 RED, T2 binary 9/9 GREEN; Pro 8/9 (phi_flow2 = divergence).
- A/B (scripts/audit_all_variables.py run, 45 s): sfnet identical; test-2026 only new tests + round-off in element_normal3d_column1; validation-suite: cd_geom, cd_method changed (still pass), materi_direct_mc target fails (233 -> 232 PASS), explained and confirmed against Pro with the same constraints (collapse at t = 0.7 in both) -> BACKLOG C17 (target decision pending, outside edit surface).
- build/tochnog = T2 clean build.

- 2026-10-03 user: fix materi_direct_mc input (vely on nodes 1 2 3 4); GNU rc 0, target sigxy 1 met; BACKLOG C17 resolved (validation-suite 233 PASS).

## Next step
Merge into dev/masin-kernel-surface (user decision); push is a user decision.
