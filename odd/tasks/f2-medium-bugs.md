# F2 — medium GNU-suspect bugs (contact, truss/beam/support, power-law elasticity)

## Objective
Close the medium-size GNU-suspect clusters of the all-variables audit, each
verified against the Professional 02-08-2026 or an analytic value.

## Tasks
- [ ] T1 K06: power-law elasticity semantics (`myoung6_e2`: GNU E=3000 → −4.0385 vs Pro E=1000+3000 → −5.3846; `myoung6_e3`; `mpower` poisson_power). Pro adds `group_materi_elasti_young` to the power law.
- [ ] T2 K17: truss / beam / support_edge_normal divergence on Pro inputs (`trubea2/3` disz 2e5 vs 0.00333; `support_edge_normal_3` disy 124.75 vs 5.01e-4; `truss8`; `dynamic3`; `truss2/3/14`, `update2` ERRDUMP).
- [ ] T3 K05 + C5: contact (`contac1` no interaction, GNU node 2 velx 0.03 vs Pro 0.02/0.01; `validation_10` sigyy ≈0 vs −21031; `contact_target_geometry` alias aborts "Data is allocated in a parallel loop" in 6 runs and drops an index-0 contact in examp9).

## Rules
Principle: manual/analytic decides; Pro is the oracle when the manual is
silent; divergences and many-target changes need a user decision. Each fix:
RED/GREEN test in test-2026 registered in build_safe.sh,
`TN_CHECK_ONLY=1 ./scripts/build_safe.sh` green (226 runs, 121+ OK, 0 FALLO),
corpus A/B with the same runner (validation-suite 233/10, test-2026, sfnet),
docs + SEGUIMIENTO + AUDIT/BACKLOG. Time-box ~1.5 h per root cause.

## Progress
- Branch `dev/f2-medium-bugs` from `e693018`.
