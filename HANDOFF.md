# HANDOFF — DEV-D: force16 / force17 / mpc9 residuals

Sprint: 4-dev parallel convergence (GNU → Professional, binary 25-10-2023).
Branch: dev/fmpc (from 19e995a), worktree /tmp/opencode/tn-dev-fmpc.
Started: 2026-02-06 (sprint slot 45 min corpus window; commands <10 min).

## State

- [x] Worktree verified clean on dev/fmpc @ 19e995a.
- [x] Read DIAG-SOLVE-MIXTO §16-17 + SEGUIMIENTO dev/segv + dev/veln entries.
- [x] Corpus inputs: test/other/{force16,force17,mpc9}.dat.
- [x] Pro reference values measured (see Findings; evidence /tmp/tn-dev-fmpc-pro/).
- [x] GNU baselines (all 3 rc=1; /tmp/tn-dev-fmpc-gnu/base/).
- [x] Root causes decoded EXACTLY (closed-form matches, see Findings).
- [ ] Implement mpc.cc eligibility extension (fixes force16 + mpc9).
- [ ] Implement calcul_force.cc _force/_inertia consumption (fixes force17).
- [ ] Blast radius: mpc*/force*/post* families before/after; mpc1/2/7/8 byte-identical.
- [ ] build_safe suite 16/16 green.
- [ ] Final commit + Engram observation.

## Findings (measured, Pro 25-10-2023 vs GNU 19e995a)

### mpc9 (target post_element_force_result = 5/11 = 0.4545454545455)
Pro state (full-dof probe): v=[0, 5/11, 5/11, 1], sigma=5/11 all nodes,
reactions n1=+5/11, n4=-5/11. bounda_time 20 1. = Dirichlet v4=1 (NOT a
force). GNU: v2=v3=0 value-constrained (slave row dropped), sigma=[0, 5/6].
Pro semantics = STATIC CONDENSATION of the explicit user tie v3=v2:
master row2+row3 fold: (k1+k2)*v2 - k2*v4 = b2+b3, with the Dirichlet
column moved to the rhs -> (11/6)v2 = (5/6)*1 -> v2 = 5/11, condensed
diag 11/6 = 1.833333 (prior attempt CONFIRMED). Same mechanism as DIAG
§16 but for EXPLICIT user records with free masters.

### force16 (target 5.034285714286 = 881/175)
GNU field [1, 3.93667, 3.93667, 5.89444] vs Pro [1, 6.0342857, 6.0342857,
7.2928571]. Decoded in closed form: momentum matrix = M/dt + K (dt=1),
lumped m=0.5/node (rho=1, L=1), gravity +4.905/node, v1=1 Dirichlet,
mpc_element_group-generated tie v3=v2. GNU = value-constrained (drops node-3
row): 1.5 v2 = 4.905+1 -> v2 = 3.93667 EXACT match of the GNU. Pro =
CONDENSED tie: 3 v2 - v4 = 2G+1, 1.5 v4 - v2 = G -> v2 = 1056/175 EXACT
match. So the Pro CONDENSES mpc_element_group-generated ties too. post then
= sigma*A = v2-v1 = 5.0342857 automatically (calcul_force.cc needs NO
change for force16).

### force17 (target 9.931723274805)
Velocities ALREADY match Pro exactly [0.07447243528581, 0.06208101510514, 0].
Pro nodal sigma = average of adjacent element stresses (GNU prints its own
averages of DIFFERENT element stresses - downstream only). Pro post result =
sum over group-1 element-1 nodes of (gravity + inertia lumped shares) =
10 - (0.5*0.0744724 + 0.5*0.0620810) = 9.931723274805 (13 digits exact).
Lumped shares = volume*h_i*rho with INITIAL geometry: bar2 -> rho*L/2 per
node. _force gates gravity, _inertia gates inertia. Elastic internal forces
self-cancel per element (no elastic term in the result). GNU 1D branch:
_force/_inertia parse-only; current GNU value = 2*N1 = -2.478284 (section
device applied at both ends).

### Blast-radius map (corpus)
- mpc1: explicit user tie, master node1 DIRICHLET -> NOT eligible (masters
  must be free) -> byte-identical.
- mpc8: explicit user tie, SLAVE node1 directly bounda'd -> NOT eligible
  (slave must be free) -> byte-identical.
- mpc7: mpc_element_group-generated ties (nodes 5,6,7 -> free masters of
  element 1) -> WILL be condensed. Conforming patch: value-constrained and
  condensed solutions coincide (slave residual ~0); MUST verify 1.0 +-1e-8
  and byte-identical output. mpc7 has no bounda on the tied slaves; watch
  the post_node_rhside_ratio/iteration path too.
- mpc6 (mpc_element_group temp ties, RUNFAIL): may change values (not in
  the neutrality list) - report.
- mpc2 (mpc_geometry record): untouched. mpc3/4/5 (lq range): already
  eliminated, unchanged. validation_8 (-veln from_bounda): own branch,
  unchanged.
- post7/force5/6: no mpc records; post7 has no _force/_inertia flags ->
  untouched. Only force17.dat uses post_element_force_force/_inertia in the
  whole corpus (grepped).

## Next step

Implement: (1) mpc.cc - extend tie-map eligibility to (a) the
mpc_element_group generated range (index-range only, mirroring the lq class)
and (b) explicit user records with slave-free AND all-masters-free; eligible
user-tie slaves join the per-pass bounded-flag clearing list (list renamed
mpc_veln_applied -> mpc_applied). (2) calcul_force.cc 1D branch - consume
post_element_force_force/_inertia: accumulate the lumped gravity/inertia
shares (share_i = sum_ip wrule*2|jac|*N_ip(i), initial coords; gravity via
force_gravity_calculate + get_materi_density; inertia via cached previous-
step velocity dofs, first call = zero initial velocities; dtime from db
DTIME VERSION_NEW), flags gate the two contributions, elastic section
device skipped when either flag active. Then build + A/B + blast radius.

## Evidence paths

- Pro binary: /tmp/opencode/tn_prof/prof25/linux_64_bit_user_supplied/tochnog_user_supplied
- Corpus: /tmp/opencode/tn_prof/tochnog_version_02-08-2026/test/
- Pro probes (print-enabled): /tmp/tn-dev-fmpc-pro/{mpc9v,force16p,force17p}_pro.out
- GNU baselines: /tmp/tn-dev-fmpc-gnu/base/*_gnu.out
- Prior killed attempt (leads only, re-derived): /tmp/fmpc_base, /tmp/fmpc_d, /tmp/tnfmpc
- Build log: /tmp/tn-dev-fmpc-build.log (binary OK; suite needs
  validation-suite/test-2014 which the worktree lacks - run the suite
  against another worktree's suite dir with MY binary).

## Evidence paths

- Pro binary: /tmp/opencode/tn_prof/prof25/linux_64_bit_user_supplied/tochnog_user_supplied
- Corpus: /tmp/opencode/tn_prof/tochnog_version_02-08-2026/test/
- Scratch: /tmp/tn-dev-fmpc-*
