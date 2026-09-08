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
- [x] IMPLEMENTED mpc.cc: explicit-user-tie condensation (slave-free + all-masters-free) + per-pass clearing list + eg-range EXCLUSION (see the mpc7 note).
- [x] IMPLEMENTED calcul_force.cc: _force/_inertia lumped-share consumption (1D branch).
- [x] mpc9 rc=0 (5/11 at all nodes = Pro). force17 rc=0 (9.93172 = Pro to 13 digits; gravity-only probe = 10, inertia-only = −0.0682767).
- [x] force16: stays RUNFAIL **byte-identical to base** — DIAGNOSED: the Pro condenses the mpc_element_group tie (closed form 1056/175; prototype reproduced the Pro field exactly and closed mpc6 too) but the GNU staggered map has a pre-existing ~1e-8/iteration inconsistency (mpc3 shows it; accepted at ±1e-3) and mpc7's ±1e-8 target with the same 2-iteration cap misses by 1.9e-9. The eg range is explicitly EXCLUDED (mpc7/mpc6 byte-identical); the future closure = fix the map inconsistency (so.cc/materi.cc) + activate in_eg_range (kept in mpc.cc).
- [x] Blast radius: mpc1-8, force1-15/18, post1-7 byte-identical (base binary from 19e995a rebuilt in /tmp/tn-dev-fmpc-base; "CPU time" lines excluded as wall-clock noise); only mpc9/force17 rc 1→0.
- [x] build_safe suite 16/16 (run manually against the main repo's validation-suite/test-2014 inputs — the worktree lacks the suite dir; command in the log below).
- [x] Full corpus with the new binary: see /tmp/tn-dev-fmpc-corpus_result.txt (final summary line).
- [x] Docs: DIAG-SOLVE-MIXTO §18 + manual-developer/force_mpc_residuals.md (new).
- [ ] Engram observation (before returning).

## Commits (dev/fmpc)

1. e11f153 docs(dev-fmpc): initial HANDOFF
2. 032c42e docs(dev-fmpc): measured Pro semantics + plan
3. 3e91abd mpc(fmpc): condense explicit user ties with free slave/masters + clearing list
4. e50e938 force(fmpc): consume post_element_force_force/_inertia (1D)
5. docs(diag): §18 + 6. docs(dev): force_mpc_residuals.md
6. (pending) final HANDOFF + any mpc.cc comment fixups

## Next step (for the orchestrator)

- Merge dev/fmpc. mpc9/force17 close; force16 documented (needs the
  §16-map inconsistency fix in so.cc/materi.cc — DEV-B territory — then
  the in_eg_range activation in mpc.cc; prototype measured: force16 rc=0,
  mpc6 rc=0, mpc7 field = the Pro's).

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

## Evidence paths

- Pro binary: /tmp/opencode/tn_prof/prof25/linux_64_bit_user_supplied/tochnog_user_supplied
- Corpus: /tmp/opencode/tn_prof/tochnog_version_02-08-2026/test/
- Pro probes (print-enabled): /tmp/tn-dev-fmpc-pro/{mpc9v,force16p,force17p}_pro.out
- GNU baselines: /tmp/tn-dev-fmpc-gnu/base/*_gnu.out
- Prior killed attempt (leads only, re-derived): /tmp/fmpc_base, /tmp/fmpc_d, /tmp/tnfmpc
- Build log: /tmp/tn-dev-fmpc-build.log (binary OK; suite needs
  validation-suite/test-2014 which the worktree lacks - run the suite
  against another worktree's suite dir with MY binary).
- Suite 16/16 run: cd <main-repo>/validation-suite/test-2014 && timeout 120
  <worktree>/build/tochnog t.dat for the 16 build_safe tests (log
  /tmp/tn-dev-fmpc-suite_*.out).
- Blast A/B: /tmp/tn-dev-fmpc-gnu/blast/{base,new} (base binary rebuilt
  from 19e995a in /tmp/tn-dev-fmpc-base).
- Corpus result: /tmp/tn-dev-fmpc-corpus_result.txt (harness
  /tmp/tn-dev-fmpc-run_corpus.sh).
