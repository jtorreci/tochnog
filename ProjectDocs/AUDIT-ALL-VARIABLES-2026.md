# All-variable audit: Tochnog GNU vs Tochnog Professional (2026-09-29)

GNU build: `dev/masin-kernel-surface` @ `4cfaa0e`, clean `build_safe.sh --clean`
build, frozen copy sha256 `f93906f2…`. Professional oracle: static binary
02-08-2026 (sha256 `97923ac4…`). Nothing in the source tree was changed; the
only repository additions are this report and the re-run scripts
`scripts/audit_all_variables.sh` + `scripts/audit_all_variables.py`.

## Summary

The GNU → Professional convergence cycle measured one target variable per test.
This audit ran each test on both binaries and compared every variable of the
final state. For 482 tests, both binaries finished and the final states were
compared:

- **261 match** every compared variable, within tolerance.
- **27 differ by less than 1 %** only.
- **194 differ by 1 % or more** in at least one variable. In 112 of these, none
  of the mismatching variables is the one the test's target checks. This is the
  blind spot that hid the damage-stress defect fixed in `4cfaa0e`.
- The **damage models are clean**: 0 of the 7 damage tests mismatch after `4cfaa0e`.

The strongest evidence comes from the Professional's own tests, which carry
targets written by the Professional's authors:

- The Professional meets every one of its own targets that the audit could
  evaluate: 0 failures in 210 tests.
- The GNU misses these targets in **30 tests**. Most of them are analytic
  checks, such as edge heat flux, force_point, K0, tension cut-off and
  -updated_linear.

The clusters below hold the differences. The main GNU-suspect findings are:

1. `-updated_linear` gives wrong answers in the GNU. It breaks the analytic
   2-material bar (−0.010 vs −0.005). The same model with `-updated` or
   `-updated_without_rotation` is right. The defect predates `c945833`.
   **Fixed 2026-09-29** (stress rotation + compressibility increment); see
   section 5.1.
2. With a static phreatic level and no pressure dof, the GNU builds the initial
   gravity stress with the wet density everywhere: −1995 against the analytic
   −1250. **Fixed 2026-09-29** (wet/dry density + static pore pressure on the
   skeleton); see section 5.2.
3. A set of feature bugs found by the Professional's analytic tests:
   - edge heat flux doubled
   - `force_point` sign flipped
   - K0 ignored
   - `tension_direct_normal_automatic` ignored
   - `plasti_bounda_factor` ignored
   - `stress_pressure_history` switch timing
   - `expansion_volume` has no effect
4. Contact: 1D contact has no interaction, and the `contact_target_geometry`
   alias crashes ("Data is allocated in a parallel loop") or silently drops the
   contact. This is the known backlog item C5, but broader than documented.
5. The power-law elastic moduli differ in meaning: the Professional binary adds
   `group_materi_elasti_young` to `young_power`. This refutes the "Professional
   semantics" claim in the GNU implementation notes.
6. Output-state defects:
   - the default of `convection_apply` differs from the Pro manual
   - nodal strain is recovered with an operator different from the Pro's (the
     stress is recovered with the Pro's)
   - strain and stress lag one iteration when only one iteration is run
   - direct plasticity does not store plastic strain

---

## 1. Method

### 1.1 Corpora and dialect

| corpus | source | how it was run |
|---|---|---|
| sfnet 2014 (197) | `external-downloads/sfnet/extracted/test` | GNU dialect: the GNU ran the original file. The GNU and the Pro both ran the mechanically converted file. |
| validation-suite test-2014 (233) | `validation-suite/test-2014` (reconstructed from `.dbs`; `print_mesh_dof.dat` excluded) | same as sfnet |
| validation-suite test-2026 (15) | `validation-suite/test-2026` | same as sfnet |
| Pro `test/other` (330), `test/validation` (15), tutorials, `easy_tochnog`, `incremental_driver` (16) | Professional distribution | as distributed, on both binaries |

**The comparison always uses the same input text on both binaries.** For the
GNU-dialect corpora this is the converted file (the `gnuconv` side). Two extra
comparisons show the effect of the conversion:

- `convcheck`: the GNU on the original file vs the GNU on the converted file.
- `origcmp`: the GNU on the original file vs the Pro on the converted file.

**Conversion rules** (`convert()` in the script). Each pure rename maps a GNU
legacy spelling onto the Professional spelling. The GNU accepts the Professional
spelling as an alias (`database.cc db_number`):

- `bounda_unknown`→`bounda_dof`
- `options_mesh`→`mesh`
- `options_inertia`→`inertia_apply`
- `options_convection`→`convection_apply`
- `options_stabilization`→`convection_stabilization`
- `options_processors`→`processors`
- `control_options_solver`→`control_solver`
- `group_truss_young`→`group_truss_elasti_young`
- `-trussbeam`→`-truss_beam`
- `groundflow_phreaticlevel`→`groundflow_phreatic_level`
- `contact_geometry`→`contact_target_geometry`
- `druckprag`→`druck_prag`
- `force_element_{edge,volume}*`→`force_{edge,volume}*`
- `group_groundflow_porosity`→`group_porosity`
- `…plasti_boundary`→`…plasti_bounda`
- `…plasti_cap`→`…plasti_cap2`
- `-sizetot`/`-sizedev`→`-size_tot`/`-size_dev`
- `contactspring*`→`contact_spring*`

Some changes remove records or add a post-processing request, without changing
the model:

- `target_*` records are dropped. They only drive the GNU exit code, and the Pro
  rejects GNU dof names in them.
- `convection_apply -no` is dropped next to `materi_displacement`. The Pro
  refuses the explicit record, `-no` is its default, and the GNU ignores
  convection when `materi_displacement` is on.
- `post_calcul -groundflow_pressure -total_pressure` is added when a pressure
  dof exists and no `post_calcul` is present. The Pro requires it.

Four approximations are tagged `!` and checked by `convcheck`:

- `-updated_without_rotation`→`-updated_linear`. The Pro has no
  `-updated_without_rotation`; the Pro's own copies of these tests use
  `-updated_linear` or `-total_linear`.
- The `number_of_integration_points` line is dropped. The Pro does not know it
  and uses its per-operation default.
- `materi_velocity_integrated` is dropped when `materi_displacement` is present,
  but only when every `node_dof` input record is zero.
- `mohrcoul`→`mohr_coul` (same `phi c phiflow` layout).

Conversion neutrality measured on the GNU (418 tests where `convcheck` ran):

- **Pure renames: 156 neutral, 3 changed.** All 3 are `contact_geometry` →
  `contact_target_geometry`: the GNU alias is broken (see K05).
- **Approximations: 227 neutral, 32 changed.** Almost all 32 come from
  `-updated_without_rotation`→`-updated_linear` (see K01).

### 1.2 Comparison

- **Final state only.** The records compared are:
  - `node_dof` (all base dofs, plus the derivative slots, reported separately)
  - `post_point_dof`
  - `node_dof_calcul` (post_calcul)
  - `node_rhside` (principal dofs only: reactions and residuals)
  - `.his` histories and `incremental_driver_result.txt` time histories, when
    both binaries wrote them
- **Dofs are mapped by name, never by column.** Both `.dbs` dumps carry a
  `dof_label` record in `node_dof` order, and each value is read through it.
  Labels present in only one binary are reported, for example the GNU `sph` and
  `epphs*` or the Pro `hshis`.
- **Nodes are paired by number when the reference coordinates agree.** The
  reference coordinates are `node_start_refined`, or `node` when that record is
  absent. Otherwise nodes are paired by coordinates; this happens with refined
  or generated meshes whose numbering differs. `node_dof_calcul` uses the same
  pairing.
- **Tolerance.** A value is flagged when
  `|g−p| > 1e-4·max(|g|,|p|) + 1e-6·F + 1e-10`. `F` is the magnitude of the
  variable's family in that test: the max over all nodes of the whole tensor or
  vector. `sig*` share one family, `ept/epe/epp` share one, and `vel/veli/dis/acc`
  share one. This keeps round-off in identically-zero components from flagging.
  Sign flips (both values above 1e-3·F) and NaN/inf are flagged too.
  "Magnitude" in the tables is `frel = max|g−p| / F`, and a variable counts as
  a mismatch at `frel ≥ 1 %`.
- **Run classification.** `OK` means a complete `.dbs` was written. `PARSE`
  means the input was rejected. `ERRDUMP` means a `.dbs` was written after a
  fatal error; such runs are never compared, and the non-fatal "Bi-CG iteration
  broke down" message does not count. Otherwise the run is `TIMEOUT` (45 s per
  test and binary) or `RUNFAIL`. The Pro returns exit code 0 even on errors, so
  the classification reads the log text.
- **Targets.** Every `target_item`/`target_value` of the original input is
  evaluated with the same code on both dumps. The Pro does not report target
  outcomes itself.
- **Hand validation of the mapping:**
  - `mk0` (K0 oedometer): `sigxx=-0.75`, `sigyy=-1.5`, analytic and equal in
    both to 2e-8.
  - `validation_2` (consolidation): the `node_dof` fields match to 1e-11.
  - `spring4` and `interac1` (analytic 1D bars) were followed by hand.
  - `axisym1`: the dof labels and nodes were checked record by record.

### 1.3 Diagnostics used for the triage

- **Five-iteration re-run.** Every mismatching test was run again on both
  binaries with `control_timestep_iterations i 5` for every timestep index. This
  separates one-iteration lag artefacts from physics differences:
  - 6 mismatching tests vanish
  - 10 are partly explained
  - 150 persist
- **Long runs.** 25 long Pro tests were re-run at 180 s. The list is in the
  script: the geotechnical models, the Masin hypo7/8/9, validation_8/9/12/14 and
  a few others.
- **Bisect of the `-updated_linear` defect.** The trees of `c945833^` and
  `c945833` were built in a scratch directory with `git archive`. Both give the
  same wrong answer, so the defect is older than that commit.

## 2. How to re-run

```bash
scripts/build_safe.sh --clean            # only clean builds are reliable
TOCHNOG_PROF_BIN=/tmp/opencode/tn_prof/tochnog_version_02-08-2026/src/tochnog \
  scripts/audit_all_variables.sh /tmp/tochnog-audit
# -> /tmp/tochnog-audit/report_tables.md   (all tables of sections 3-4)
#    /tmp/tochnog-audit/runs/summary.json  (per test: status, every mismatch,
#                                           worst node, targets on both dumps)
#    /tmp/tochnog-audit/runs/clusters.json (cluster membership)
```

- The whole run takes about 45 min with 16 cores.
- Knobs: `AUDIT_TIMEOUT`, `AUDIT_JOBS`, and `AUDIT_LONG`/`AUDIT_LONG_TIMEOUT` for
  the long re-runs.
- Single test: `python3 scripts/audit_all_variables.py compare GNU_DIR PRO_DIR NAME`.
- The cluster rules in `assign()` are the triage snapshot of this report. They
  include explicit lists of the representatives verified by hand, so update
  them when a cluster is fixed.

## 3. Coverage

| corpus | tests | run on both | GNU only | Pro only | neither |
|---|---|---|---|---|---|
| sfnet 2014 (GNU dialect, converted) | 197 | 100 | 77 | 0 | 20 |
| validation-suite test-2014 (reconstructed, converted) | 233 | 116 | 110 | 2 | 5 |
| validation-suite test-2026 (own tests, converted) | 15 | 14 | 1 | 0 | 0 |
| Pro distribution test/other (Pro dialect, as is) | 330 | 237 | 18 | 28 | 47 |
| Pro distribution test/validation | 15 | 10 | 1 | 2 | 2 |
| Pro tutorial / easy_tochnog / incremental_driver | 16 | 5 | 0 | 7 | 4 |
| **all** | **806** | **482** | **207** | **39** | **78** |

**The Pro fails** where the GNU ran in 207 tests: 143 PARSE, 61 RUNFAIL and
3 TIMEOUT. The main causes are:

- GNU-only initialisation records (`residue`, `materi_density`,
  `materi_history_variables` for hypo, the GNU per-model plastic-strain initia)
- `group_interface`
- record-arity differences
- `node_support_edge_normal_force`
- `group_element_selective_reduced_integration`
- legacy `condif_convection`
- `tendon`
- `materi_displacement` checks
- Pro-only requirements

**The GNU fails** where the Pro ran in 39 tests: 19 PARSE, 7 RUNFAIL,
10 ERRDUMP and 3 TIMEOUT. The main causes are:

- the `easy_*` macros (not implemented)
- `-quad9` not available for axisymmetric analysis (4 tests)
- the `contact_target_geometry` crash (6 tests)
- solver divergence on truss models (truss2/3/14, update2)
- a handful of missing keywords

**Neither binary finished** in 78 tests. After the 180 s re-run, the large
geotechnical models still did not finish on either binary:

- `foundation1`
- `dam_building`
- `slope_classical_numerical`
- `vajont_slope_fs`
- `ground11_*`
- `heat_exchanger_pile_*`
- `validation_8`
- `validation_9`
- `tutorial_3`

## 4. Results

### 4.1 Match summary

| corpus | compared | match | minor only (<1 %) | mismatch ≥1 % | …only in untargeted vars | …node_dof fine (post_point/calcul only) |
|---|---|---|---|---|---|---|
| sfnet 2014 | 100 | 50 | 6 | 44 | 22 | 3 |
| validation-suite 2014 | 116 | 50 | 3 | 63 | 34 | 10 |
| test-2026 | 14 | 10 | 2 | 2 | 2 | 0 |
| Pro test/other | 237 | 146 | 12 | 79 | 50 | 7 |
| Pro test/validation | 10 | 1 | 3 | 6 | 4 | 0 |
| Pro tutorial/easy/incremental_driver | 5 | 4 | 1 | 0 | 0 | 0 |
| **all** | **482** | **261** | **27** | **194** | **112** | **20** |

### 4.2 Targets, evaluated identically on both dumps

| targets | GNU | Pro | tests |
|---|---|---|---|
| Pro-authored (Pro distribution) | PASS | PASS | 173 |
| Pro-authored | **FAIL** | PASS | **30** |
| Pro-authored | NA | PASS | 7 |
| GNU-authored (sfnet/validation-suite) | PASS | PASS | 161 |
| GNU-authored | PASS | **FAIL** | **35** |
| GNU-authored | FAIL | PASS | 11 |
| GNU-authored | FAIL | FAIL | 3 |

"NA" means the item type is not supported by the evaluator (for example
`-post_node_result`), or the dof name does not exist in that binary.

- **Pro-authored targets.** The Pro never fails its own targets, which confirms
  both the oracle and the evaluator. The 30 GNU failures are the most direct
  GNU-suspect evidence.
- **GNU FAIL / Pro PASS on GNU-authored targets (11).** These are the
  `-updated_linear` cases: the GNU fails its own old target once the input is
  converted, and the Pro meets it.
- **GNU PASS / Pro FAIL on GNU-authored targets (35).** These are mostly GNU
  extensions or GNU-dialect semantics (clusters K06, K07, K09, K12, K16, K03).

### 4.3 Which variables and features differ

Tests with a `node_dof` variable off by ≥1 % of its family scale:

| family | tests |
|---|---|
| stress | 97 |
| total strain | 70 |
| displacement / velocity | 70 |
| plastic strain | 23 |
| groundflow | 19 |
| kappa / f | 10 |
| temperature | 7 |

Mismatch rate by feature (tests with ≥1 % / tests compared):

| feature | rate |
|---|---|
| druck_prag | 6/7 |
| hypo | 8/11 |
| viscous | 17/24 |
| cap | 3/5 |
| hardsoil | 3/5 |
| tension cut-off | 11/19 |
| groundflow | 20/35 |
| `-updated_linear`/`-updated_without_rotation` | 79/141 |
| axisym | 18/36 |
| refine/remesh | 34/65 |
| contact | 8/17 |
| mohr_coulomb | 6/14 |
| memory `-total*` | 66/206 |
| vonmises | 2/10 |
| condif | 14/102 |
| **damage** | **0/7** |
| camclay | 0/2 |

## 5. Clusters (ranked)

Each mismatching test is assigned to one primary cluster by `assign()`. The
score is tests × median magnitude (capped at 1) × relevance. Relevance is
1 for post-processing, 2 for general mechanics, and 3 for the geotechnical core:
soil plasticity, groundflow, consolidation and the initial stress state. The
verdict is one of:

- **GNU-suspect**
- **PRO-suspect**
- **BOTH-plausible** (a formulation or convention choice)
- **INPUT-artifact** (the conversion, or GNU/Pro dialect semantics)
- **UNKNOWN**

The priority list in section 6 re-orders by judgement.

| # | id | cluster | tests | median magnitude | verdict | representative: evidence |
|---|---|---|---|---|---|---|
| 1 | K01 | GNU `-updated_linear` memory model | 13 | 0.80 | **GNU-suspect → FIXED 2026-09-29** (section 5.1) | `interac1`: 2-material 1D bar, analytic σ = −0.005. **Pro −0.00500** (exact). The GNU gives −0.0100 with `-updated_linear`: element 2 has strain 0 but stress −0.01, so equilibrium is broken. The GNU gives −0.00501 with `-updated` or `-updated_without_rotation`. Same result on the `c945833^` build, so the defect is older. Other cases: `druckpr1` target 0.11547 (GNU `-updated_without_rotation` 0.11547, Pro 0.11547, GNU `-updated_linear` 0.0693); `mdp_shear` DP analytic 34.64 (Pro 34.641, GNU `-updated_linear` 34.145). |
| 2 | K09 | Groundflow `-pres` dof datum | 11 | 1.0 | BOTH-plausible / INPUT-artifact | `ground4`: GNU `pres` 0 vs **Pro 500**. The Pro stores the head h, as its manual says ("Primary dof is the hydraulic pressure head h"). The total pore pressure (−500) and `sigyy` (−1247.5) are identical. `ground14`: safety factors identical (2.3333 / 2.0). Also in this cluster: the GNU-only `groundflow_total_pressure_limit` and `pressure_atmospheric`, which the Pro ignores. |
| 3 | K14 | Hypoplasticity: history slots and stress path | 8 | 1.0 | UNKNOWN (partly GNU-suspect) | `hypo7/8/9` (Masin): stress equal to 1e-4, but the GNU leaves `hyhis1..5` = 0 where the Pro fills them (e.g. `hyhis3` 22519). `hypo3`: Pro target −0.1618±0.02, **GNU −0.1034**, Pro −0.1629. `hypo10`: `disx` GNU 4.3e-9 vs target 9.84e-8. |
| 4 | K15 | Other plasticity (MC, DP+refine, von Mises+refine) | 8 | 1.0 | UNKNOWN (pre-existing MC difference) | `validation_14` (DP, refine, delete): `velx` differs 48 %, `eppxx` 25 %. `validation_3` (von Mises plate): `kap` 11 %. `mohrcou2`: GNU −6.932 (its target) vs Pro 0. `mmc_tension`: `sigzz` 2 %. |
| 5 | K12 | post_calcul item conventions | 22 | 1.0 | BOTH-plausible | `calcul1`: the sfnet target (GNU) is [1, 0, 0.707, 0.577] and the Pro's own `calcul1` target is [0, 1, 1, 0.816]. Each binary passes its own. `examp10`/`validation_5`: the order or sign of the `-prival` principal values (0.002 vs −0.00402). Also here: `force10`/`force13`/`gforce10` (known C3), the `msf_*` section forces, and the `qsri` patch calcul records. |
| 6 | K05 | Contact | 8 | 1.0 | **GNU-suspect** | `contac1` (1D contact): **GNU node 2 keeps `velx` 0.03**, so there is no interaction. The Pro gives 0.02/0.01 (momentum exchange). `validation_10`: GNU `sigyy` ≈ 0 vs Pro −21031 (target −20538 ± 1e4). The `contact_target_geometry` alias aborts with "Program error… contact_geometry. Data is allocated in a parallel loop" in 6 runs (`vs/contact*`, `examp6/9/24`). In `examp9` it silently drops an index-0 contact. This is the known backlog item C5, which only documents index > 0. |
| 7 | K13 | Direct / visco plasticity output | 7 | 1.0 | **GNU-suspect** (state output) | `visc_pl4`: the stress is 1 in both (the target). The **GNU `eppxx` = 0 and `f` = −1e10**, while the Pro stores `eppxx` 0.49998 = ε − σ/E. So the GNU output violates σ = E(ε − εp). **2026-10-06 (F1 T4): fixed, the direct cut-offs store C⁻¹(σ_trial − σ) and the direct-law f (`visc_pl4` eppxx 0.5, f 0); shear components tensor (Pro twice, BACKLOG C34); `_visco` stress path C33.** |
| 8 | K02 | One-iteration lag of non-principal dofs | 6 | 1.0 | **GNU-suspect** | `spring4` (rigid translation `velx`=1): with 1 iteration the **GNU gives `eptxx` = `sigxx` = −0.01**, while the Pro gives 0 (exact). With 3 iterations both give 0. The post-solve sync in `top.cc` refreshes only the stress slots, and only when `stress_recovery_absolute` is on, so it misses the strain dofs and models with a strain unknown. 6 tests vanish at 5 iterations and 10 improve. **2026-10-06 (F1 T5): parked as design (BACKLOG C35): the integration-point state itself lags (GNU ip velx 0.5, eptxx -0.01; Pro ip velx 1, 0), so the fix is a post-solve material pass, not a nodal sync.** |
| 9 | K16 | GNU extensions the Pro ignores | 12 | 1.0 | INPUT-artifact | `strain_settle`, `mfactor`, `mdirect_gate`, `transv1`, `temp4/5` (`-empty` group: the Pro does not move its nodes), `mvoid*`, `cmat_gate`. The GNU meets its own analytic targets; the Pro returns the unmodified state. |
| 10 | K17 | Truss / beam / support divergence on Pro inputs | 6 | 1.0 | **GNU-suspect** | `trubea2/3`: `disz` **GNU 2e5** vs target 0.00333 (Pro exact). `support_edge_normal_3`: `disy` **GNU 124.75** vs 5.01e-4. `truss8`: `disx` 0.2 vs 0.1. `dynamic3`: `disx` 32.7 vs 60.4. `truss2/3/14` and `update2` also diverge (ERRDUMP) in the GNU. **2026-10-06 (F2 T2): fixed `support_edge_normal_3` (damping missing from the matrix, = Pro), `truss14` (-total was linear) and `truss2/3` (no geometric stiffness), all = Pro. Parked: `trubea2/3` (planar GNU beam vs Pro 3D beam, BACKLOG C37); `truss8` node 1 (inactive-element node, both plausible), `update2` (not a truss: -updated robustness), `truss13` (truss_distribute_mpc not implemented) in C39; `dynamic3` = C11.** |
| 11 | K04 | Pro analytic targets the GNU misses | 11 | 0.5 | **GNU-suspect** | See the list below the table. |
| 12 | K19 | Nonlinear kinematics residual | 4 | 1.0 | BOTH-plausible / UNKNOWN | `elasti1` (tria3 stretched 300 %): `eptyy` 3.4375 vs 3.0, a finite-strain measure. `incnav2`: GNU σ −0.142 vs Pro 0. `update1`: initial stress in the `-ra -from -to` form, rotated 90°. |
| 13 | K06 | Power-law elasticity semantics | 4 | 0.56 | **GNU-suspect (semantic)** | `myoung6_e2` (E clamped to 3000): **GNU −4.0385 (E=3000)** vs Pro −5.3846 (E = 1000 + 3000). `myoung6_e3`: GNU E=800 vs Pro 1800. The Pro adds `group_materi_elasti_young` to the power law. The GNU implementation notes attribute "the law IS the modulus" to the Pro; the binary refutes this. `mpower` (`poisson_power`) behaves the same way. **2026-10-06 (F2 T1): `mpower` was the opposite case, the Pro REPLACES the poisson ratio and the GNU added C(young,0.3) + C(young,nu(p)); fixed, GNU -0.4389/-1.4101 vs Pro -0.4417/-1.4123. `young_power`: without the `young` record Pro = GNU; with it the Pro adds C(young) + C(E_law) while its history variable reports E_law (3000 on `myoung6_e2`) and the manual says the law is the modulus: needs a decision (BACKLOG C36).** |
| 14 | K21 | Viscous fluid / Navier-Stokes | 9 | 0.68 | UNKNOWN | `incnav1/4/5`, `tube1`, `examp1`: the stresses differ 6-100 %. The sfnet and Pro variants disagree with each other. |
| 15 | K07 | Positional `node_dof` initial state over different dof layouts | 6 | 1.0 | INPUT-artifact | `mhardsoil_gp0`: the initial `-2` sits at slot 72, which is `sigxx` in the GNU but another dof in the Pro (the GNU has `sph`/`epphs*`, the Pro has `hshis`). |
| 16 | K22 | Phreatic level: wet/dry density in the gravity stress | 2 | 1.0 | **GNU-suspect** | `ground18` (= `ground4` with `groundflow_phreatic_level_static` and no pressure dof). Analytic total `sigyy` at the bottom = −750 − 500 = −1250 (the `ground4` comment). **Pro −1247.5, GNU −1995**, which is the wet density over the whole column. With the pressure dof (`ground4`) the GNU is right. `ground8` (multiple phreatic levels): `sigyy` differs 19 %. **2026-09-29: `ground18` fixed, GNU −1247.5 = Pro (section 5.2); `ground8` unchanged (it has a pressure dof, other mechanism).** |
| 17 | K10 | Groundflow velocity | 5 | 1.0 | PRO-suspect | `ground1`: manual v = k ∂h/∂y = 0.1·(−0.5) = −0.05. **GNU −0.05, Pro 0**, at the nodes and at the integration points. This is post-processing only. |
| 18 | K11 | Nodal recovery of total strain | 19 | 0.12 | **GNU-suspect** | `tslv_part`: the velocities and stresses are equal to 1e-6. At the clamped edge, where v ≡ 0, the Pro gives `eptyy` = 0 (exact) and the GNU −2.35e-6; `eptxy` differs by a factor of 3. The GNU recovers the stress with the Pro operator (`06e175d`) but still recovers the strain with the old one. **2026-10-06 (F1 T5): fixed, the strain tensors use the stress recovery weight; `tslv_part` strains = Pro within 4e-11 (`manual-developer/nodal_strain_recovery.md`).** |
| 19 | K18 | Refinement, higher-order and special elements | 8 | 0.27 | UNKNOWN | `disk1`, `force6`, `interface_patch`, `validation_11` (GNU fails the Pro target: 6.09e7 vs 8.59e7 ± 3e6). |
| 20 | K23 | Round-off in zero fields | 3 | 1.0 | BOTH-plausible (false positive) | `dpoint1`: stresses 6e-8 vs 5e-14 in a rigid motion. |
| 21 | K03 | Default `convection_apply` | 2 | 0.5 | **GNU-suspect** | `condif2`: GNU T = 1, Pro T = 0.5. With an explicit `convection_apply -yes`, **both give 1**. The GNU default is `-YES` (`general.cc:54`); the Pro manual 6.395 says the default is `-no`. **Resolved 2026-10-06 (F1 T2)**: default `-no`, see recommendation 6. |
| 22 | K08 | `bounda_time_increment` tables | 2 | 1.0 | PRO-suspect | `mcap1_elast`: the Pro prescribes −1.4 on dofs whose table is all zero. The GNU follows manual 6.146. |
| 23 | K20 | Axisymmetric hoop (zz) nodal values | 5 | 0.05 | BOTH-plausible | `axisym2` (Akin ring, analytic e_t(ID) = 9.964e-4, S_t = 10910): GNU 1.000e-3 / 10000 (element-constant), Pro 1.0498e-3 / 10498 (linear). Moved to F4 (design + corpus A/B), user decision 2026-10-06. |
| 24 | K24 | 3D coupled consolidation | 1 | 0.05 | UNKNOWN (high relevance) | `large2`: pore pressure differs 4.5 %, displacements 4 %. The 1D consolidation (`validation_2`) matches to 1e-11. |
| — | K99 | Unclassified | 10 | 0.23 | UNKNOWN | `axisym2`/`axisym5`/`force17`/`relax2`/`silo1`/`slide2-4`/`strain_volume1`/`interface_bar2_hex8` |

**K04: Pro analytic targets the GNU misses** (GNU value vs Professional value
and target):

| test | quantity | GNU | Pro / target |
|---|---|---|---|
| `condif7` | T(0), 1D edge heat flux; analytic T = 10 − x | **20** (flux doubled) → 10 after the 2026-10-06 fix (F1 T3: the `_node` restriction now selects whole sides, `area.cc` `side_node_restriction_ok()`) | 10 |
| `force18` | `force_point` at 0.8 | **`velx` −0.8** → +0.8 after the 2026-10-06 fix (F1 T3: the force was subtracted from `node_rhside`) | +0.8 |
| `elasti7` | K0 = 0.1 → `sigxx` | **0** → −0.0500000000025 after the 2026-09-29 fix (§5.3) | −0.05 (Pro −0.04999999884) |
| `tension_direct1` | `sigyy` with `…_normal_automatic` | **0.667** → 0 after the 2026-09-30 fix (new `element_normal` record; `disy` middle nodes 0.42 → 0, top 1; every variable = Pro except the nodal `eptyy` recovery 0.667 vs 0.5, K11; that residual is gone since F1 T5, 2026-10-06) | 0 |
| `mohr_coul_direct8` | `sigxy` with `plasti_bounda_factor` 0 | **1788.9** → 0 (PASS on the 2026-10-06 re-run, fixed by the C14/C16 work) | 0 |
| `mohr_coul_direct4` | `disx` | **0.3007** (GNU = manual 6.632 strength c + |Fn| tan phi; the Pro stays elastic at t = 1 and drifts below the strength by t = 100: deliberate divergence, user decision 2026-10-06, BACKLOG C30) | 0.1024 |
| `stress_pressure_factor` | `sigxx` | **−3.10** → −3.05 after the 2026-10-06 fix (F1 T3: factored trial pressure, strict tie, = Pro rule) | −3.05 |
| `incnav3` | `tvely` (`expansion_volume` effect) | **−2.49** → −4.9904 after the 2026-10-06 fix (F1 T3: rho(1 − beta T) in the inertia too; sfnet copy target −2.5 stale, BACKLOG B11) | −5 |
| `tet10_2` | `sigzz` | **−0.894** (TET10 nodal quadrature only degree 1, `polynom.cc`; needs design, moved to F4 (design + corpus A/B), BACKLOG C31) | −1 |
| `ground9` | interface groundflow `pres` | **0.333** → 0.5 after the 2026-10-06 fix (F1 T3: interface permeability/capacity/conductivity per unit length, capacity implicit) | 0.5 |
| `axisym9` | `post_point_dof 10 -sigxx` | 0.00158902 (uniform field; 1-point axisymmetric quad4 rule = K20, moved to F4 (design + corpus A/B), not the `post_point_dof` convention; BACKLOG C32) | 0.0015915 (±1e-6) |
| `distri2` | `post_point_dof 1 -temp` | 0.998933 (different random generator, both valid N(1, 2e-3); not a defect, BACKLOG C32) | 0.999774 (±1e-4) |

The `post_point_dof` of recovered fields follows a different convention in the
two binaries. It accounts for part of the 20 tests whose `node_dof` agree while
the `post_point` records differ. In `validation_2` the Pro reports the
element-local recovered value, taken from `element_dof`: −0.036441. The GNU
interpolates the averaged nodal field: −0.036558. Both are legitimate, but
targets on `post_point_dof` of stress or strain inherit the difference.

Known backlog items seen in this audit:

- **C5** is broader than documented: see K05.
- **C3** is inside K12.
- **The pre-existing Mohr-Coulomb return-mapping difference** is inside K15.
- **C7 (Mazars gamma, damage + plasticity)** is not exercised by any corpus
  test. The 7 damage tests match after `4cfaa0e`.
- **The `-updated` strain 5.99823e-4 case** did not show up as a separate
  cluster.
- **The force7 crash on the Pro** did not reproduce: `force7` ran and matched
  on both binaries.

### 5.1 K01 resolution (2026-09-29)

Two defects, both fixed in the working tree on top of `5568fb3` (commit
pending):

1. **Spurious stress rotation.** `materi.cc` "rotate to new configuration"
   listed `-UPDATED_LINEAR` (since `22801c3`) and `-UPDATED_AREA` (since
   `a1d9d12`) and rotated the stress by `inc_rot`. `set_deften_etc()` sets
   `inc_rot` to the identity for the linear memories and then overwrites it
   with `set_deften_u_rot(inc_deften, inc_u, inc_rot)`, so the stress was
   rotated by the real incremental rotation. Manual Pro 2.2.1:
   `-updated_linear` adds the extra stresses to the previous ones with no
   rotation step. Both memories removed from the list.
2. **Compressibility increment.** `stress.cc` used `tr(new_epe - old_epe)`
   for the linear memories; `old_epe` is 0 unless `materi_strain_elasti` is
   initialized, so the total strain was added every step (legacy 2014). Now
   `tr(inc_epe)`.

All 13 K01 tests re-run on the Professional, the pre-fix GNU (`5568fb3`) and
the fixed GNU, same converted text, every variable compared with the audit
machinery. Worst node per variable:

| test | variable | Pro | GNU before | GNU after |
|---|---|---|---|---|
| Pro `interac1` / sfnet `interac1` | sigxx n1 | -0.005 | -0.01 | -0.005 |
| | eptxx n3 | -0.005 | 9.7e-11 | -0.005 |
| | sigyy/sigzz n3 | -0.005 | -0.01 | -0.005 |
| Pro `matrix7` / sfnet `matrix7` | sigxy | 1 | 0.26593 | 1 |
| | sigyy | 0 | -2.0802 | 0 |
| | velx n5 | 0.5 | 0.31837 | 0.5 |
| sfnet `druckpr1` | sigxx / sigxy / sigyy | 0 / 0.11547 / 0 | 0.092376 / 0.069282 / -0.092376 | 0 / 0.11547 / 0 |
| sfnet `mohrcou1` | sigxx / sigxy | 0 / 0.1 | 0.08 / 0.06 | 0 / 0.1 |
| sfnet `tria6` | sigyy n8 | 1 | 1.8211 | 1 |
| | vely n5 | 0.5 | 0.67788 | 0.5 (all 13 nodes match) |
| vs `mshf` | sigxy / sigxx | 50 / 0 | 48.312 / 7.0268 | 50 / 0 |
| vs `mshf_nof` | sigxy / sigxx | 33.333 / 0 | 33.384 / 4.3608 | 33.333 / 0 |
| vs `mdp_shear`, `mstrain_druckprag` | sigxy | 34.641 | 34.145 | 34.640 |
| | sigxx | 34.631 | 45.824 | 34.654 |
| | sigyy | 34.643 | 34.031 | 34.538 |
| vs `mstrain_druckprag_elast` | sigxx | 0 | 823.50 | 0 |
| | sigxy | 666.67 | 356.56 | 666.93 |
| | eptyy | 0.66667 | 1.1746 | 0.66335 |
| vs `mesh_act_grav` | sigxx | 0 | 850.60 | 0.11 |
| | sigxy | 200 | -104.89 | 211.94 |
| | sigyy | 600 | 368.66 | 439.48 |
| | vel x/y | 1.2 / 0.2 | 2.0199 / 1.3988 | 1.2 / 0.2 |

- **9 of 13 now match in every `node_dof` and `post_point_dof` variable**
  (`interac1` x2, `matrix7` x2, `druckpr1`, `mohrcou1`, `tria6`, `mshf`,
  `mshf_nof`). Only exception: `interac1` `velx` at the middle node,
  1.2e-7 (GNU) vs 5.6e-14 (Pro), a solver-tolerance artifact on a value that
  should be 0 (pre-fix it was -5.1e-8).
- **The residual of the 4 other validation-suite tests is not the memory.**
  They all use `mesh -fixed_in_space` (convected). The fixed GNU with
  `-updated_linear` and the fixed GNU with the original
  `-updated_without_rotation` give the same dump for all four, so the
  remaining difference (0.07–0.3% on `mdp_shear`/`mstrain_druckprag*`,
  larger on `mesh_act_grav`) lives in the fixed-in-space/convection path
  and belongs to another cluster.
- **Other memories.** `-updated_area`: same fix; the Professional rejects
  it in an ordinary `group_materi_memory` (only accepted with
  `incremental_driver`), so there is no oracle. `-total_linear` +
  compressibility: the fix makes the GNU analytic (-0.005 for every dt)
  while the Professional is time-step dependent (-0.0075 / -0.0099 /
  -0.009999 for dt 1e-2 / 1e-3 / 1e-4, the pre-fix GNU digits): deliberate
  divergence, BACKLOG C9. `-updated`, `-total`: unchanged.
- **Reference vs current gradient.** With an explicit `mesh
  -follow_material` the Professional builds the `-updated_linear` increment
  on the current geometry (1D bar +50%: Pro 0.40714 = its `-updated` value;
  GNU 0.5 via the `c945833` reference gradient). The match without a mesh
  record comes from the different mesh defaults (Pro `-fixed_in_space`, GNU
  `-follow_material`), which also explains the `interac1` `-updated` residual
  (-0.0050117 vs -0.005). Recorded as BACKLOG C8, not changed.
- **Regression.** validation-suite 233 PASS / 1 RUNFAIL, sfnet 187/10,
  test-2026 17/17 (with the 2 new self-checking inputs
  `updated_linear_compress1.dat`, `updated_linear_shear1.dat`);
  classification identical to the pre-fix binary; every `.dbs` identical
  after masking uninitialized `element_dof` denormals and timestamps. The
  GNU-dialect corpora use `-updated_without_rotation`, so the fix only shows
  through the converted text. The conversion-neutrality count of section 1.1
  (32 changed approximations) was not re-measured.

### 5.2 K22 resolution (2026-09-29)

`ground18` is `ground4` with `groundflow_phreatic_level_static -yes` and no
`groundflow_pressure`. The stress dofs hold the effective stress, so the
analytic bottom `sigyy` is −1.5·10·50 − (2.0−1.0)·10·50 = −1250 (the
comment of `ground4` calls it total, but its numbers are effective). The
nodal value is the bottom element centre value, −1247.5.

Two gaps in `materi.cc`, both fixed in the working tree on top of `886d582`
(commit pending):

1. `get_materi_density()` chose wet/dry from `NODE_PHREATICLEVEL`, which is
   written only with a pressure dof, so every element was wet. The rule of
   the pressure-dof branch (average nodal total pressure of
   `groundflow_phreatic_coord()` < 0) now also serves the static groups.
2. The pore pressure was added to the force stress only with a pressure
   dof. A new branch adds the static total pressure at the integration
   point (times `groundflow_pressure_factor`).

Both apply to groups with `-groundflow` in `group_type`
(`groundflow_static_pressure_group()` in `groundfl.cc`), as the Pro does.

| test | variable | analytic | Pro | GNU before | GNU after |
|---|---|---|---|---|---|
| Pro `ground18` | `sigyy` bottom | −1250 (−1247.5 discrete) | −1247.5 | −1995 | −1247.5 |
| | `sigyy` at the level | −750 (−749.375) | −749.375 | −1000 | −749.375 |
| | `disy` top | — | −0.06875 | −0.1 | −0.06875 |
| | total pressure bottom | −500 | −500 | −500 | −500 |
| Pro `ground4` (pressure dof) | `sigyy` bottom | −1250 | −1247.5 | −1247.5 | −1247.5 |
| test-2026 `phreatic_static1` | `sigyy` bottom / level | −1247.5 / −749.375 | same | −1995 / −1000 | same as Pro |
| test-2026 `phreatic_static_submerged1` | `sigyy` bottom / y=−50 | −997.5 / −500 | same | −2095 / −1100 | same as Pro |
| Pro `excavate1` | `vely` node 13 (target) | — | −0.00210987 | −0.00210987 | −0.00253891 |

- `ground18` after the fix equals the Pro in every `node_dof` and
  `node_dof_calcul` value (worst `sigyy` difference 1e-9).
- **`excavate1`: new divergence (BACKLOG C10).** The Pro applies the
  static pore pressure only when `post_calcul -groundflow_pressure
  -total_pressure` is present. `excavate1` has no such record, so the Pro
  (and the pre-fix GNU) carry no pore pressure there. With the record
  added, the Pro gives −0.00253891, the GNU value. Other probes on
  `ground18`: without the record Pro −1745 (density only), GNU −1247.5.
  The GNU follows manual 6.574; the decision is left to the project owner.
- `ground8` (multiple levels with a pressure dof) is not touched by this
  fix: its `.dbs` is identical before and after. Its 19 % `sigyy`
  difference belongs to the multiple-level coupling
  (`groundflow_multilevel_mechanics.md`).
- **Regression.** Clean build: suite 16/16, the same 70 file-check failures
  by name. validation-suite 233 PASS / 1 RUNFAIL, sfnet 187/10, test-2026
  19/19 (the 17 earlier inputs + `phreatic_static1`,
  `phreatic_static_submerged1`); every `.dbs` identical to the pre-fix
  binary after masking timestamps and denormals (sfnet `examp14` times out
  at 45 s in both). No corpus test uses the static switch without a
  pressure dof.

### 5.3 K04 `elasti7` resolution (2026-09-29)

`elasti7` sets `control_materi_elasti_k0 20 -yes`. The hook in
`set_stress()` was guarded by `db_active_index( CONTROL_MATERI_ELASTI_K0,
0, ... )`, which tests that exact index, so the switch was ignored. The
record is now read at the ICONTROL index (fix in the working tree on top
of `28fde5d`, commit pending). Probes on the Professional 02-08-2026 also
showed that K0 applies without `group_materi_elasti_poisson` and to
`group_materi_elasti_hardsoil` (nu50/nuur); both are implemented.

| test | variable | Pro | GNU before | GNU after |
|---|---|---|---|---|
| Pro `elasti7` | `post_point` `sigxx` = `sigzz` | −0.04999999884 | 0 | −0.0500000000025 |
| | `post_point` `sigyy` | −0.4999999999521 | −0.499999999959 | −0.500000000025 |
| | top `disy` (nodes 3, 4) | −4.909090913e-7 | −5e-7 | −4.90909090909e-7 |
| | `node_rhside` `velx` | −0.025 | ~1e-17 | −0.0250000000012 |
| K0 without poisson (young / young_power) | `sigxx` | −0.05 | 0 | −0.05 |
| hardsoil nu50 = nuur = 0.2 + K0 | `sigxx` | −0.05 | −0.125 | −0.05 |
| K0 = 2 (clamp 0.95) | `sigxx` | −0.475 | 0 | −0.475 |

- After the fix every `elasti7` variable agrees with the Pro within 2.4e-8
  relative.
- The same literal-index-0 pattern was reviewed at every
  `db_active_index( CONTROL_xxx, 0, ... )` site:
  `control_materi_plasti_hardsoil_gammap_initial` (`top.cc`, crash for a
  non-zero index) and `control_solver_bicg_stop` (`so_bicg.cc`, the Pro
  uses index 0 only and lets it overrule `solver_bicg_stop`) were fixed;
  `control_materi_dynamic` in `so.cc` is not the cause of its divergence
  (BACKLOG C11); the other sites use internal records declared
  `no_index` (index 0 by design).
- **Regression.** Clean build: suite 16/16, the same 70 file-check failures
  by name. validation-suite 233 PASS / 1 RUNFAIL, sfnet 187/10, test-2026
  23/23 (19 + `elasti_k0_index1`, `elasti_k0_nopoisson1`,
  `elasti_k0_hardsoil1`, `hardsoil_gammap_initial_index1`, all rc=1 on
  the pre-fix binary); every other `.dbs` identical to the pre-fix binary
  after masking timestamps, `time_calculation` and denormals. The Pro
  `heat_exchanger_pile_dt` / `_load_dt` (K0 at index 45) time out at 45 s
  on both binaries, so they are not compared.

## 6. Recommended next fixes (priority by judgement)

1. **K01, `-updated_linear` in the GNU.** DONE 2026-09-29 (section 5.1).
   - Use the analytic `interac1` bar as the regression test: 2 materials, −0.005.
   - Start from the difference between the `-UPDATED_LINEAR` and
     `-UPDATED_WITHOUT_ROTATION` paths in `materi.cc`/`stress.cc`, with the
     compressibility material first.
   - This affects any model on the Pro's recommended memory.
2. **K22, phreatic level without a pressure dof.** DONE 2026-09-29 (section 5.2).
   - The wet/dry density split in the initial gravity stress; use the analytic
     −1250 of `ground4`/`ground18` as the check.
   - This is basic geotechnical initial state.
3. **K04, the Pro analytic target failures.** These are cheap, isolated bugs
   with exact answers, in order:
   - `condif_heat_edge_normal` doubled with `_geometry` + `_node`
   - `force_point` sign
   - `group_materi_elasti_k0` in the `-total_linear` + gravity path —
     DONE 2026-09-29 (§5.3: the switch was read for control index 0 only)
   - `plasti_tension_direct_normal_automatic`
   - `group_materi_plasti_bounda_factor`
   - `mohr_coul_direct4` — deliberate divergence (GNU = manual 6.632, BACKLOG C30)
   - `stress_pressure_history` switch
   - `expansion_volume`
4. **K05, contact.**
   - Move the `contact_target_geometry` alias out of the parallel loop and
     handle every index (C5).
   - Investigate why the 1D `-follow_material` contact (`contac1/2/4`,
     `axisym3`) transfers no momentum.
5. **K06, power-law elastic moduli.**
   - Decide the semantics against the Pro binary: the Pro adds the base
     `young`/`poisson`.
   - The current GNU choice was justified by a claim the binary refutes.
   **2026-10-06 (F2 T1):** `poisson_power` fixed (the Pro replaces nu;
   the GNU added both stiffnesses). `young_power` + `young`: the Pro adds
   both stiffnesses, against its manual and its own history variable;
   parked for an owner decision (BACKLOG C36).
6. **K03.** Align the `convection_apply` default to `-no` (Pro manual 6.395),
   the same deliberate alignment as `inertia_apply` in `0560dcf`. Audit the
   corpus inputs that relied on the old default.
   **Done 2026-10-06 (F1 T2):** default `-no` via `convection_apply_switch()`
   (general.cc, one source for general/beam/truss/spring/conspr). The Pro
   applies the switch only to MATERIAL convection wrt the mesh: the
   `group_condif_flow` and groundflow-velocity heat transport stay on (sfnet
   condif1/condif5 give temp 1 in the Pro with `-no`, `-yes` or absent), and
   the GNU now does the same. condif2 now equals the Pro (0.5). Inputs that
   relied on the old default: validation-suite `cmat_gate`, `cd_method`,
   `materi_direct_wall` got an explicit `convection_apply -yes`; sfnet
   `condif2/3`, `examp1`, `incnav5`, `strconv1` are BACKLOG B10.
7. **K11 and K02, output state.**
   - Recover the strain dofs with the same operator as the stress.
     **Done 2026-10-06 (F1 T5, `nodal_strain_recovery.md`).**
   - Extend the post-solve sync to every non-principal dof: strain, plastic
     strain, and models with a strain unknown. **Parked 2026-10-06
     (BACKLOG C35): the integration-point state lags too; design item,
     scheduled for F4 (project decision 2026-10-06).**
   - Cheap, and it removes about 25 tests of noise from future audits.
8. **K13.** Store `materi_strain_plasti` and `plasti_f` for the direct and visco
   cut-off models. **Done 2026-10-06 (F1 T4)**, see
   `manual-developer/direct_plastic_strain.md`; follow-ups C33 (visco
   time step, scheduled for F4) and C34 (shear convention, deliberate
   divergence, project decision 2026-10-06).
9. **K17.** Truss / beam / `support_edge_normal` divergence on the Pro inputs
   (`trubea2/3`, `truss2/3/14`, `support_edge_normal_3`).
10. **Investigate before touching code.** These are UNKNOWN with high
    relevance:
    - K14: `hypo3`/`hypo10` stress path, and the empty hypo history slots
    - K15: `validation_14` DP, `validation_3` von Mises
    - K24: `large2` 3D consolidation

## 7. What makes this audit less reliable

- **Half of the GNU corpora cannot run on the Pro.** 207 tests run only on the
  GNU, mostly because of GNU extensions and GNU-dialect records. The
  reconstructed validation-suite is the worst case: 110 of 233 run on the GNU
  only. The features the GNU added itself are therefore the least cross-checked.
- **The conversion is not always neutral.**
  - `-updated_without_rotation`→`-updated_linear` changes the GNU result in
    32 tests; this is what exposed K01.
  - Dropping `number_of_integration_points` and `materi_velocity_integrated`
    are approximations, flagged per test in `summary.json`.
  - The comparison itself always uses the same text on both binaries.
- **No integration-point comparison.** The GNU does not write
  `element_intpnt_dof`, and its `element_dof` layout differs from the Pro's
  (integration-point values vs nodal values). The audit therefore compares
  nodal recovered fields. A recovery-operator difference (K11, K20) looks like
  a constitutive difference.
- **Almost no time-history comparison.** The Pro writes no `.his` files; its
  histories go into GiD `.res` files, which the runner deletes. Only the
  `incremental_driver` histories (3 tests, 1000 steps each) were compared.
  - Drained oedometer and triaxial: MC match to 1e-4.
  - Undrained: 0.04 %.
- **Timeouts hide the large geotechnical models.** The slope, foundation, dam,
  Vajont, pile and consolidation examples do not finish in 180 s on either
  binary.
- **Node pairing.** When the refined-mesh numbering differs, nodes are paired by
  coordinates: 43 tests. Four of them are only partly paired: `disk1`
  (101/105), `validation_14` (561/580), and `delete4`/`examp15` (467/507). In 3
  tests (`elasti6`, `hex20`, `quad8`) the node counts differ; the common nodes
  are compared by number with identical coordinates.
- **Target evaluation is partial.** `-post_node_result`, `-post_line_dof` and
  `element_intpnt_*` items are "NA".
- **Clusters are rule-based.** Hand lists cover the representatives, and 10
  tests remain unclassified. Magnitudes are relative to the family scale, so a
  small absolute difference in a nearly zero family can rank high (K23).
- **The Pro cannot be inspected.** For "the Pro ignores X" (K16, K10) the audit
  cannot tell a Pro defect from a different Pro syntax without the Pro source
  or manual examples.
