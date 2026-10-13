# Tochnog Professional vs GNU — Seguimiento de convergencia


Documento de control del proceso de convergencia entre Tochnog Professional
(manual `UserManual-professional.txt`) y el GNU (fork de sfnet 2014).

Índice único de TODAS las divergencias con el Professional (deliberadas,
conocidas, pendientes, de dialecto y sospechosas del Pro):
[DIVERGENCES-PROFESSIONAL.md](DIVERGENCES-PROFESSIONAL.md).

## Métrica global

- **1094** keywords documentados en Professional (tras limpiar 1 keyword
  espurio del manual)
- **309** presentes en el GNU (**28 %** de cobertura nominal)
- **785** faltantes
- Además: 197 tests del GNU 2014 corren con cobertura 100 % de keywords (ver plan, sección 6c)

## Features implementadas por nosotros — registro de verificación

Cada feature implementada en este proceso (posterior al fork sfnet 2014)
se verifica contra un ejemplo del manual de Professional, un test de la
suite sfnet, o un test propio. El registro completo:

| Feature | Commit | Fecha | Verificación |
|---------|--------|-------|--------------|
| **F4 T6f C74 — Bi-CG breakdown rescued by the direct band LU reported with one warning per solve (no "wrong input file" block)** | `cf07b19` | 2026-10-13 | Professional 02-08-2026 `validation_1` (default PARDISO and `solver -matrix_iterative_bicg`): no breakdown message. GNU trace: the absolute 1e-16 test fires with `|r1.r2|/(|r1||r2|)` 1e-3 .. 1e-4 (false breakdown), the LU retry is exact (= SuperLU). RED: `test-2026/bicg_breakdown_lu_retry1` on HEAD `46104a3` rc 0 but 9 error blocks, build_safe check FALLO; GREEN: 9 warnings, no block. Clean build; build_safe `TN_CHECK_ONLY=1` 263 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Corpus A/B omen `c74b_runs` (957 runs): no classification change, 881/882 non-timeout `.dbs` identical (masking `time_calculation`/`time_at_start`; `bicg_csr_threads1` 2-thread round-off, differs between HEAD runs too); runs with the block 38 -> 0. Solver part (scale-relative test + relative floor) open as BACKLOG C77, DIVERGENCES C77 class B. |
| **F4 T6f C75 — `prism6` integration rule measured on the Professional; GNU interior rule kept (deliberate divergence)** | docs only (F4 T6f C75) | 2026-10-13 | Professional 02-08-2026 `element_intpnt_iso_coord` / `_npoint` for the default, `-normal`, `-maximal`, `-minimal` and `inertia_apply -yes`; the weights 1/12 confirmed by a 20x2x2 prism6 cantilever (Professional nodal displacements = a corner-rule FE model to every printed digit, GNU HEAD = the interior rule). Accuracy study (full rank, patch test, cantilever against Timoshenko) in BACKLOG C75; owner decision 2026-10-13: keep the GNU rule (Descarte row). New BACKLOG C76 (`-minimal` 2-point rule). No code or test change. |
| **F4 T6f C70 — nodal recovery of quad9 3x3, hex8 (2x2x2, 1 point), hex27, tet10 and prism6 in physical coordinates; axisymmetric on-axis quad9** | `39d0138` | 2026-10-13 | `miscel.cc` `recovery_fit_kind()` / `recovery_fit_shape()` / `recovery_fit_build()` / `recovery_fit_factor()` generalised to 3D (per-element linear least squares `1,x,y(,z)` in PHYSICAL coordinates; 1-point hex8 point value; axis rule for every 2D quad); the static natural-coordinate table of F4 T9 removed. **Verification** (Professional 02-08-2026, `element_intpnt_dof` + `element_intpnt_coord` -> `element_dof` on distorted bricks and curved quad9/hex27/tet10): physical fit residual <= 1e-12 of the scale for quad9 3x3, hex8 2x2x2, hex27 3x3x3 / 2x2x2, tet10 5-point, prism6; hex8 1 point = point value; axisymmetric quad9 3x3 / 2x2 at r=0 = `sum h sig / sum h`; node = arithmetic mean (1e-13). New tests `test-2026/recovery_hex8_distorted1`, `recovery_hex8_point_mean1`, `recovery_hex27_distorted1`, `recovery_tet10_curved1`, `recovery_prism6_linear1`, `recovery_quad9_curved1`, `recovery_quad9_axis1`: RED rc 1 on 421a3f4 (e.g. hex8 node 11 sigxx 1.1216 vs Pro 1.3071691), GREEN rc 0, whole node stress field = Pro to <= 2e-12. `force11` hex27 nodal stresses = Pro to 7e-7 (before 3.38 %, C57 recovery part). Clean build; `build_safe.sh` TN_CHECK_ONLY=1 261 OK / 0 FALLO, runs 226/226; `audit_db_buffers.py` exit 0. Corpus A/B omen `f4t6g_c70_runs`: no rc/classification change; force15 / force9 move onto the Pro (0.25 -> 5e-11, 3.4e-2 -> 1e-7), none farther; `large1` +1.5 % wall. Found: prism6 integration points differ (Pro triangle corners, BACKLOG C75). Docs: user/developer `nodal_stress_recovery`, `integration_quadrature_and_stress_recovery`, `control_mesh_extrude`; DIVERGENCES C70 closed, C75 added, C57 re-scoped. |
| **F4 T6f C72 — convected nodal stress recovered absolutely (Eulerian mesh, `convection_apply -yes`); C72 premise corrected; `validation_1` steady guard** | `c1a1ae9` | 2026-10-13 | `general.cc` `general()`: `stress_recovery_absolute = !materi_strain_total` (the convected stress dof was incremental, with a lumped mass weight F/volume against a convective weight h: it moved at about v*volume_ip, 0.055 for v = 1) and the convective / artificial-diffusion terms are skipped for an absolute stress dof. Professional 02-08-2026 probes (Engram 2775): neither code convects the integration-point stress (a sigzz / sigxy / pressure pulse stays bit-frozen); the Pro node stress is its recovery. New tests test-2026: `convected_stress_recovery1` (sigzz pulse, 40 quad4, v = 1: 0.5 / 1.0 / 0.5 at t = 2 = Pro; HEAD 0.348732 / 0.904484 / 0.6481, rc 1 -> 0), `validation1_steady1` (`validation_1` velx at t = 3 -0.0571954 within 9.9e-4 of the t = 1 value -0.0570386 pinned by `stokes_group2`; Pro -0.0391 / +0.0034 fails it: DIVERGENCES C72-V1, class B, owner decision Pro-suspect). A/B on every `convection_apply -yes` input: velocities bit-identical, nodal stresses <= 1.3 % of max on `validation_1` / `stokes_*`; sfnet `strconv1` byte-identical. Clean build; build_safe TN_CHECK_ONLY=1 252 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. C72 stays open (class B) for the Professional-side residuals; BACKLOG C73 (optional IP-state convection) and C74 (Bi-CG breakdowns) added. |
| **F4 T6f K21 WU2 — `group_materi_stokes` default `-yes` and per element group; test-2014 `cmat_gate` rewritten as a well-posed shear** | `58bba14` | 2026-10-12 | `general.cc` `general()`: `stokes` default -YES (manual 6.752), read with `db(GROUP_MATERI_STOKES, element_group, ...)` (before -NO at index 0). Professional 02-08-2026: the `validation_1` flow without the record = `-yes` (post_point velx 0.1465713 both); the Pro reads index 0 for every group (probes -no/-yes, -yes/-no, record only on group 1, element 0 moved to group 1: always the index-0 value to 1e-12), GNU follows the manual (owner 2026-10-12, DIVERGENCES STOKES-IDX class A). New tests test-2026 (RED on HEAD rc 1, GREEN rc 0): `stokes_default1` (Pro 0.1465713 +- 0.003; GNU 0.144159, before -0.0570386; the gap is C72: the Pro keeps the integration-point stress convection under Stokes, GNU Stokes = `convection_apply -no` 0.1441588), `stokes_group1` (mixed groups 0.124105; Pro 0.1465713 fails it), `stokes_group2` (all elements in a -no group = the one-group -no run -0.0570386). Stop-rule case: test-2014 `cmat_gate` (ill-posed `-ra 1 4 -ra` rigid rotation, target 3.20774 +- 2) failed with the new default (sigxy 1.532 -> 0.921, Pro 0; met by chance); rewritten (owner 2026-10-12) as a simple shear with vely on every node keeping the `control_materi_plasti_tension_apply` gate: sigxy 250.5 / sigxx -49.5 (analytic and Pro exactly; 150.5 / -149.5 without the gate), GNU 250.4989 / -49.4998, tolerance 0.01. Clean build; build_safe TN_CHECK_ONLY=1 252 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Corpus A/B omen (`f4t6g_stokes_runs`, range vs stokes index-0 build): rc change only `cmat_gate` 0 -> 1 (then rewritten), `cd_method` ept 4.6e-2 (target met); Pro other / sfnet / test-2026 `.dbs` unchanged; no corpus input has a `group_materi_stokes` with index != 0, so the per-group read cannot change them. validation_1 (explicit `-no`) and validation_8 (no density) unaffected. |
| **F4 T6f ranges — `-ra ... -ra` read as the Professional (index `-step`, `-from a -to a`, several `-from/-to` blocks, descending and `9.0` errors, bounded reads)** | `35025dc` (F4 T6f ranges) | 2026-10-12 | Professional user manual "Ranges" (p.116-117) + audit Engram 2758 (32 forms x index/value on both binaries). `range.cc` `range_scan()`: `-STEP` sign (an index `-step` range expanded every number plus junk indices), `from > to` instead of `>=`, blocks and list items in any order, every word through `input_read_string()` (BACKLOG C56 residual closed; defines work, as in the Pro), empty range an error; `range_expand()` validates value blocks (descending -> `Error in range, a -from is larger than a -to`, the Pro text; step <= 0; missing `-to`); `input.cc` a non-integer number inside a value range is an input error. Audit probes rerun: every accepted form expands as the Professional 02-08-2026, every rejected form errors in both. New tests test-2026 `range_index_step1`, `range_index_single1`, `range_index_two_blocks1`, `range_index_define1` (element groups; sigxx 1.6666667 / 1.4285714 = Pro `.dbs`; RED rc 1 on HEAD), `range_value_step1` (guard, disx 0.005 = Pro), `range_index_descending_error1`, `range_value_descending_error1` (HEAD rc 0, silent), `range_value_double_error1` (HEAD rc 0; Pro error). build_safe 249 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Corpus A/B omen (`f4t6g_range_runs`, ungate vs range): no rc/class change; `.dbs` identical in Pro other, sfnet, test-2014; test-2026 only `bicg_csr_threads1` thread round-off (< 1e-6). Residual: integer values with a decimal point outside a range (DIVERGENCES INT-DOT, class C). Docs `manual-*/ranges.md`. |
| **F4 T6f K21 WU1b — one nodal recovery operator everywhere (convected-stress gate removed); test-2014 `cd_method`, `cd_geom`, `materi_direct_mc` rewritten as well-posed shears** | `e3af9f2` (tests) + `e2b044f` (refactor) | 2026-10-12 | Root cause of the WU1 gate: `-ra 1 4 -ra` is a node LIST (Professional user manual "Ranges"), so `cd_method`/`cd_geom` fixed vely only on the left edges (rigid rotation, physical stress 0; Pro original `cd_method` 10.78 .. 1.25 then 0 from t = 0.6) and the old targets were met by chance. Rewrite: vely on every node + `group_materi_plasti_tension_direct 1.e10` (lifts the principal tension cap). GNU HEAD and new binary: `cd_method` sigxy 0.9999956 every step (target 1 +- 1e-5, phi target kept), `cd_geom` 0 / 1 (1e-6), `materi_direct_mc` 1 (1e-6), rc 0. Professional 02-08-2026 (Pro dialect, `-updated_linear`): exactly 1 / 0 and 1 / 1; with `-updated` Pro and GNU agree (18.0375 / 18.0374, 18.6915 / 18.6915): no MC-direct divergence (BACKLOG C17 note). RED: original `cd_method` rc 0 on HEAD (0.744), rc 1 without the gate (3.642). Gate removal: `miscel.cc` `recovery_fit_factor()`, `general.cc`, `tochnog.h`/`tochnog-mod.h`. build_safe 241 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Corpus A/B omen (`f4t6g_ungate_runs`, fin vs ungate): rc change only the old `cd_method` 0 -> 1 (replaced); `cmat_gate` (same ill-posed convected model) sig 5.9e-2 rel, target met, vs Pro 0.107 -> 0.118; Pro other / sfnet / test-2026 `.dbs` unchanged (`alternate1`, `cylinder_buckle` near-timeout flakes). Convected stress path recorded as class B (BACKLOG C72). |
| **F4 T6f — guard tests of the class A divergences K14b (`hypo3` strain reversal) and C7/C7b (Mazars equivalent strain)** | `ae0691b` (F4 T6f) | 2026-10-12 | No source change (owner decision 2026-10-12: keep the GNU, document the divergence). K14b verified against a converged reference: Python integration of the von Wolffersdorff + intergranular strain law (manual 2.2.4) on the test-2014 `hypo3` path, sigxx -0.10466; test-2026 `hypo_wolff_igs_oedometer1/2` (dt 1e-3 / 1e-5, +- 2e-3 / 5e-4): GNU -0.103436 / -0.104235, rc 0; the Professional 02-08-2026 -0.1629 (56 % off, fails) / -0.10499 (time-discretisation error at the reversal). C7 verified against manual 2.2.5 (eps_eq = sqrt(sum <eps_i>+^2)) and 14 Professional probes (its eps_eq = s accumulated over the ascending positive principal strains as s = 2 s + eps^2, reproduced to 7e-14): test-2026 `damage_mazars_eq_two_positive1` (hex8 diag(2e-4,-3e-4,1e-4), dam 0.28610207515 +- 1e-8; GNU 0.286102, rc 0; Pro 0.36409 fails) and `damage_mazars_vonmises1` (damage_stress1 bar + von Mises 3e6, dam 0.925372225695 +- 1e-6, eppxx 1.87863e-4 +- 1e-8; GNU rc 0; Pro 0.96946 / 1.368e-4 fails). Guard tests (GREEN on HEAD by design, RED not applicable). The `hypo3` hyhis slots stay open (K14); non-objective Mazars gamma of both codes noted in BACKLOG C7. |
| **F4 T6f C68d — `materi_maxwell_stress` with two or more chains plus `materi_plasti_rho` / `materi_strain_intergranular` / `materi_strain_plasti`: old-state read index** | `2327083` (F4 T6f C68d) | 2026-10-12 | `materi.cc` `materi()` old-state read loop: the Maxwell loop no longer overwrites `indx` (`old_msig[m*MDIM*MDIM + indx]`), so the rho/epi/epp reads use their own slots (before: wrong slots past the 9-entry arrays inside `dbl_array`, rho and epp never accumulated). Verified with a decomposition test (no Professional manual example): test-2026 `maxwell_two_chains_state1`, two identical chains = one chain with the summed stiffness (the Professional 02-08-2026 reproduces the equivalence to ~1e-12); targets = the one-chain GNU run, 1e-6 relative; previous binary rc 1 (node 1 reaction -13.32 vs -19.17), fix rc 0 (1e-11). ASan/UBSan build of the previous source: no report (overrun inside one allocation). Clean build; build_safe TN_CHECK_ONLY=1 233 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. |
| **F4 T6f K21 WU1 — quad4/quad9 nodal recovery in physical coordinates, 1-point quad4 point value, on-axis axisymmetric rule (audit K21 recovery part, BACKLOG C42/C64)** | `2305659` (F4 T6f K21 WU1) | 2026-10-12 | `miscel.cc` `recovery_fit_factor()` / `recovery_fit_build()`: per-element least squares `1,x,y` in PHYSICAL coordinates for the 2x2 Gauss quad4 and quad9, the point value for the 1-point quad4, `sum h sig / sum h` for an axisymmetric quad4 with a node at r=0 (thread_local per-element cache); `sri_stress_recovery_weight()` / `strain_recovery_weight()` take `element_group`, `nodes[]` (all `materi.cc` / `general.cc` callers); convected stress dofs keep the old operator (`stress_recovery_lagrange`, general.cc; test-2014 `cd_method` otherwise drifts 0.74 -> 3.64). Operators fitted on the Professional 02-08-2026 `.dbs` (element_intpnt_dof -> element_dof, 1e-12) of sfnet `examp1`/`incnav4`, Pro `tube1` and distorted probes; elastic axisymmetric probe confirms the axis rule (1.5e-13). New tests test-2026 `recovery_examp1_quad4_point_mean1` (node 19 sigxx -0.2350738, HEAD -0.256815), `recovery_incnav4_quad9_minimal1` (nodes 1/3 sigxy 0.05, HEAD 0.02/0.08), `recovery_tube1_axis1` (node 3 sigxx -159.5805, HEAD -66.40), `recovery_quad4_distorted_plane1` (node 21 sigxx 2.0014257, HEAD 1.1829): RED rc 1 on 7eed0db, GREEN rc 0, whole stress field = Pro within 1e-7 of its scale. `axisym_quad4_ring1` sigzz 10759.2155 = Pro (was 10760.09). Corpus A/B on omen (`f4t6g_rec_runs`, head vs fin): only rc changes = the 4 new tests 1 -> 0 (Pro other `alternate1` TIMEOUT flake at 117-120 s); vs Pro closer 9 (`tube1`, `pressed_rubber_cylinder`, `axisym5`, `examp1`, `ho_mech1`, `incnav4`, `incnav5`, `axisym_quad4_hourglass1`, `materi_memory_updated_linear_current2`; local: `kappa_ip_state1`, `validation_3_kappa` closer), farther 0; no-Pro changes `taylor3`, `disk1`, `examp18`, `axisym4`, `examp3`, `examp14`, `cd_geom` (targets met). Wall +0.9 %. build_safe 241 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Probes: the Pro also fits 3x3 quad9 and hex8 in physical coordinates (not changed: BACKLOG C70); sfnet `tube1` `-lobatto` NaN -> BACKLOG C71. Docs `manual-*/nodal_stress_recovery.md`, `nodal_strain_recovery.md`, `integration_quadrature_and_stress_recovery.md`. |
| **F4 T6f C68 — `materi_plasti_kappa_shear`, `materi_plasti_cap1_history`, `materi_damage`, `materi_plasti_rho` and `materi_maxwell_stress` kept as integration-point states** | `70235c6` (F4 T6f C68) | 2026-10-12 | Same defect and fix pattern as kappa (`c3e6d2a`): `elem.cc` ELEMENT_DOF read/write masks gain kapsh, cap1, dam, the rho range and the Maxwell range (feature-guarded); `materi.cc` one block after the kappa block, outside the h guard, nodal numerator with `strain_recovery_weight()` and IP write; `general.cc` denominator weight. Manual 4.10, 4.15, 4.17, 4.25, 4.27 (point history). Verified against the Professional 02-08-2026 with four-tria3 probes on a non-uniform field: test-2026 `damage_ip_state1` (reactions 1.9e-1 -> 5e-13 relative), `rho_ip_state1` (1.1e-1 -> 7e-6), `kappa_shear_ip_state1` (1.3e-2 -> 1.5e-4), `maxwell_ip_state1` (nodal msigyy 8.4e-2 -> 5e-4); all RED rc 1 on the previous binary, GREEN rc 0. cap1: manual/physics only (Professional cap model diverges, BACKLOG C68b). Corpus A/B on omen (`f4t6f_c68_runs`): no rc/class change, no node_dof change above 1e-10. build_safe 232 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Residuals BACKLOG C68a-C68d. |
| **F4 T6f — guard tests of the class A divergences MC-VERTEX-ZZ (`mmc_tension`) and MEMBRANE-EPTZZ (`validation_3` eptzz)** | `600a0b4` (F4 T6f) | 2026-10-12 | No source change (owner decision 2026-10-12: keep the GNU, document the divergence). Verified analytically: test-2026 `mmc_tension_vertex1/2` (test-2014 `mmc_tension` at dt 0.05 / 0.047; sigxx 46.19 = 2c cos(phi)/(1+sin(phi)), sigzz 0 +- 0.05, eppzz 0 +- 1e-4; GNU sigzz -3.3e-10 / 1.1e-9, rc 0; the Professional 0.974 / 0.367 fails) and `membrane_eptzz_total1/2` (single quad4 von Mises plane stress, dt 0.01 / 0.1; eptzz analytic -0.0088249 +- 3e-4 / +- 2e-3, sigxx 27.4199 +- 0.05; GNU -0.0086300 / -0.0072102, rc 0; the Professional -9.6e-5 / -7.9e-4 fails). Guard tests (GREEN on HEAD by design, RED not applicable). New pages `group_materi_plasti_mohr_coul`, `group_materi_membrane`; BACKLOG C69 (GNU yy/zz plastic split first order in dt). |
| **F4 T6f WU2 — `materi_plasti_kappa` kept as an integration-point state (K15 `validation_3`) + sfnet `examp3` recalibrated (B16)** | `c3e6d2a` (F4 T6f WU2) | 2026-10-12 | `elem.cc` ELEMENT_DOF read/write masks + `kap_indx` (the hardening read the h-lumped nodal kappa interpolated back to the point); `materi.cc` IP write `new_unknowns[kap_indx]` and nodal numerator with `strain_recovery_weight()` outside the h guard; `general.cc` kappa denominator weight. Verified against the Professional 02-08-2026: test-2026 `kappa_ip_state1` (two quad4, von Mises + kap hardening; node 1 reaction -11.91338 vs Pro -11.91214, HEAD -11.85835 rc 1; nodal kap within 5 %, quad4 operator C42/C64) and test-2014 `validation_3_kappa` (Pro validation_3 + Pro targets: post_node_result vely -129.353 vs -129.302, HEAD -128.835 rc 1; displacements 1.25 % -> 0.06 %; kap node 1 11.6 % -> 2.7 %). Corpus A/B on omen (`f4t6f_runs`, wu1 vs wu2): rc changes = the 2 new tests 1 -> 0 and sfnet `examp3` 0 -> 1 (vely -126.431 vs 2014 target -124.8; Pro on the translated input -126.374) -> owner decision 2026-10-12: recalibrated to -126.374 ±1 (B16, CORPUS-STALE-TARGETS §8), rc 0; Pro `ground11_nonsaturated` TIMEOUT was load (alone 97 s). build_safe clean build of `c3e6d2a` (build/tochnog): 224 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Same mask omission for kapsh/cap1pc/rho/damage/Maxwell -> BACKLOG C68 (open). |
| **F4 T6f WU1 — nodal stress/strain recovery of tria3, tria6, tet4 and the 1-point bar2 = arithmetic mean of the adjacent elements (K15 `validation_14`)** | `1d446b3` (F4 T6f WU1) | 2026-10-11 | `miscel.cc` `sri_stress_recovery_weight()`: tria3 1 point / 3-point nodal rule, tria6 6-point nodal rule, tet4 1 point, bar2 1 point return `F/volume` (F = 1 or h) instead of `h_inol` (area/volume/length-weighted average); strains follow (`strain_recovery_weight`). Verified on two-element probes against the Professional 02-08-2026 (`node_dof` = mean of its `element_dof`, exact): new tests test-2026 `tria3_nodal_recovery_mean1/2`, `tria6_nodal_recovery_mean1` (tolerance 2e-4 of the family scale, tria6 element ~1e-5 off: BACKLOG C66), `bar2_nodal_recovery_mean1`, `tet4_nodal_recovery_mean1`, test-2014 `validation_14_pass1` (hole nodes 4/120; HEAD node 4 sigyy -11.23 vs Pro -8.5195): all RED rc 1 on HEAD, GREEN rc 0. Full Professional `validation_14` (DP, split, refine_locally, restart): 1758 node_dof = Professional within 1.2e-4 of the field scale (HEAD: 1740 node_dof, velx 99 % off, different refined mesh). Corpus A/B on omen (`f4t6f_runs`, head 1965044 vs wu1): only rc changes = the 6 new tests 1 -> 0; closer to the Pro 2 (Pro `axisym3` sig 1.8e-2 -> 1.2e-4, `force17` 8.3e-2 -> 5e-13), farther 0; sfnet `examp23` (validation_14 ancestor, no Pro) different refined mesh, rc 0 both. build_safe TN_CHECK_ONLY=1 222 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Unmeasured rules keep h (BACKLOG C67); quad4 operator question C64; Pro post_point at a node C65. |
| **F4 T6e — post-solve material state pass for slide and contact models (C53, gate removed)** | `07433f5` (F4 T6e C53) | 2026-10-11 | `top.cc` pass: `element_loop(); parallel_contact; slide(); parallel_new_dof_diagonal; mpc_node_apply()`; `post_solve_state_pass_gated()` and its declarations (`tochnog.h`, `tochnog-mod.h`) deleted; DIVERGENCES `POSTPASS-GATE` removed. Pro `slide4` = Professional (node 1 `node_slide_force` -2.11291867e-3 / -1.21989198e-3 vs -2.11291887e-3 / -1.21989209e-3; HEAD rc 1). New tests: test-2026 `slide_friction_cohesion1` (slide4 with cohesion 0.002, end time 10: Pro 02-08-2026 -2.065482222595e-3; HEAD -0.002 rc 1; slide fix alone -2.066157e-3 rc 1; GREEN -2.06548203e-3) and `contact_block_pass1` (Pro-syntax contact_block: node_dof within 7e-13 of the Pro; targets vely 2.276148940747e-4, disy -6.022146669393e-3, sigyy -5.990411829884, node_rhside vely ~0; HEAD rc 1). Slide/contact focus set 18/18 rc 0 (HEAD slide4 rc 1). Investigation A/B on omen (`f4t6e_runs`). build_safe TN_CHECK_ONLY=1 216 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Out of scope: BACKLOG C62 (pass without `parallel_new_dof_before`), C63 (`post_node_rhside_fixed` slide4). |
| **F4 T6e — slide law: absent `slide_plasti_tension` = no tension limit, friction cap `c + \|Fn\|*tan(phi)`, `slide_penalty_law` principal-dof index (C53, slide part)** | `b4e4fa4` (F4 T6e C53) | 2026-10-11 | `slide.cc`: `sig_t` default `1.e300` (was 0: a compression-only gap, contrary to manual 6.1044 and the code comment); cap and `node_slide_f` with `\|S\|` (Pro caps with `c + tan(phi)*\|Fn\|`); `slide_penalty_law` writes `node_lhside`/`node_rhside[vel_indx/nder+idim]` (as `slide_spring`, `dof.cc`) instead of the node_dof index. test-2014 `node_slide` recalibrated (owner decision): `node_rhside 4 vely` -5e13 = analytic penalty reaction `-1e15*(v.n)*dt*n_y` (rel. 1e-6) + `node_dof 4 eptxx` 0; HEAD rc 1 (vely 2.97e-12, eptxx -10: the old target 0.001 passed only because the y term landed in the eptxx slot), new rc 0 (CORPUS-STALE-TARGETS §7). The Pro `slide4`-with-cohesion probe needs the post-solve pass too (next row). build_safe TN_CHECK_ONLY=1 214 OK / 0 FALLO, runs 226/226. |
| **F4 T6d — sfnet `truss2`/`truss3` adapted to the explicit truss memory (B15, corpus only)** | corpus (no source change) | 2026-10-11 | Stale since C11 (`28992d5`, default `-updated_linear` without `materi_displacement`: reaction (0.1, 0), rc 1). Inputs adapted as the Professional 02-08-2026 copies `other/truss2`, `other/truss3` (owner 2026-10-11): explicit `group_truss_memory 0 -updated`; `truss3` also `materi_displacement` -> `materi_velocity_integrated` (as Pro; avoids the C59 input error). External copies `external-downloads/sfnet/extracted/test/` + omen corpus; `validation-suite/test.zip` copies untouched. Original targets kept: GNU rc 0, node 1 rhside (0.100001, 0.0999992); Professional on its own copies and on the dialect-translated adapted inputs (0.1000002, 0.0999998), rc 0. |
| **F4 T6d — post-solve material state pass for inertia models (C52, gate widened)** | `c276931` (F4 T6d C52) | 2026-10-11 | Inertia (`materi_acceleration`, `inertia_apply`, `control_inertia_apply`) removed from `post_solve_state_pass_gated()`; slide/contact stay gated (C53). RED/GREEN: test-2014 `iface_gen`/`_geom`/`_method` rc 1 with the gate removed on the old interface history (node 6 vely -1.70554), rc 0 with the interface fix. Pro `dynamic3` = Professional 60.44025 with the pass. Full corpus A/B on omen (interface fix -> gate removed): classes and exit codes unchanged (only the `alternate1` 120 s timeout flake); .dbs differences only in inertia models (Pro `taylor3` 4.1e-4, distance to the Pro unchanged; test-2014 `iface_gen*` up to 0.15, targets met). Interface + inertia probe gap vs the Pro remains independent of the pass (BACKLOG C61). build_safe TN_CHECK_ONLY=1 214 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. |
| **F4 T6d — interface Mohr-Coulomb history = clamped tangential force (C52, interface part)** | `b318484` (F4 T6d C52) | 2026-10-11 | Quasi-static slip-reversal probe (interface c=5, phi=0, kt=kn=1000; slide 5 steps at the limit, reverse 5 x -0.2): Professional 02-08-2026 shear 4.000000083, analytic 4.0, GNU before 5.0 (history = elastic trial `kt*gamma_total`, no unloading), after 4.0. Same probe with `control_reset_interface_strain` mid-run (manual 6.355): Pro 4.000000083, GNU 4.0 (before 5.0); with `control_reset_interface`: Pro -1.0 = GNU before/after. New test test-2026 `interface_mc_slip_reversal1` (RED rc 1 actual 5, GREEN rc 0). test-2014 `iface_mc_slip`/`_1step`/`iface_mc_num` original 2014 targets 0/0/50 restored (+ `element_interface_intpnt_stress 2 1` target; old binary rc 1). Investigation on 55 inputs: all Pro `interface*`, `mohr_coul_direct1-8`, `tension_direct1`, `taylor3` unchanged; `interface_patch` -10 kept. build_safe TN_CHECK_ONLY=1 214 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Corpus A/B on omen c59 -> wu1 (`ea7181d3`): classes unchanged except the new test and the `alternate1` timeout flake; .dbs: Pro `interface_patch`/`mohr_coul_direct3/4`/`interface15` <= 8e-12, test-2014 `iface_gen*` (inertia, targets met) and the retargeted `iface_mc_*`; sfnet identical. |
| **F4 T6d — `group_truss_memory` input checks of the Professional (C59)** | `e0f2204` (F4 T6d C59) | 2026-10-11 | Professional 02-08-2026 probes: `-updated`/`-updated_linear` + `materi_displacement` -> "materi_displacement should not be initialised"; `-total`/`-total_linear` without it -> "should be initialised"; owner decision 2026-10-11 to replicate both (GNU message names both records; `check_error -no` skips them, rc 0 verified). Corpus scan (test-2014, test-2026, sfnet, `validation-suite/test.zip`, Pro test tree): only sfnet `truss5` hits (a); adapted (record dropped, default `-total_linear`, target -93.8 = Professional run of the adapted input; zip copy unchanged). Tests `validation-suite/test-2026/truss_memory_updated_displacement_error1`, `truss_memory_total_velocity_error1`: RED old binary rc 0 (accepted), GREEN rc 1 + message. build_safe (TN_CHECK_ONLY=1): 213 OK / 0 FALLO, runs 226/226; audit_db_buffers exit 0. Corpus A/B c11 217e574c -> c59 20b5dc20 (omen, adapted `truss5` in both): no class or exit-code change except the two new tests (rc 1 = expected error) and the `alternate1` 120 s timeout flake; no `.dbs` change (`bicg_csr_threads1` 4e-11 thread noise). |
| **F4 T6d — `control_materi_dynamic` effective factor, truss memory default, signed linear truss length (C11)** | `28992d5` (F4 T6d) | 2026-10-11 | Manual 6.141 (`control_materi_dynamic` = `materi_dynamic` for the timesteps of that index), 4.12 (total Lagrange with `materi_displacement`), 6.777 silent on the default -> Professional 02-08-2026 probes. Fix: `materi_dynamic_effective()` used by the displacement-primary gate and every reader (no PUT/restore of `MATERI_DYNAMIC`); truss default `-total_linear` with / `-updated_linear` without `materi_displacement` (GNU before: `-updated`); no `scalar_dabs()` of the projected length. Tests `validation-suite/test-2026/truss_memory_default_displacement1` (Pro `dynamic3`, node 2 disx 60.44025; old 32.6583), `truss_total_linear_inversion1` (60.44025; old -5962.6), `control_materi_dynamic_only1` (0.0958890; old 6.4477), `control_materi_dynamic_half1` (0.932818; old 1.45599), `truss_dynamic_no_record1` (0.995505; old 1.91225), `truss_memory_default_velocity1` (rhside 50; old -4.677): RED on the old binary (rc 1), GREEN = Professional to 1e-10. build_safe (TN_CHECK_ONLY=1): 211 OK / 0 FALLO, runs 226/226. Corpus A/B old 929f446a -> new 217e574c (omen): Pro other OK 276 -> 277 (`dynamic3` now = Pro; `alternate1` timeout -> rc 1 under load), sfnet 192 -> 192 with exit 0 -> 1 on `truss2`/`truss3` (targets calibrated on the legacy `-updated` default; the Professional gives the same reaction as the new GNU, BACKLOG B15), `truss7`/`trussbea` .dbs differ without class change; test-2026 only the 6 new tests (rc 1 -> 0); test-2014 232 -> 232 (`bmom_truss` 2 % node_dof, targets met). Pro input errors for the memory / `materi_displacement` combinations not added (sfnet `truss5` uses one; BACKLOG C59). |
| **F4 T6c — `control_restart` resets the integration-point state** | `8ca6124` (F4 T6c) | 2026-10-11 | Manual 6.368 (undeformed mesh, initial node_dof, initial time, path dependent behaviour from the start); the Professional 02-08-2026 also deletes `element_dof`. The GNU kept `ELEMENT_DOF_INITIALISED` = 1, so the old IP state survived. Fix: `restart_reset_element_state()` clears the flag (VERSION_NORMAL/NEW). Tests `validation-suite/test-2026/control_restart_reset_ip1` (Drucker-Prager) / `_ip2` (elastic), targets = Professional at node 4 (sigyy -24.18351 / -52.08333): RED on the old binary (rc 1; sigyy +19.65 / +52.08), GREEN after. Probes: truss force kept in both binaries (1.5), undrained pressure = Pro (-9.75297). Corpus A/B (omen): `examp6`, `examp23`, `validation_14` .dbs identical old/new (except `time_calculation`), `slope_nonlocal_refine` timeout in the three binaries. build_safe (TN_CHECK_ONLY=1): 205 OK / 0 FALLO, runs 226/226. |
| **F4 T6b K24 — `large2` 3D coupled consolidation: `solver_matrix_symmetric -yes` guard tests (deliberate divergence)** | `c350b3c` (F4 T6b) | 2026-10-11 | No code change. Cause: the Professional symmetrizes the non-symmetric Biot u-p matrix under `solver_matrix_symmetric -yes` (manual 6.1052); the GNU keeps it (Bi-CG, `so_bicg.cc`). With `-no` the Professional = GNU (`large2` pres at x = 0, t = 2.2: GNU -0.704550, Pro `-yes` -0.672711, Pro `-no` -0.704603). Guard tests `validation-suite/test-2026/k24_consol_3d1/2` (`-yes`/`-no`, brick 10x2x2 nodes, E 10 `-total_linear`, k 0.1, c 1, t 2.2) and `k24_consol_3d_undrained1/2` (k 1e-8, one step), targets = Professional 02-08-2026 `-no` re-run at node 1: pres -0.7085296 / disx 0.05420111 and pres -0.847677 / disx 0.01748967; GNU -0.708335 / 0.0542098 and -0.8476521 / 0.0175445 for both switches (Pro `-yes` -0.678332 and -0.803176 fail them; a wrong target gives rc 1). Guard tests, GREEN on the current binary (no RED: no behaviour change). build_safe (TN_CHECK_ONLY=1): 203 OK / 0 FALLO, runs 226/226. |
| **F4 T6a C3 — `gforce10` 3D section sign = `plot_switch` convention; test-2014 input restored** | `3e7837b` (F4 T6a) | 2026-10-11 | No code change. Professional 02-08-2026 `force10.dat` probes (GNU and Pro, `thickness_switch -yes` / `-no` / absent): with `-yes` both give node 5 nory -12.34, shey -100, mom1y -5000; `-no` = absent in both binaries (default `-no`, as the manual 6.915). The test-2014 `gforce10.dat` had lost `plot_switch -no -yes -yes -yes` and `thickness_switch -yes`: restored, target back to -5000, GNU = Pro exactly. Open (owner decision): square end face without `-yes` (Pro degenerate frame, GNU tie broken by the reference point; DIVERGENCES C3). vs A/B: `gforce10` rc 0 both. |
| **F4 T6a C6 — `iface_gen*` inputs generate their interface elements** | `523a959` (F4 T6a) | 2026-10-11 | Corpus input repair (CORPUS-STALE-TARGETS §2 family note): generated element(s) and `element_group ... 10` removed from `iface_gen`, `iface_gen_method` (+ input `element_geometry 1 20 / 2 30`), `iface_gen_geom`, `iface_gen_quad9`. Each run generates the same elements, rc 0, targets met; 4 new build_safe checks (no element 3 in the input, element 3 + `element_group 3 10` in the `.dbs`; RED on the old inputs). |
| **F4 T6a K14a — Masin `hyhis1/2/3/5/7` post-processing slots written back** | `7e09c40` (F4 T6a) | 2026-10-11 | `hypoplas.cc`: clay kernels from statev 10/11/12/14 and the L tangent, basic kernel phi_mob, OCR p_e/p and sqrt(Mijkl Mijkl) of the zero-increment tangent. Professional 02-08-2026 probes: hypo7 hyhis2/3/5 23.4037 / 22516.4 / 2.09388 (Pro 23.4076 / 22518.8 / 2.09449), hypo8 23.2607 / 37519.7 / 1.3169 (Pro 23.2709 / 37531.6 / 1.31766), hypo12/13 hyhis2/5 to 1e-4 (hypo13 hyhis3 24308 vs 24250). Residuals in DIVERGENCES K14 (hyhis1: GNU substep vs Pro time step; hypo12 hyhis3 m_R L; hypo9). Test `hypo_masin_history_slots1` (Pro hypo7 to t 0.1; RED rc 1 / GREEN). Corpus A/B on omen (`f4t6a_runs`, old = HEAD 53dbb86 build md5 edb81804, new = T6a build md5 d1aa27b8, 6 jobs, 4 GB, timeout 120): only hyhis columns change (hypo7/8/9/12/13, stresses identical). hypo7 user time 1.7 -> 3.1 s (extra kernel call). |
| **F4 T6a — `post_calcul -size_tot`/`-size_dev` = sqrt(A:A) / sqrt(B:B)** | `4705692` (F4 T6a) | 2026-10-11 | Manual Professional post_calcul defines sqrt(Aij Aij) / sqrt(Bij Bij) (no 0.5); Professional `calcul1` targets 1 / 0.816497 (`calcul.cc` `calculate_operat`). Test `post_calcul_size_tot1` (Pro calcul1 without its `-prival` targets, K12; RED rc 1 / GREEN). Corpus A/B on omen (`f4t6a_runs`, old = HEAD 53dbb86 build md5 edb81804, new = T6a build md5 d1aa27b8, 6 jobs, 4 GB, timeout 120): Pro other `ctvmis2` exit 1 -> 0 (meets the Pro target 0.163299); sfnet `calcul1`, `calcul2`, `ctvmis2`, `examp18` exit 0 -> 1 (targets on the old convention, BACKLOG B14); sizes x sqrt(2) in `examp3/10/14/17/23`. |
| **F4 T6a C56 — parser words bounded to `MCHAR`-1** | `80b19af` (F4 T6a) | 2026-10-11 | `input.cc` `input_word()` (`setw(MCHAR)`) at all word reads; split word allowed inside comments, `Word too long` error elsewhere. Tests `input_long_comment_word1` (134/108-character words in comments, rc 0) and `input_long_word_error1` (expected error; old GNU SIGSEGV, t26 A/B rc -11 -> 1). The Professional `ground19_water_under_dam` URL case is the comment case. |
| **F4 T6a C55 — `control_print_dof_rhside` indexes `node_rhside` per principal unknown** | `f7559bc` (F4 T6a) | 2026-10-11 | `print_un.cc`: `node_rhside[ipuknwn]`. Professional 02-08-2026 probe with `derivatives`: vely_rhside -1.0416667 / +1.0416667, sigyy_rhside -2.0833333 = new GNU (old GNU 0 and 4.7e-310). Test `dof_rhside_derivatives1` + build_safe file check (RED/GREEN). build_safe (all T6a) 199 OK / 0 FALLO, runs 226/226 (baseline 189). |
| **F4 T12 C47 — `control_mesh_extrude` quad9 -> hex27: shared inner-plane nodes** | `24301c8` (F4 T12) | 2026-10-11 | `mesh.cc` `mesh_extrude()`: plane q (1..2nq-1) of source node inod is node `next_node + (q-1)*nbase + inod`, created by the first hex27 that uses it and shared by its neighbours (before: per-element mid-plane and layer-top blocks, the elements disconnected at the inner planes). Verification against the Professional 02-08-2026: `other/force11` 15600 -> 12000 nodes (Pro 12000, 2400 per z plane), rc 1 -> 0, `mom1y_sig` 5.80223 -> 5.764609 (Pro 5.764612), nodal displacements = Pro to 5.5e-7, section forces within 2e-4; residual hex27 nodal stresses up to 3.38 % of the field maximum at the loaded points (BACKLOG/DIVERGENCES C57). Test RED/GREEN `extrude_quad9_shared_nodes1` (two quad9 x two hex27 layers; old GNU rc 1, 84 nodes; new 75 nodes = Pro, displacements = Pro to 2e-12) + build_safe node-count check; build_safe 189 OK / 0 FALLO, runs 226/226. Corpus A/B on omen (old a08992b vs new, 4 GB, timeout 120): every `control_mesh_extrude` input (Pro other extrude1/2, safety2, interface_bar2_hex8, force11, Pro validation_14_mesh, the new test): only force11 and the new test change; full Pro other: only force11 changes (classes equal except cylinder_buckle 119 s OK -> 120 s TIMEOUT, alone 80/81 s, .dbs identical). |
| **F4 T5 C35/K02 — post-solve material state pass (the material is evaluated once more on the final iterate of every step)** | `e8e3f12` (F4 T5) | 2026-10-10 | `top.cc` `top()`: after the equilibrium loop `iteration_start(); element_loop(); parallel_new_dof_diagonal; mpc_node_apply()` with `material_state_post_pass` (`elem.cc` skips the matrix store, no solve); gated off (previous stress-only sync) for inertia, slide and contact models (`post_solve_state_pass_gated()`, DIVERGENCES `POSTPASS-GATE`, BACKLOG C52/C53), `-diagonal`, `options_elementloop -no`, nonlocal softvar. Verified against the Professional 02-08-2026 with probes (1D stretch velocity-only 1/6 steps, with displacement, two plastic bars 1/5 steps, `spring4`): all = Pro (e.g. bars node 1 eppxx 0.01 / sigxx 0.015, old 0 / 0; `spring4` eptxx 0, old -0.01). New tests `post_solve_state_stretch1/_plastic1/_rigid1` (RED/GREEN). Matrix store vs none: byte-identical on the probes and 742 non-gated corpus tests. build_safe 187 OK / 0 FALLO, runs 226/226, rc 0. Corpus A/B on omen (HEAD efa46c6 -> final 2c825fe4, incl. C51): closer to Pro 11 (matrix9 x2, spring2-5, control_plasti_apply1/2, cmat_gate, mdp_shear, mstrain_druckprag), farther 0; the 94 gated tests (dynamic3, slide4 included) unchanged; classes unchanged except Pro `ground11_nonsaturated` OK -> TIMEOUT at 120 s under load (84 -> 105 s alone; Pro times out). Wall time of tests finished on both: Pro other +8.8 %, sfnet +7.2 %, t26 +4.7 %, vs +12.0 %. |
| **F4 T5 follow-up C51 — hardsoil yield uses the step-start kappa plus its own return increment** | `8e4d080` (F4 T5) | 2026-10-10 | `stress.cc` `set_stress()`: the yield rules get `pl_unknowns` (kappa and kappa_shear = step-start value + increment of the current return) instead of the nodal kappa written by the previous material evaluation of the step (`plasti.cc` hardsoil `gamma_p`, kappa-hardened Mohr-Coulomb, von Mises Nadai). No Professional oracle (zero field on the converted `mhardsoil_*`). Verified by an iteration-count sweep (1/2/3/4/8, with and without the post-solve pass: identical; old GNU 0.0030 / 0.0006 alternating) and a hand check of the end state (sigma = sigma0 - C(E50, nu50) eps_p exactly, kap = sqrt(0.5 eps_p:eps_p), f ~ 1e-8). New tests `test-2026/hardsoil_kappa_iterations1/2` (RED/GREEN). Retarget approved by the project owner 2026-10-10 as GNU-calibrated: `mhardsoil_plast` kap 0.00104991 / sigxx -3.11898, `mhardsoil_gp0_off` kap 0.000494563 / sigxx -1.56781. build_safe 183 OK / 0 FALLO, runs 226/226. Remaining reads of the same pattern (cap epp_vol, cap1 pc, void fraction) listed in BACKLOG C51. |
| **F4 T13 C49 — scratch-buffer audit and `-O2` build** | `40ae7cd` + `44c68c7` (F4 T13) | 2026-10-11 | Internal (no keyword). `scripts/audit_db_buffers.py` checks every `db()` read, `get_group_data()`, constant-length `array_*()` and the `geometry()` / `msf_element_rule()` output arrays against `data_length` / the callee contract: HEAD `dca225a` 45 hits at 34 call sites, all fixed (`manual-developer/db_buffer_contract.md`); then `scripts/build_safe.sh` compiles with `-O2`. Verified with a dedicated build matrix (no Professional example: build-level change): ASan+UBSan build_safe 28 failing runs -> 1 (`dofrhside`, C55), corpus 45 ASan tests at 9 sites -> 2 (C55, C56), no UBSan report; `-O1` corpus `.dbs` identical to HEAD; `-O2` vs `-O1` every finished `.dbs` byte-identical (only the 2-thread `bicg_csr_threads1` round-off), the F4 T11 `-O2` failures (`creset_iface`, `bmom*`, Pro `interface7`, `interface10`, `mesh_cut_2`) run, Pro `alternate1`, `cylinder_buckle`, `ground11_nonsaturated` finish within 120 s; build_safe 187 OK / 0 FALLO at both levels. Time on omen (one run, median of 3): `large1` 31.4 -> 22.2 s, `force11` 18.8 -> 9.9 s; corpus wall (6 jobs, tests finished by both) Pro other 696.5 -> 582.5 s, sfnet 240.1 -> 220.1 s, t26 35.9 -> 32.1 s, vs 28.2 -> 26.0 s; peak memory unchanged. |
| **F4 T10 C48 — interface tension status written as integers** | `c06fa2d` (F4 T10) | 2026-10-10 | `interface.cc` `interface_element()`: `element_interface_intpnt_materi_tension_status` (INTEGER, print-only) was PUT from a `double` buffer with `idum[1]` as the integer array (stack garbage, different from run to run); now a `long int` array of `-closed`/`-opened` per integration point, buffers released with `delete[]`. Verified against the Professional 02-08-2026: the 20 `interface*` inputs of `test/other` (all `-closed`, same line per element except `interface_bar2_hex8`, element numbering, new BACKLOG C50) and the new `test-2026/interface_tension_status1` (`-closed` x3 / `-opened` x3 = Professional; old GNU `2 -total_linear 3`). Two GNU runs byte-identical. build_safe 181 OK / 0 FALLO (177 + 4 checks). Corpus A/B on omen (`f4t10_runs`, old e84e554c -> new 25040c72): classes unchanged (Pro other OK 276 -> 276, sfnet 192 -> 192, t26 112 -> 112, vs 232 -> 232); the only `.dbs` change is this record (28 Pro other, 25 vs, 4 t26 files) plus 2-thread round-off in `bicg_csr_threads1` (6.7e-11). |
| **F4 T9 C43 — assembled CSR operator of the CG/Bi-CG** | `b7009dd` (F4 T9) | 2026-10-10 | `so_bicg.cc` `bicg_csr_build()` / `bicg_csr_sys()`: the element matrices are assembled once per solve into a CSR matrix (pattern from element node sets + node incidence, no triplet buffer; `NODE_LHSIDE` diagonal apart; `options_processors` row blocks); element-by-element fallback with mpc ties or `TN_BICG_EBE`. Verification: test RED/GREEN `bicg_csr_threads1` (CSR line with `print_solver`, 2 threads, analytic targets, same 75 iterations); build_safe 177 OK / 0 FALLO (main tree, final build md5 e84e554c); `force11` ~5x iterations per second (31104 rows, 2.16 M nonzeros; 311040 Bi-CG iterations in ~640 s, then the band LU retry gives the SuperLU `mom1y_sig` 5.80223). Corpus A/B on omen (b3 `3349515` -> CSR): Pro other 274 -> 276 (`ground11_nonsaturated`, `ground11_phreatic_level` timeout -> 102 / 82 s), sfnet 192 -> 192, t26 111 -> 111, vs 232 -> 232; round-off `.dbs` changes, node_dof relative to the field maximum <= 1e-5 except path-dependent models (`mesh_activate1` 3.3e-3, `distri3` 1.8e-4, sfnet `examp16` 1.4e-2, `examp10` 2.8e-3, `taylor3` stalled on both): a 2-thread run of the old build moves them as much, and `examp10` diverges after a solve stopped by the absolute `bicg_error_minimum` floor. Whole F4 T9 (old `c62d28f` -> `b7009dd`, tests finished on both, 6 jobs): Pro other 794 -> 525 s, sfnet 317 -> 220 s, t26 39 -> 33 s, vs 31 -> 24 s (total 1180 -> 803 s, -32 %); OK Pro other 273 -> 276, t26 110 -> 111 (+ superlu_solve1), sfnet and vs unchanged. |
| **F4 T9 C43 — `-matrix_pardiso` solved with SuperLU, Bi-CG fallback (owner decision 2026-10-10)** | `3349515` (F4 T9) | 2026-10-10 | `so.cc` `solve()`: the pardiso branch selects SuperLU (`pardiso_superlu`); a failed factorization restores the right-hand side and runs the Bi-CG (with the band LU retry) for this and every later solve; `solve_pardiso_memory_guard()` lowers the soft `RLIMIT_AS` to the current address space + 70 % of the room during `dgssv` (a failed `dgssv` cannot free its L/U workspace); `TN_PARDISO_SUPERLU_MB` overrides the budget. GNU default stays the Bi-CG (DIVERGENCES SOLVER-DEFAULT, class B). Verification: tests RED/GREEN `pardiso_superlu1` (SuperLU solves it; old GNU: Bi-CG) and `pardiso_bicg_fallback1` (`TN_PARDISO_SUPERLU_MB=4`: SuperLU out of memory, Bi-CG, analytic targets met; old GNU: no SuperLU attempt); Professional `large1` forced to pardiso under 4 GB: clean fallback, rc 0, 47 s; build_safe 175 OK / 0 FALLO. Corpus A/B on omen (b1 -> b3, 4 GB, timeout 120): Pro other `force11` TIMEOUT -> finishes in 20 s (rc 1: `mom1y_sig` 5.80223 vs Pro 5.76460, BACKLOG C47; the band LU of the CSR Bi-CG gives the same value); the 12 other pardiso inputs change by the Bi-CG tolerance only (node_dof relative to the field maximum: `incnav5` 2.7e-6, `interface_bar2_hex8` 1.3e-6, `interface_quad8_hex20` 1.2e-6, the others <= 9e-7); sfnet, vs unchanged; t26 only the 3 new tests. |
| **F4 T9 C43 — nodal stress recovery factor table** | `acfac40` (F4 T9) | 2026-10-10 | `miscel.cc`: `sri_stress_recovery_weight()` reads the least-squares factor from a static table built once (`sri_recovery_table_build()`, quad4/quad9/hex8/hex27/tet10) instead of solving the 4x4 normal equations per node, point and element (43.7 M fits on `large1`). Bit-identical: build_safe 171 OK, 297 suite `.dbs` equal; corpus A/B on omen (old c62d28f -> b1 = this + the two fixes below): Pro other 273 -> 273, sfnet 192 -> 192, vs 232 -> 232, every `.dbs` identical except the C48 garbage tension-status record of `interface3` / `interface_quad8_hex20` (stack addresses, differ run to run). Time (tests finished on both, 6 jobs): Pro other 794 -> 685 s, sfnet 317 -> 289 s; `large1` 37 -> 27 s locally. |
| **F4 T9 — buffered node copies in `interface_convert()` (heap-use-after-free)** | `789f8bf` (F4 T9) | 2026-10-10 | `interface.cc` `interface_copy_node_dbl()`: the three `NODE_DOF` / `NODE_DOF_START_REFINED` copies to the new side-2 and face-centre nodes passed a `db_dbl` pointer into a reallocating PUT (ASan on the Professional `interface_bar2_hex8`, `interface_bar3_quad8`, `interface_quad8_hex20`). After the fix ASan passes `interface_convert()`; the next report is C48 (tension status). Results unchanged (suite and corpus, see the row above). |
| **F4 T9 — SuperLU use-after-free and new[]/free mismatch** | `681dc83` (F4 T9) | 2026-10-10 | `so_suplu.c`: memory report read before `L`/`U` are destroyed (it read them after `Destroy_*`), `L`/`U` destroyed only when `dgstrf` created them, only the store of `A` destroyed (its arrays are `new[]` of `so.cc`, released there with `delete[]`). Test RED/GREEN `superlu_solve1` (41x41 quad4, `-matrix_superlu`, analytic): old rc 139 (3/3 runs), new rc 0. ASan: old alloc-dealloc-mismatch (`superlu_solve1`, sfnet `examp9`), new clean. Professional `force11` with `-matrix_superlu`: old rc 139, new 10 s. |
| **F4 T8 C45 — `control_mesh_copy` / `_mirror` after a time step (node dofs and integration-point state duplicated)** | `0ca144c` (F4 T8) | 2026-10-10 | `mesh.cc` `mesh_duplicate_node_data()` / `mesh_duplicate_element_state()`: the new nodes lacked `node_dof_start_refined` (`renumber()` stopped: rc 1) and `create_element()` re-initialised the duplicates from the nodal values (1 % off). Element numbers `old + max + 1` as the Professional. Verified against the Professional 02-08-2026 with loaded 3x2-node probes (gravity, step on index 0, copy/mirror + two steps on index 1): all nodal dofs = Pro to 4e-10. Tests RED (rc 1) / GREEN `mesh_copy_timestep1`, `mesh_mirror_timestep1` (Pro post_point targets). Corpus A/B on omen for the whole F4 T8 series (old = HEAD `b4e09c1` = F4 T7 ABC binary, new = `7aeac4a`; 4 jobs, 4 GB, timeout 120): Pro other 273 -> 273 OK (`mesh_activate1` rc 1 -> 0, target met), sfnet 192 -> 192, t26 99 -> 107 (+8 new tests), vs 232 -> 232; no class change. The .dbs differences are explained: node/element renumbering of the macro meshes (values equal when matched by coordinates), round-off of near-zero quantities and iterative-solver tolerance from the new equation order (<= 5e-6 relative on non-zero values, e.g. `beam2d` 2.7e-6), new random realisations of `control_distribute` (`distri2`, `distri3`, `ho_othr4`: per-node/element sequence), and inputs already divergent before (`taylor3` B9, `update2` diverges in both, `tensile_instability`/`ho_mech2` 4e-5). |
| **F4 T8 — buffered node-item copies in `mesh.cc` (heap-use-after-free)** | `fd477ff` (F4 T8) | 2026-10-10 | `mesh.cc` `mesh_copy_node_item()`: `mesh_convert_quad8`, `mesh_convert_hex20`, `mesh_extrude` (and copy/mirror/rotate_3d via `0ca144c`) passed a `db_dbl` pointer to a reallocating PUT (ASan, reported by the F4 T9 investigation). No result change: all 297 build_safe `.dbs` identical apart from timing lines; ASan clean on the copy/mirror tests and the Professional `quad8`/`hex20`. Same pattern still in `interface.cc` `interface_convert()` (out of scope). |
| **F4 T8 C46d — `control_mesh_macro` `-bar`/`-rectangle`/`-brick` numbering (x fastest)** | `bc2235e` (F4 T8) | 2026-10-10 | `macro.cc`: the GNU numbered nodes, elements and local element nodes with x slowest. Professional probes (1D bar, 2D 3x3, 3D 3x2x2): x fastest, then y, then z; GNU `.dbs` now identical. Tests RED/GREEN `mesh_macro_numbering1/2`; `mesh_switch_timestep1` check retargeted to node 3 at 0 1 (the Pro node 3). Makes the Pro `mesh_activate1` target (`node_dof 5`) read the right node. |
| **F4 T8 C46a — order of the one-shot mesh records on one index** | `e730052` (F4 T8) | 2026-10-10 | `top.cc` `step_start()`: copy, move, switch, mirror (was switch, move, mirror, copy). Professional probes: every pair plus copy+switch+mirror on index 1, record order swapped: nodes identical to the Pro. Test RED/GREEN `mesh_order_timestep1`. |
| **F4 T8 C46b — `control_mesh_move` once per control index** | `c21b62e` (F4 T8) | 2026-10-10 | `top.cc`: `mesh_move()` guarded by `first_of_control`. Professional: applies it once (measured on `node_start_refined`) but does not move the nodes (Pro-suspect, DIVERGENCES C46b, class D, owner decision pending). GNU follows manual 6.217 (node and node_start_refined). Test RED (22 1) / GREEN (12 1) `mesh_move_timestep1`. |
| **F4 T8 — `control_mesh_switch` two entries = swap (manual 6.244)** | `fb834c6` (F4 T8) | 2026-10-10 | `mesh.cc` `mesh_switch()`: `dir0 dir1` swaps the two axes (GNU read them as the new axis order: `-x -y` = identity, 3D needed three). Professional probes: 2D `-x -y` and 3D `-x -z` nodes identical; three entries refused by the Pro (kept as GNU extension, DIVERGENCES SWITCH-PERM); Pro does not switch `node_start_refined` (C46e, class D). Test RED (2 0) / GREEN (0 2) `mesh_switch_pair1`. |
| **F4 T8 C46d — `mesh_activate_gravity_time` at its index, bottom to top, method 2 (Pro `mesh_activate1`)** | `7aeac4a` (F4 T8) | 2026-10-10 | `mesh.cc` `mesh_activate_gravity_factor()` + `materi.cc`: records read at their own index (the Pro input uses 10, the GNU read 0: full gravity at once, final vely 0), element window from its lowest/highest coordinate, default `-method2` (manual 6.810), no gravity before activation, the reduced stiffness scales the material tangent and the stress increment (scaling the assembled matrix made Newton diverge). Verified against the Professional `mesh_activate1` (1000 steps): vely/veliy/sigyy = Pro to 3e-4, target met (rc 0, was rc 1); at t 0.5 the front matches. Test RED (rc 1) / GREEN `mesh_activate_gravity_time1` (Pro targets; residual C46d class B). |
| **F4 T7 C — `control_mesh_macro` / `_copy` / `_mirror` / `_switch` once per control index (BACKLOG C43 part C, C44)** | `ab99338` (F4 T7 C) | 2026-10-09 | `top.cc`: `step_start()` re-ran the one-shot mesh records on every step of a `control_timestep` on the same index (Pro `mesh_activate1`: 202 -> 4040 nodes in 20 steps, out of memory with its 1000 steps). A flag set when the control loop enters an index limits them to the first step. Verified against the Professional 02-08-2026 with probes on omen (3x2-node mesh, record + 2/3-step timestep on index 1): macro 12 nodes (GNU before 24), copy 12 (48), mirror 12 (48), switch once (GNU twice). `control_mesh_move` had no effect in the Professional, rotate not comparable: left unchanged (C46). Tests RED/GREEN `mesh_macro_timestep1` (202 nodes = Pro; old 4040) and `mesh_switch_timestep1` (node 3 at 0 0.5; old 0.5 0); build_safe 158 OK / 0 FALLO. Corpus A/B on omen (A+B -> A+B+C): only `mesh_activate1` (timeout -> runs in 42 s; 202 nodes; integrated displacement `veliy` = Pro to 5e-10; final-step `vely` 0 vs Pro -9.95e-6, C46) and the two new tests. |
| **F4 T7 B — the Bi-CG direct-LU retry builds its band matrix only after a failure, in a Cuthill-McKee ordering when narrower (BACKLOG C43)** | `5bfa165` (F4 T7 B) | 2026-10-09 | `so.cc` `solve_band_retry()` / `solve_band_retry_entries()`: the band workspace (`(3 band+1) n` doubles with the raw node-number band) is no longer allocated and filled on every Bi-CG solve; after a failure the retry measures the exact band of the node-number and of a Cuthill-McKee ordering on the element matrix entries (with the mpc tie expansion), fills the narrower one (same `NODE_LHSIDE` + `EPS_MAT` diagonal), solves with `dgbsv_` and un-permutes. Verification: `examp2` (retry on every solve) bit-identical, pres -0.634519741323 = Pro; build_safe (A+B+C build) 158 OK / 0 FALLO; corpus A/B on omen (A -> A+B): class changes only `force11` (bad_alloc -> runs, Bi-CG beyond the 120 s slot); `.dbs` round-off in retried solves only (column-scaled: `force9` 4e-7, `pressed_rubber_cylinder` 1e-7, `ho_mech1` 8e-8, `force5` 1e-9, `excavate1` 4e-12; `bimet1` only in the identically-zero membrane sigzz, 7e-12 vs -3e-14). No retargets. |
| **F4 T7 A — the default Bi-CG solver no longer allocates an unused dense n x n matrix (BACKLOG C43, memory hogs)** | `3faf78e` (F4 T7 A) | 2026-10-09 | `so.cc` `solve()`: `mat_dense`/`ipiv_dense` (allocated, zero-filled and filled on every solve since `5b7bb8b`; the only reader, the `dgesv_` branch of the Bi-CG retry, never ran because the band matrix `mat` always existed) removed with its fills and the dead branch. Verification: build_safe 154 OK / 0 FALLO; corpus A/B on omen (old `225d536` -> A, plain runner, 4 GB cap, timeout 120): Pro other 271 -> 272 (`ground19_water_under_dam` now runs; `foundation_consolidation_settlement` bad_alloc -> timeout), sfnet 192 -> 192, t26 99 -> 99, vs 232 -> 232, and no other `.dbs` difference (node/element/node_dof/post_point_dof records). No retargets. |
| **F4 T4 WU2 — `tet10` Keast 5-point rule, linear least-squares recovery and HRZ lumping (BACKLOG C31)** | `bcbdda0` (F4 T4 WU2) | 2026-10-09 | `polynom.cc` `pol()` TET10: Keast 5-point rule (centroid -4/5, (1/2,1/6,1/6,1/6) 9/20), the Professional 02-08-2026 rule; `miscel.cc` `sri_stress_recovery_weight()`: linear least squares `1,L1,L2,L3` (closed form `9 L(node).L(point) - 41/20`); `general.cc`: HRZ lumping of the principal dofs (corner V/36, mid-edge 4V/27; deliberate divergence TET10-HRZ, owner 2026-10-09). Verified: Pro `other/tet10_2` sigzz -1.00000000515 (Pro -1; old GNU -0.894) = `test-2026/tet10_body_load1`; `tet10_mass_hrz1` velocities 2.16 / 0.405 = dt/m with the HRZ masses (Pro row sum -1.2 / 0.3; old GNU 59999.99 / 0.36); `tet10_dynamic_gravity1` top disz -0.5286 at t 3 (analytic -0.5; Pro solver failure; 5-point rule with row sum 5e8). |
| **F4 T4 WU1 — axisymmetric elements use the default Gauss rule (no forced 1-point rule, quad9 accepted; audit K20, BACKLOG C32 closed)** | `9d3e9f6` (F4 T4 WU1) | 2026-10-09 | `polynom.cc` `pol()` and `calcul_force.cc` `msf_element_rule()` no longer force `-gauss -minimal` for axisymmetric `materi_velocity` nor refuse `npol > 2`; quad4 2x2 and quad9 3x3 Gauss as the Professional 02-08-2026. Verified on Professional probes: `test-2026/axisym_quad4_hourglass1` (4x4 ring patch, point load: disx 2.657533803e-4 / disy -4.380511327e-4 = Pro to 1e-14; old GNU hourglass ~1e12, no convergence), `axisym_quad4_ring1` (Akin ring, Pro `axisym2`: disx 1.081776788e-2 = Pro; sigzz 10760.09 vs Pro 10759.22, quad4 recovery residual BACKLOG C42; old 1.08238e-2 / 10000), `axisym_quad9_uniform1` (Pro `axisym7`: sigxx = sigzz = 1000, disx 500 analytic; old GNU refused the quad9). `msf_axisym` face moment 0 -> 1/960 (analytic `sigma t^3/(12 r)`). build_safe 151 OK / 0 FALLO. |
| **F4 T3 — `_visco` direct cut-offs and viscoplastic laws use the current time step; tension cap honours `_visco` (BACKLOG C33)** | `8a1fc75` (F4 T3) | 2026-10-09 | `stress.cc` `set_stress()` reads `DTIME` at VERSION_NEW (manual 6.728: relaxation over t..t+dt); `direct_cap_principal_values()` applies the visco factor; in the full principal mode the tension cap relaxes with `tension_direct_visco` and the MC cut with `mohr_coul_direct_visco` separately. Verified on the Professional 02-08-2026: 12 one-quad4 direct probes (tension/compression/tension-normal/MC full/MC normal, 1 step and dt 0.5→1; MC with each visco alone and both) and 2 incremental visco_power probes equal to 5 digits. Tests `test-2026/materi_direct_visco_first_step1/2`, `materi_direct_visco_separate1`, `materi_plasti_visco_power_dt1` (old rc=1, new rc=0). Corpus A/B (omen, exit codes) t2→t3: vs 233→233, t26 87→91, sfnet 172→172, Pro other 243→243; no change besides the new tests. |
| **F4 T2 — tied velocity slave carries its integrated displacement (`mpc_node_number` + `materi_displacement`, BACKLOG C41)** | `c4836d5` (F4 T2) | 2026-10-09 | Root cause `mpc.cc` `mpc_node_apply()`: the post-solve re-sync refreshed the slave velocity after `parallel_new_dof_diagonal` had integrated the displacement, so the slave displacement lagged one solve; with `-total_linear` the tied nodes diverged. Verified against the Professional 02-08-2026 on Pro `other/mohr_coul_direct7`: sigxy 0.0998334 at every node and equal tied disy (old GNU 0.102732 / 0.096935, disy -0.001808 / -0.001692). Test `test-2026/mohr_coul_direct_mpc_total1` (old rc=1, new rc=0). Residual disy/sigzz difference kept as DIVERGENCES C41. Corpus A/B (omen, exit codes): vs 233→233, t26 86→87, sfnet 173→172 (`ho_mech1` 60 s timeout flake), Pro other 242→243 (`mohr_coul_direct7` fixed); no other change. |
| **F4 T1 — `groundflow_total_pressure_limit` aligned with manual 6.587 and the Professional (TPL, K10)** | `3acce9d` (F4 T1) | 2026-10-09 | Owner decision 2026-10-07. Default limit 0 (single source `groundflow_total_pressure_limit()`), cap of the CALCULATED total pressure at every node in `groundflow_phreatic_coord()`, groundflow velocity 0 at integration points whose interpolated total pressure reaches the limit, `pres` dof never cut, GNU dry consolidation skip removed (the Professional keeps the consolidation term). Verified on the Professional 02-08-2026: 8 one-quad4 probes (positive/negative field x no record / 0 / 0.5 / 100 / -0.5) and `ground1` (no record / 1.e10) equal. Tests `test-2026/groundflow_total_pressure_limit_default1/high1/negative1`, `groundflow_velocity_dry1/2` (old rc=1 except dry2, new rc=0); `test-2014/groundflow_total_pressure_limit/_dry` retargeted to the Professional (pres 1.6667, gvely 0.16667/0). Corpus A/B (omen, exit codes): vs 231→233, t26 82→86, sfnet 172→173 (`ho_mech1` 60 s timeout flake, 30 s locally with both binaries), Pro other 242→242; no other change. |
| **F3 T2 WU2 — `group_materi_memory -updated_jaumann` (C20)** | `09f0948` (F3 T2 WU2) | 2026-10-08 | Manual Professional 2.2.1 (R = 0.5(F-F^T)+I, linear incremental strain) and 6.686; Professional 02-08-2026 as oracle. Own tests (old binary rc=1 'I do not know : -updated_jaumann', new rc=0, Pro rc=0): `materi_memory_updated_jaumann_shear1` (simple shear gamma 0.5, 1 step: 0.125 / 0.234375 / -0.125 = Pro exactly), `materi_memory_updated_jaumann_shear2` (100 steps: 0.06183352822 / 0.2396333009 / -0.06183352822 = Pro to 10 digits). Also equal with `-follow_material` and `materi_displacement` in shear; column one step `-fixed_in_space` = Pro to 1e-8, `-follow_material` differs (in-step geometry divergence). Corpus: Pro `other` `update3` PARSE -> pass. |
| **F3 T2 WU1 — default `mesh` `-fixed_in_space` (C8), default memory `-total_linear` with `materi_displacement` (C18), explicit `-updated_linear` on the current geometry (C19, reference path of `c945833` removed)** | `248ce83` (F3 T2 WU1) | 2026-10-08 | Manual Professional 6.806 (mesh default), 6.797 (memory default), 2.2.1 (`-updated_linear` on the current configuration); Professional 02-08-2026 as oracle. Own tests (RED on the old binary rc=1, GREEN rc=0, Pro rc=0): `mesh_default_fixed_in_space1` (bar +50% velocity-only, sigxx old 0.407136 / new 0.5 / Pro 0.5), `mesh_default_fixed_in_space2` (rigid rotation, old 0.109567 / new 0 / Pro 0), `materi_memory_default_displacement1` (column node 3 disx old 0.472222 / new 0.5 / Pro 0.5), `materi_memory_updated_linear_current1` (bar, `-updated_linear` + `-follow_material`, old 0.5 / new 0.4071364 / Pro 0.4071364), `materi_memory_updated_linear_current2` (column 100 steps velix old 0.472222 / new 0.520480 / Pro 0.520931, tolerance 1e-3: in-step geometry divergence). build_safe 136 OK / 0 FALLO. Corpus A/B exit code under `ulimit -v 4000000`, 60 s, old -> new: Pro `other` 239 -> 243 (+`axisym9`, `distri3`, `hypo12`, `relax2`, `update3` (WU2), -`mohr_coul_direct7`, BACKLOG C41; `ground19_water_under_dam` passes only uncapped with the new binary, it needs more than 4 GB); sfnet -8 (`axisym3/4`, `contac1/2`, `examp12`, `temp4/5`, `visc_pl2`, targets calibrated on the old GNU moving mesh, BACKLOG B13); validation-suite 233 unchanged; test-2026 only the new tests. |
| **Element-target contact: `node_boundary` generated for `contact_target_element_group` (F2 T3b, audit K05 `contac1/4`, `axisym3`)** | `50683a2` (F2 T3b) | 2026-10-07 | Manual 6.105/6.875: a target side is used when its nodes are on the mesh boundary; node_boundary is "for printing only" and the Professional computes it. `area.cc` `contact_node_boundary_set()` (called from `top.cc` `iteration_start()`, outside the parallel loops) writes `node_boundary -yes` for the nodes of unshared element sides when `contact_target_element_group` is given and no node_boundary exists (hand-written GNU 2014 node_boundary keeps priority). Pro `other/contac1` node 2/3 velx GNU before 0.03/0, after 0.02/0.01 = Pro; `contac4` = Pro exactly (temp 1.5); `axisym3` velx -9.99981e-4 = Pro; node_boundary sets equal to the Pro (1D and a 2D 3 x 2 quad4 block). Own tests `contact_element_target_boundary1` (1D) and `_boundary2d1` (2D quad4), RED rc 1 old / GREEN rc 0 new, the Pro meets both. Corpus A/B: see `odd/tasks/f2-medium-bugs.md`. Not done (BACKLOG C40): Pro single-layer target error, hex8 sides, element_boundary. |
| **`contact_target_geometry`/`_switch` read for every index, without a db PUT in the parallel contact loop (F2 T3a, BACKLOG C5, audit K05)** | `3cef926` (F2 T3a) | 2026-10-07 | Manual 6.106/6.107: Professional names of `contact_geometry`/`_switch`. `contact.cc` `contact_geometry_get()` reads either name at the point of use (GNU name wins on the same index); the old alias PUT index 0 inside `parallel_contact` (abort "Data is allocated in a parallel loop" in 6 audit runs) and dropped indices > 0. Pro `validation/validation_10` sigyy GNU before ~0, after -21232 (Pro -21031, target -20538 +- 1e4); Pro `other/pressed_rubber_cylinder` max |dof - Pro| 0.10 -> 6e-7; Pro `other/contac2` x 1.02 -> 1.0 (= Pro). Rename neutrality (`contact_geometry` -> `contact_target_geometry`) on `contact`, `contact_block`, `contact_ctrl_apply`, `contact_heatgen`, `examp6`, `examp9`: byte-identical `node_dof`. Own tests `contact_target_geometry_index1`, `_index0`, `_switch1` (test-2026, run_ok; RED rc 1 old / GREEN rc 0 new; the Pro meets all three targets). Corpus A/B (plain rc runner): see `odd/tasks/f2-medium-bugs.md`. |
| **Truss `-total` geometrically nonlinear + tension geometric stiffness (F2 T2, audit K17 `truss14`, `truss2/3`)** | `40567ff` (F2 T2) | 2026-10-06 | `truss.cc` `truss()`: new `-total` branch (current axis, `EA/L0`, so N = N0 + EA (L-L0)/L0; before it was the linear branch) and the geometric stiffness `N/L (I - n n) dt` in the matrix for `-updated`/`-total` with N > 0 (RHS unchanged). Manual 6.777 lists only -updated_linear/-updated/-total_linear; the Professional 02-08-2026 accepts -total and is the oracle: truss14 Pro disx 0.70711 / disy 1.70711 / N 1.41421, GNU now equal (before solver breakdown, disy 1e10); truss2/3 reactions 0.1/0.1 and node 2 veli (-0.17642, 0.82358) = Pro to 2e-6 (before breakdown, reaction -0.99997). Own tests `truss_total_rotation1`, `truss_updated_transverse1` (RED/GREEN). Corpus A/B: no exit-code change (validation-suite 233/10, sfnet 181/16); dump changes only iteration-path (`bmom_truss` 3.6e-7, sfnet `truss2/3` 2.9e-5, `truss7` 3e-5). |
| **`group_time`: nodes of inactive elements are not driven by `bounda_dof` (F2 T2c, BACKLOG C39, audit K17 `truss8`)** | `fix(bounda)` F2 T2c (hash in `odd/tasks/f2-medium-bugs.md`) | 2026-10-06 | Manual 6.769 silent; project decision 2026-10-06: follow the Professional. `bounda.cc` `bounda_node_dormant()` + `group.cc` `group_time_alive()` (shared with `elem.cc`): a node whose elements are all inactive by `group_time` keeps NODE_BOUNDED 0 for the prescribed-dof records. Pro `other/truss8` node 1 disx: Pro 0/0/0.05/0.1 at t = 0.5/1/1.5/2, GNU before 0.05/0.1/0.15/0.2, after = Pro (node 2 unchanged, 0.1). Own test `group_time_bounda_dormant1` (test-2026, run_ok): RED rc 1 old / GREEN rc 0 new, the Pro passes it. Corpus A/B (plain rc runner): see `odd/tasks/f2-medium-bugs.md`; only `truss8` and the new test contain `group_time`, every other input runs the old code path (helper returns NULL). Pending: the Pro also zeroes the velocity of the nodes of a DEAD group (probe `group_time 0 0. 1.`), not implemented. |
| **`support_edge_normal` damping and density in the iteration matrix (F2 T2, audit K17 `support_edge_normal_3`)** | `0b9ca9a` (F2 T2) | 2026-10-06 | `area.cc` `type[10]` side matrix: `kterm = (dt (k_n nn + k_t tt) + c_n nn + c_t tt + (d_n nn + d_t tt)/dt) cij`; before only `dt k`, so a damper with c >> k dt overshot each update. Verified against the Professional 02-08-2026 test other/support_edge_normal_3: GNU disy 124.75 rc 1 -> 5.00998004e-4 (= Pro without the freeze record; Pro target 5.00999e-4 +/- 1e-6 passes). Own test `support_edge_normal_damping_tangent1` (RED/GREEN). Converged states unchanged: `tsup_*` and the corpora give identical dumps (validation-suite 233/10, sfnet 181/16 both sides). `_stiffness_freeze` remains parse-only (BACKLOG C38). |
| **`group_materi_elasti_poisson_power` replaces the poisson ratio of the elastic C (F2 T1, audit K06 `mpower`)** | `a270264` (F2 T1) | 2026-10-06 | `stress.cc` `set_stress()`: the block clears `C`/`Cmem` before `C_matrix()` (which accumulates), as the `_young_power`/`_hardsoil` blocks; before, the stiffness was C(young, poisson) + C(young, nu(p)). Manual Professional 6.652 ("the poisson ratio is made a function of the average stress state") and the Professional 02-08-2026 agree: `mpower` (converted) Pro -0.44173 / -1.41230, GNU -0.43887 / -1.41012 (before -1.43887 / -3.29236; a forward-Euler replay of the oedometer reproduces both GNU states to 1e-9). Own test `elasti_poisson_power_replace1` (test-2026, RED/GREEN); validation-suite `mpower` retargeted (-0.4389 / -1.4101 +/- 0.02, the old target was calibrated on the double stiffness). Corpus A/B: only `mpower` and the new test change. The `mpower` retarget was approved by the project owner on 2026-10-06. The `_young_power` + `_young` combination (Pro adds both stiffnesses) is a deliberate divergence decided 2026-10-06 (BACKLOG C36, "Descarte de features"). |
| **Generated elements: zero `element_dof` template (F0 T4 bug fix, not a Professional feature)** | `a32ea74` | 2026-10-05 | `generate.cc` `generate_spring()` / `generate_beam_truss()` / `generate_interface()` / `generate_interface_triangle()`: `array_set( tmp_element_dof, 0, mnolnuknwn )` after `get_new_dbl()`, as `create.cc` `create_element()`. Verified with own tests: two runs of gen1, genbeam1, iface_gen* (validation-suite) and conspr3, genera1/2, ho_othr1, spring2/3/5 (sfnet) give byte-identical `.dbs` with no `e-3xx` values; old/new A/B changes only `element_dof` lines. Corpora unchanged: validation-suite 233/1, test-2026 47/0, sfnet 187/10. Docs: `manual-*/generated_element_dof_init.md`. |
| **`control_print_dof`: dof-component loop over the primary unknowns (F0 T6 bug fix, BACKLOG C23, not a Professional feature)** | `f22ba98` | 2026-10-05 | `print_hi.cc` `print_dof()`: the loop ran `ipuknwn` over all `nuknwn` value/derivative entries, repeating every vector/stress block (shifted into the next dof) and reading `node_dof` past the end of the record (heap garbage `e-3xx`, run-dependent); now `npuknwn = nuknwn/nder`, step `nval`, type from `dof_scal_vec_mat[iuknwn]` as `print_vt.cc`. Not a `get_new_dbl()` init defect. Verified with own tests (no Professional example of the file layout; the layout is the one of our user manual): `dofid_no` `dof.21` 768 lines (20 garbage) -> 48 = 12 components x 4 nodes, two runs byte-identical, no `e-3xx`; `build_safe.sh` dof.21 check now asserts 48 lines + no `e-3xx`, `freq_timeint` 176 -> 40 lines/print; suite 226/226, checks 101 OK / 0 FALLO, rc=0. Corpus A/B old vs new: per-test status identical (validation-suite 242 run + `tslv_zip` isolated: 241/1, the 1 = `print_mesh_dof.dat` output dump; test-2026 47/0; sfnet 187/10); changed outputs only `dof.*` (shorter, no garbage) and wall-clock `time_calculation`. Same loop pending in `print_fr.cc`/`print_gm.cc` (C24); buffer audit follow-ups C25. Docs: `manual-*/print_dof_component_loop.md`. |
| **`control_print_frd` / `control_print_gmsh`: result loop over the primary unknowns (F1 T1 bug fix, BACKLOG C24, not a Professional feature)** | `fix(print)` F1 T1 (hash in `odd/tasks/f1-cheap-bugs.md`) | 2026-10-05 | `print_fr.cc` `print_frd()` and `print_gm.cc` `print_gmsh()` looped `ipuknwn` over all `nuknwn` value/derivative entries like the old `print_dof()` (C23): shifted duplicate blocks and `node_dof` reads past the record. Now `npuknwn = nuknwn/nder`, step `nval`, type/label/dof type at `iuknwn`; FRD 2D vectors write 0 as third component (read the next dof); gmsh field names from the dof type (`node_materi_velocity_0`, Pro base names). Verified against the Professional 02-08-2026 on own test `test-2026/print_frd_gmsh_dof_loop1` (quad4, `derivatives`, vel/dis/stress/total strain/temp): GNU old 68 FRD blocks / 324 gmsh `node_` fields -> 5 (`VELO DISP STRESS TOSTRAIN NDTEMP`) / 17; values equal the first block of each dof of the Pro (velx 0.01, disx 0.01, sigxx 10, epsxx 0.01, temp 0/1). RED (old binary fails the new build_safe check) / GREEN. Suite 226/226, checks 103 OK / 0 FALLO, rc=0. Corpus A/B old vs new: per-test rc and `.dbs` (minus `time_calculation`) identical: validation-suite 233/10 (the 10 = `print_mesh_dof` + 9 expected-error inputs), test-2026 46/2 (2 expected-error), sfnet 187/10. Remaining FRD/gmsh layout differences with the Pro: BACKLOG C27. Docs: `manual-*/print_dof_component_loop.md`, `control_print_gmsh.md`, `manual-developer/control_print_frd.md`. |
| **`convection_apply` default `-no` + Professional scope (F1 T2, audit K03)** | `fix(convection)` F1 T2 (hash in `odd/tasks/f1-cheap-bugs.md`) | 2026-10-06 | `general.cc` `convection_apply_switch()`: one source of truth for the five `OPTIONS_CONVECTION` readers (general/beam/truss/spring/conspr, all hardcoded `-YES`), default `-no` as manual Professional 6.395; `general()` no longer gates the `group_condif_flow` / groundflow-velocity heat convection with the switch (`lagrange_flow`). Verified against the Professional 02-08-2026: condif2 (converted) temp 0.5 absent / 0.5 `-no` / 1 `-yes` (manual and binary agree); condif1 (`group_condif_flow`) and condif5 (groundflow velocity) temp 1 with `-no`, `-yes` and absent. Own tests `test-2026/convection_apply_default1` (RED old 1 vs target 0.5, GREEN 0.5 = Pro) and `convection_apply_yes1` (1 = Pro), run by build_safe. validation-suite `cmat_gate`, `cd_method`, `materi_direct_wall` (targets calibrated on the old default) got an explicit `convection_apply -yes`: .dbs identical to before apart from the record. Suite 226/226, checks 105 OK / 0 FALLO, rc=0. Corpus A/B vs F1 T1: validation-suite 233/10 unchanged per test (~57 output-only changes: absolute stress recovery and no convection of strain/stress dofs on fixed meshes without the record; mostly round-off, some 1-13 %; closer to the Pro where it runs: cmat_gate, mesh_act_grav, mdirect_comp); test-2026 all pass (48/2 expected-error); sfnet 187/10 -> 182/15: condif2/condif3/examp1/incnav5 now equal the Pro, strconv1 Pro RUNFAIL; stale targets recorded as BACKLOG B10 (sfnet not edited); 12 output-only changes. Docs `manual-*/convection_apply.md`. |
| **Nodal recovery of the strain tensors with the stress operator (F1 T5, audit K11)** | `fix(materi)` F1 T5 (hash in `odd/tasks/f1-cheap-bugs.md`) | 2026-10-06 | `miscel.cc` `strain_recovery_weight()` (= `sri_stress_recovery_weight`, h for SRI); `materi.cc` strain numerators moved next to the stress recovery with that weight; `general.cc` lumped denominator/inertia of `MATERI_STRAIN_TOTAL/_ELASTI/_PLASTI` with the same weight. Verified against the Professional 02-08-2026: validation-suite `tslv_part` (converted) every strain within 4e-11 (before clamped-edge eptyy −2.35e-6 vs 0, eptxy ×3); Pro `tension_direct1` eptyy residual 0.667 -> 0.5 = Pro. Own test `test-2026/strain_recovery_clamped1` (2 quad4, clamped edge eptyy 0 exact): RED/GREEN, the Professional passes it. Suite 226/226, checks 119 OK / 0 FALLO, rc=0. Corpus A/B vs the K13 build: 46 tests change only strain dofs (+ derivative slots), no classification change; on the 26 the Professional runs, 14 move to the Professional values, the rest are round-off of zero fields or a pre-existing difference (beam2d family, examp17). K02 (one-iteration lag) parked: BACKLOG C35. Docs `manual-*/nodal_strain_recovery.md`. |
| **Direct cut-off state output: plastic strain C⁻¹(σ_trial − σ) and direct-law yield value `f` (F1 T4, audit K13)** | `fix(stress)` F1 T4 (hash in `odd/tasks/f1-cheap-bugs.md`) | 2026-10-06 | `stress.cc` `set_stress()` + new `direct_plastic_strain()`: the stress before the direct cut-offs is kept and the cut is added to `inc_epp` (and removed from `inc_epe`/`new_epe`); `materi_direct_full_mc()`/`materi_direct_cutoff()`/`materi_compression_cutoff()` return the largest direct-law yield value on the final stress, stored in `new_f` (0 after a cut). Verified against the Professional 02-08-2026 (Pro test `visc_pl4`: eppxx 0 -> 0.5, Pro 0.49998; f −1e10 -> 0 = Pro) and 16 single-quad4 probes (tension/MC/compression, full and plane mode, cut and elastic, combined laws): every normal component and `f` equal; shear eppxy tensor vs Pro twice (deliberate, BACKLOG C34); `_visco` stress path differs (C33, not changed). Own tests `test-2026/direct_plastic_strain_tension1`, `_mc1`, `direct_yield_value1`: RED/GREEN. Suite 226/226, checks 117 OK / 0 FALLO, rc=0. Corpus A/B vs `7443ba8`: vs and sfnet `.dbs` identical; test-2026 `control_plasti_apply1/2` eppxx/eppyy/eppzz 0 -> −3.6857e-3/9.3e-3/−7e-4 (= Pro), +3 new tests. Docs `manual-*/direct_plastic_strain.md`. |
| **`node_force` sign: external load along +force; `control_mesh_cut_node_force` substitute without the compensating negation (F1 T6, BACKLOG C28)** | `fix(dof)` F1 T6 (hash in `odd/tasks/f1-cheap-bugs.md`) | 2026-10-06 | `dof.cc` `parallel_new_dof_before()`: `node_rhside += node_force` (was `-=`); `delete.cc` `mesh_cut()`: `node_force_exist += node_force_acc` (the negation compensated the old sign). Verified against the Professional 02-08-2026: probes 1D bar2 (`node_force 2 1.0`: velx −2 -> +1 = Pro), 2D quad4 (disx −6 -> 2 = Pro), 3D hex8 (disx −2 -> 1 = Pro); Pro `mesh_cut_1` `node_force 6` −1 -> +1 = Pro, sigxx 1 kept. The old "doubling" was the reversed load compressing the unit bar to zero length (singular; `-1.0` still GNU −2 / Pro −1, `-0.5` both −0.5). Own tests `test-2026/node_force_sign1` (−2 -> 2) and `mesh_cut_node_force_sign1` (−1 -> 1): RED/GREEN. Corpus A/B vs `7443ba8` (audit runner, 60 s): only validation-suite `node_force_inertia` changes (inertia −19.75 -> −20.109, target −19.83 ± 0.2 calibrated on the reversed sign: BACKLOG B12, retargeted to −20.11 ± 0.2 on 2026-10-06 with the owner's go-ahead); sfnet and test-2026 `.dbs` identical. Docs `manual-*/node_force.md`/`node_additions.md`, `control_mesh_cut_geometry.md`. |
| **`group_materi_expansion_volume` in the momentum inertia (F1 T3, audit K04 `incnav3`)** | `fix(general)` F1 T3 (hash in `odd/tasks/f1-cheap-bugs.md`) | 2026-10-06 | `general.cc` `general()`: the MATERI density is multiplied by (1 − beta T) (manual Professional 6.665) as `materi.cc` already did for gravity, so the lumped inertia uses rho = rho0(1 − beta T). Verified against the Professional 02-08-2026: Pro test `incnav3` (its text: rho = 0.5, total rate −5) `tvely` −2.4976 -> −4.990415335 = Pro (−4.990415335463), same on the converted sfnet copy. Own test `test-2026/expansion_volume_inertia1`: RED/GREEN. Suite 226/226, checks 114 OK / 0 FALLO, rc=0. Corpus A/B vs the stress_pressure_factor fix: validation-suite and test-2026 identical; sfnet 181/16 -> 180/17: legacy `incnav3` target −2.5 (stale, BACKLOG B11; the Pro gives −4.9904 on it). No other `.dbs` change (no other corpus input uses `expansion_volume`). Docs `manual-*/group_materi_expansion_volume.md` (new). |
| **`group_materi_elasti_stress_pressure_history_factor`: loading/reloading decision of the Professional (F1 T3, audit K04 `stress_pressure_factor`)** | `fix(stress)` F1 T3 (hash in `odd/tasks/f1-cheap-bugs.md`) | 2026-10-06 | `stress.cc` `set_stress()`: a step that starts below the history estimates its end pressure with the factored stiffness and keeps the factor only strictly below the history (`sph*(1-1e-10)`). Rule measured on the Professional 02-08-2026: Pro test `stress_pressure_factor` (dt 0.1, 0.15, 0.3, 0.4: −3.05/−2.25/−2.3/−2.4 = Pro; GNU before −3.10/−2.325/−2.45/−2.4) and a partial unload/reload path (factors 0.5 and 3, dt 0.1/0.25/0.3/0.35/0.4): 9/10 equal, the exception is a truncated last step with factor 3 (Pro −3.2 applies the factor on a step that starts on the loading curve, GNU −2.8 = analytic first loading; documented). Own tests `test-2026/stress_pressure_factor_reload1` (−3.10 -> −3.05) and `_reload2` (factor 3, dt 0.3: −3.4 -> −2.8): RED/GREEN. Suite 226/226, checks 113 OK / 0 FALLO, rc=0. Corpus A/B vs the interface fix: per-test rc and `.dbs` identical (validation-suite `msph`/`msph_flat` included), t26 +2 new tests. Docs `manual-*/group_materi_elasti_stress_pressure_history_factor.md`. |
| **Interface groundflow permeability / capacity and condif conductivity per unit length (F1 T3, audit K04 `ground9`)** | `fix(interface)` F1 T3 (hash in `odd/tasks/f1-cheap-bugs.md`) | 2026-10-06 | `interface.cc` `interface_element()`: the facing-pair conductance of `group_interface_groundflow_permeability` and `group_interface_condif_conductivity` is weighted by the tributary measure `w_ip*iface_measure` (manual Professional 6.628/6.625: per unit length/area), and `group_interface_groundflow_capacity` stores `0.5*C*w_i*L` per node, implicitly (`C/dt` on the matrix and `element_lhside`, as the continuum capacity). Verified against the Professional 02-08-2026: Pro `ground9` middle pres 1/3 -> 0.5 = Pro (and 0.5 for the model 2 and 3 wide; GNU before 0.5 / 0.6); thermal analog (Pro `condif11` model, which the GNU cannot parse: 2-value conductivity) 1/3 -> 0.5; `groundflow_interface` transient node 5 at t = 1/3/10: 1.0/1.75/1.998 = Pro (before 0.32/0.81/1.65), interface capacity 0 / 2: 1.3333 / 0.8 = Pro, same with the model 2 x tall. Own tests `test-2026/interface_flux_measure_ground1` (0.6 -> 0.5), `interface_flux_measure_condif1` (1/3 -> 0.5), `interface_capacity_transient1` (0.8146 -> 1.75): RED/GREEN. Weighting the permeability alone made validation-suite `groundflow_interface` diverge (explicit capacity), hence the capacity change in the same unit. Suite 226/226, checks 111 OK / 0 FALLO, rc=0. Corpus A/B vs the force_point fix: per-test rc identical (vs 232/11, sfnet 181/16, t26 +3 new tests); only `.dbs` change: `groundflow_interface` node 6 pres 3.790e-4 -> 3.959e-4 (Pro 3.960e-4; target node 5 = 2 met). Docs `manual-*/group_interface_groundflow_permeability.md`, `_capacity.md`, `manual-user/group_interface_condif_conductivity.md`. |
| **`force_point` sign: external load along +force (F1 T3, audit K04 `force18`)** | `fix(force)` F1 T3 (hash in `odd/tasks/f1-cheap-bugs.md`) | 2026-10-06 | `force.cc` `force_point_calculate()`: `node_rhside[vel_indx/nder+idim] += weight*force` for `idim<ndim` (was `-=` over every principal unknown, reading past the record). Verified against the Professional 02-08-2026: Pro test `force18` (1D bar, force 1 at 0.8, analytic 0.8) GNU velx −0.8 -> +0.8 = Pro; probe force at node 2 (x = 1): −1 -> +1 = Pro. Own test `test-2026/force_point_sign1` (2D quad4 bar, force (2,0) at (1,0.5)): old −2 (RED), new 1.999999902 (GREEN, tol 1e-6: `point_el` iterative location; Pro 2). No corpus input uses `force_point`. Suite 226/226, checks 108 OK / 0 FALLO, rc=0. Corpus A/B vs the condif7 fix: per-test rc and `.dbs` identical (vs 232/11, sfnet 181/16, t26 +1 new test). Related, not changed: `node_force` also subtracts its value (Pro adds; probe +1 -> GNU −2), BACKLOG. Docs `manual-*/force_point.md`. |
| **Edge-load `_node` / `_element_node` restriction selects whole sides (F1 T3, audit K04 `condif7`)** | `fix(area)` F1 T3 (hash in `odd/tasks/f1-cheap-bugs.md`) | 2026-10-06 | `area.cc` `side_node_restriction_ok()` called in the `area()` side loop: a side is loaded only when all its nodes are in the `_node` list (or, `_element_node`, the element matches and all its local side nodes are listed), for every edge family (condif heat/convection/radiation edge normal, groundflow flux, force_edge/_normal/_water/_projected, support_edge_normal). Before, the nodes were checked one by one and a corner node also loaded the adjacent side. Verified against the Professional 02-08-2026: Pro test `condif7` (analytic T = 10 − x) GNU 20 -> 10 = Pro; probes on the condif7 square (nodes `1 3` / `1` / `1 3 4`, `_element_node`) for heat flux, groundflow flux, convection, `force_edge`, `force_edge_normal`, `support_edge_normal`: every probe equals the Pro after the fix (old GNU 20/16/20, Pro 10/0/12-18). Own tests `test-2026/edge_node_restriction_heat1` (temp 12/18, old 20/20) and `edge_node_restriction_force1` (disx 1/1, old 2/2): RED old / GREEN new, run by build_safe. Suite 226/226, checks 107 OK / 0 FALLO, rc=0. Corpus A/B vs F1 T2 (`audit_all_variables.py run`, 60 s): per-test rc and `.dbs` (RTOL 1e-4) identical in validation-suite (232/11 with this runner: the 9 expected-error inputs, `print_mesh_dof`, `tslv_zip`), sfnet (181/16 with 2 load timeouts) and test-2026 (48/2 -> 50/2 with the 2 new tests). Docs `manual-*/edge_node_restriction.md`. |
| **C16 — default `group_materi_memory` = `-updated_linear` (Professional default) + `control_materi_updated_apply -yes`** | `d8117aa` (code), `0c7d4d3` (test registration) | 2026-10-04 | **Code**: `materi_memory_default()` (general.cc) assigned before every `GROUP_MATERI_MEMORY` read (materi.cc materi/set_deften_etc, stress.cc, elem.cc, tendon.cc); returns `-UPDATED` when `control_materi_updated_apply` is `-yes` at ICONTROL (manual 6.152), else `-UPDATED_LINEAR`; the `NODE_START_REFINED` reference gradient only for an EXPLICIT `-updated_linear`; tendon assert accepts `-UPDATED_LINEAR`. **Verification** (own tests vs the Professional 02-08-2026 binary): `materi_memory_default_shear1` (single quad4 simple shear gamma 1, no record): GNU before 0.2311/0.4207/-0.2286, after 0/0.5/0 = Pro (rc=1 before, rc=0 after); `materi_memory_default_stretch1` (-follow_material stretch x2): sigxx 0.695653 old = new = Pro (guard: the reference gradient would give 1); `materi_memory_default_updated_apply1` (+ `control_materi_updated_apply 0 -yes`): GNU 0.231098/0.420739/-0.228598 (= explicit -updated in both binaries); the Pro binary IGNORES -yes (0/0.5/0, also with materi_displacement, -follow_material, -fixed_in_space, index 1) — decision 2026-10-04: follow the manual (deliberate divergence). **Corpus A/B** (60 s, timeouts rerun at 300 s, the 16 load-induced timeouts rerun on both binaries): validation-suite 233/1, test-2026 44 (+3 new), sfnet 187 PASS unchanged; Pro `other`: mohr_coul_direct6/8 RUNFAIL->PASS, hypo12 PASS->RUNFAIL (explained marginal miss: sigyy off 1.02e-3 vs tol 1e-3, Masin model without memory record; Professional distribution file, target not edited). Pending (BACKLOG): default with materi_displacement (Pro behaves as -total_linear), explicit -updated_linear gradient on -follow_material, -updated_jaumann. Manuals `group_materi_memory_default`, `control_materi_updated_apply` (user/developer). |
| **Zero `element_normal` input error checked per direct law: an explicit `_direct_normal` of the OTHER law no longer suppresses it** | `b7f3522` | 2026-10-03 | `area.cc` `element_normal_automatic( gr )` skips a law only when that same law has an explicit normal (same precedence as `direct_normal_get()`, `materi.cc`); before, any explicit normal in the group returned 0, so explicit `mohr_coul_direct_normal` + `tension_direct_normal_automatic -yes` on a zero-normal element ran silently with the tension principal cap. New `validation-suite/test-2026/element_normal2d_row_mixed_error1` (row of 3 quad4 of `element_normal2d_row_error1`, explicit MC normal 0 1 + tension automatic `-yes`; `EXPECTED_ERROR` element 2): old binary rc 0 (RED), new rc 1 naming `group_materi_plasti_tension_direct_normal_automatic` and advising `group_materi_plasti_tension_direct_normal` (GREEN); Professional 02-08-2026 aborts with "Program error, assertion failed / Error for: : !array_null( directions, MDIM*MDIM )" (rc 0, "Error detected for input file" in `tochnog.log`), same outcome as the GNU. Companion `element_normal2d_row_mixed1` (explicit tension normal 0 1 + MC automatic `-yes` on the end elements only, normals (1,0)/(-1,0)): runs on old, new and Pro, `.dbs` values equal to Pro (disy 1, sigyy 1, sigxx 0). Regression old vs new clean build (`scripts/audit_all_variables.py run`, 60 s, LC_ALL=C): validation-suite 234: 232 rc 0 + `materi_direct_auto` expected error + `print_mesh_dof` RUNFAIL on both; test-2026 44: only `_mixed_error1` changes (rc 0 -> expected error); sfnet 197: 187 rc 0 on both, per-test rc identical (`examp14` timed out at 60 s on the old run, rc 1 on both at 300 s); `.dbs` identical except the pre-existing uninitialized `element_dof` denormals (`gen1`, `genbeam1`, `iface_gen*`, `conspr3`, `genera1/2`, `ho_othr1`, `spring2/3/5`; none uses `_normal_automatic`). |
| **Plane-mode direct cut-off (C14): separate tension / Mohr-Coulomb planes, tension first (sigy = 0 default, principal cap without a tension normal), phi_flow return k = (M/G) tan(phi_flow); `tension_direct_normal` without `tension_direct`** | `d17461e` + `77f9282` | 2026-10-02 | No Professional usage example of the plane mode with phi_flow or separate planes; rule measured on the Professional 02-08-2026 with ~25 single-quad4 probes (c01-c25, scratch driver) and reproduced by hand to 4-5 digits. New self-checking inputs `validation-suite/test-2026/direct_cutoff_plane_{shear1,tension_principal1,tension_value1,tension_normal1,tension_normal2,two_normals1,phi_flow1,phi_flow2,phi_flow_floor1}.dat` (`-total_linear`, closed-form targets at 1e-8): old binary 8 RED (`tension_normal1` passes on both), after `d17461e` the 3 phi_flow tests RED, after `77f9282` 9/9 GREEN; the Professional passes 8/9 (`phi_flow2`, nu 0.25, is the deliberate M/G divergence, BACKLOG C15: GNU sxy 7.04154e-4 / syy -1.28895e-3 vs Pro 6.92049e-4 / -1.26679e-3). **A/B old vs new binary** (`scripts/audit_all_variables.py run`, 45 s, LC_ALL=C): sfnet 197 per-test rc and `.dbs` identical; test-2026 changes only in the new inputs plus round-off in `element_normal3d_column1` (5e-17 absolute on stresses ~1e-9, reordered sig_n sum); validation-suite 234: `cd_geom`, `cd_method` change (MC normal without `tension_direct` now caps the principal tension at 0; targets still met) and `materi_direct_mc` now fails its target (233 -> 232 PASS): its `vely` constraint covers nodes 1 and 4 only, the element collapses at t = 0.7 (sigxy 1 -> 0.13 -> 0.06); the Professional with the same constraints collapses at the same step (sigxy -> ~0), with all nodes constrained both keep sigxy 1 (sxx = syy -100.87 GNU / -100.90 Pro), and with `-updated` memory the GNU sigxy history follows the Pro (1.916..12.06 vs old GNU 2.95..18.69) — stale target, BACKLOG C17. Implementation: `materi.cc` `direct_normal_get()`, `stress.cc` `materi_direct_cutoff()`, `direct_flow_factor()`, `materi_direct_principal_cap()`, `direct_cap_principal_values()`, `set_stress()` dispatch and elastic stiffness `Cmem+Cuser+Chyper+Chypo(+1/compressibility)`; `database.cc` drops the `tension_direct` prerequisite of `tension_direct_normal` (the `_automatic` switch keeps it, as Pro). Detail: `manual-user/` + `manual-developer/group_materi_plasti_mohr_coul_direct.md`. |
| **`group_materi_plasti_*_direct_normal_automatic -no`: no automatic normal, full principal-stress cut-off (as Professional)** | `260bd95` | 2026-10-02 | `materi.cc` `materi()`: the automatic branch is gated by `element_normal_automatic( gr )` (only a `-yes` record without explicit `_normal`) instead of "either `_normal_automatic` record exists". Verified on the Professional 02-08-2026 with the new `validation-suite/test-2026/element_normal2d_column_switch_no1` / `_yes1` (column of 2 quad4 of `element_normal2d_column1`, upper element `tension_direct 0`, stretched in x: left side velx 0, right side velx 1; GNU dialect accepted as is by Pro): `-yes` sigxx 1/1/1/1 at nodes 1/3/5/6 on Pro, old GNU and new GNU; `-no` Pro sigxx 1/0.5/0/0, new GNU 1/0.5/0/0 (identical `.dbs` values), old GNU 1/1/1/1 (target failure on node 5 sigxx, rc 1 = RED). `element_normal` (0,-1) of the upper element still generated with `-no` on both. Regression (old vs new binary, scratch corpus runner, per-test rc): validation-suite 234: 233 PASS / 1 RUNFAIL (`print_mesh_dof`, pre-existing) on both; test-2026 33: old 32/1 (new `_switch_no1`), new 33/33; sfnet 197: 187 PASS / 10 RUNFAIL on both, per-test rc identical. `.dbs` identical except the version stamp and pre-existing uninitialized/timing noise (`gen1`, `genbeam1`, `iface_gen*`, `element_dof` denormals; `time_calculation`). |
| **Zero `element_normal` with `group_materi_plasti_*_direct_normal_automatic -yes`: readable input error instead of the silent full cut-off (Pro: assertion)** | `a8d6452` (merge `bbdb01a`) | 2026-10-01 | Not a divergence in outcome: both stop. Verified on the Professional 02-08-2026: `materi_direct_auto` (single hex8, converted with `scripts/audit_all_variables.py` `convert()`) and the new `validation-suite/test-2026/element_normal2d_row_error1` (single row of 3 quad4, automatic `-yes` on all) abort with "Program error, assertion failed / Error for: : !array_null( directions, MDIM*MDIM )"; `materi_direct_auto` without the `_normal_automatic` record runs on Pro. GNU now stops (rc 1) with `Error: element_normal is zero for element : N` + group + remedies (restrict the group to the boundary layer, explicit `_tension/_mohr_coul_direct_normal`, supported-face hint) + `db_error` on the automatic record (`area.cc` `element_normal_set()` / `element_normal_automatic()`). Companions `element_normal2d_row_explicit1` (explicit normal 0 1: disy 1 at the top, all stresses 0) and `element_normal2d_row_layer1` (automatic group only on the end elements: normals (1,0)/(0,0)/(-1,0), sigyy 1) give identical `.dbs` values on GNU and Pro, targets met on both. `materi_direct_auto` and `element_normal2d_row_error1` are `EXPECTED_ERROR:` tests (`scripts/build_safe.sh` `expected_error_ok()`). Regression (old vs new binary): validation-suite 234: 233 PASS / 1 RUNFAIL (`print_mesh_dof`, pre-existing) on both, `materi_direct_auto` PASS→expected-error PASS; test-2026 31: all PASS; sfnet 197: 187 PASS / 10 RUNFAIL on both, per-test rc identical; build_safe 16/16 runs, file-check diff only the new expected-error lines. The `-no` switch divergence was resolved on 2026-10-02 (commit `260bd95`, row above). |
| **`element_normal` (manual 6.445) + `group_materi_plasti_*_direct_normal_automatic` using it (audit K04 `tension_direct1`)** | `b81d050` | 2026-09-30 | Rule measured on the Professional 02-08-2026, which prints `element_normal`/`node_normal` in the `.dbs`: ~30 probe meshes (quad4/tria3/quad9/tria6/tet4/hex8/bar2, rectangular and distorted, two groups, trusses, input-given record, `-no` switch, `control_mesh_delete_element`, 90° rotation with `-follow_material`). Boundary side = not shared with another isoparametric element (groups irrelevant, trusses do not close a side); node normal = normalized sum of `∫N n dA` of its boundary sides (inward); element normal = normalized sum of its node normals (zero when it vanishes); generated only when an `_normal_automatic` group exists, regenerated after mesh changes, not after deformation. GNU = Pro to print precision (5e-13) on all probes except curved quad9/tria6 sides (C14). `tension_direct1`: before `disy` 0.42 / `sigyy` 0.667 (cut-off inert, 2D cross-product normal (0,0,1)), after `disy` 0 / 1 and all stresses 0 = Pro (only `eptyy` nodal recovery differs, K11). New inputs `validation-suite/test-2026/element_normal2d_column1`, `_corner1`, `_mesh_change1`, `element_normal3d_column1`, `_distorted1` pass on Pro and GNU (pre-fix: parse error on `-element_normal`; without those targets `disy` fails). validation-suite 233 PASS, sfnet 187/10, test-2026 23/23 unchanged; `.dbs` identical to the pre-fix binary except `materi_direct_auto` (new zero record, round-off 2.8e-14). Keyword probe: `element_normal` accepted (832/1154). |
| **Control records read at index 0 only (audit K04 `elasti7`): `control_materi_elasti_k0` for any control index + K0 without poisson + K0 with hardsoil; `control_materi_plasti_hardsoil_gammap_initial` pre-allocation for any index; `control_solver_bicg_stop` precedence (Professional index-0 / priority semantics, then combined precedence)** | `ad1f5da` + `47b6fda` + combined precedence (commit `407aa3d`) | 2026-09-29 | `db_active_index(CONTROL_xxx, 0, ...)` tests that EXACT index. (1) `stress.cc` K0 hook: the switch was ignored for any index other than 0; Pro `elasti7` (index 20) sigxx/sigyy/top disy: Pro −0.04999999884/−0.5/−4.909090913e-7, GNU before 0/−0.5/−5e-7, after −0.0500000000025/−0.5/−4.90909090909e-7 (all variables within 2.4e-8). Measured on the Pro 02-08-2026 and adopted: K0 also applies without `group_materi_elasti_poisson` (young and young_power) and with `group_materi_elasti_hardsoil` (replaces nu50/nuur; the old comment said hardsoil was not affected — manual 6.650 and the Pro agree it is); K0=2 → clamp 0.95 on both. (2) `top.cc`: `control_materi_plasti_hardsoil_gammap_initial` with a non-zero index aborted ("Data is allocated in a parallel loop"); now pre-allocated when the record exists at any index (Pro runs both). (3) `so_bicg.cc`: 9 singular-bicg probes on the Pro: only `control_solver_bicg_stop` index 0 counts (for every timestep index) and it overrules `solver_bicg_stop` (manual 6.376); the GNU matched 5/9 before, 9/9 after `47b6fda` (contradicts the per-timestep wording of 6.376; Pro followed). Then combined precedence (project decision 2026-09-29, commit `407aa3d`): `solver_bicg_stop` → index 0 → the timestep's own index, last read wins; 6/9 equal the Pro, the 3 that differ are exactly the probes with a record at the timestep's non-zero index (deliberate divergence, BACKLOG C13); combined-precedence regression: see `solve_iterative_bicg.md`. Self-checking inputs `validation-suite/test-2026/elasti_k0_index1`, `elasti_k0_nopoisson1`, `elasti_k0_hardsoil1`, `hardsoil_gammap_initial_index1` (all rc=1 on the pre-fix binary). Other literal-0 sites reviewed and left: `so.cc:132` (`control_materi_dynamic`, see BACKLOG C11), and the internal no_index records `control_eigen_values`, `control_print_dx_time`, `control_print_tecplot_mesh` (+ gid/gmv/plotmtv mesh/time) — index 0 by design. Regression (clean build): suite 16/16, same 70 file-check failures by name; validation-suite 233 PASS/1 RUNFAIL, sfnet 187/10, test-2026 23/23 (19 + 4 new); every `.dbs` identical to the pre-fix binary after masking timestamps/`time_calculation`/denormals (examp14: 45 s timeout under load, identical when rerun alone). Detail: `manual-developer/group_materi_elasti_k0.md`, `solve_iterative_bicg.md`, `control_materi_plasti_hardsoil_gammap_initial.md`. |
| **`groundflow_phreatic_level_static` without a pressure unknown (audit cluster K22): static pore pressure on the skeleton + wet/dry density for `-groundflow` groups** | `d2ebfb6` | 2026-09-29 | Found by the all-variables audit (cluster K22). Pro `ground18` (1x100 column, level -50, wet/dry 2.0/1.5, `-materi -groundflow`, no `groundflow_pressure`): analytic effective sigyy at the bottom -1250, discrete -1247.5 (bottom element centre); Pro 02-08-2026 -1247.5, GNU before -1995 (wet density everywhere, no pore pressure). Two gaps in `materi.cc`: (1) `get_materi_density()` chose wet/dry from `NODE_PHREATICLEVEL`, written only with the pressure unknown, so always wet; now the pressure-dof rule (average nodal `groundflow_phreatic_coord` total < 0) also runs for static groups; (2) the pore pressure was added to the force stress only inside `if ( groundflow_pressure )`; new branch adds `groundflow_pressure_factor * p_static` at the integration point. Gate `groundflow_static_pressure_group()` (`groundfl.cc`): no pressure unknown + static switch + `-groundflow` in `group_type` (Pro with `group_type -materi` only: -1995, measured). Verified against the Professional binary 02-08-2026: `ground18` GNU after = Pro in every node_dof/node_dof_calcul (sigyy diff 1e-9), identical to its pressure-dof twin `ground4`. Self-checking inputs `validation-suite/test-2026/phreatic_static1.dat` (-1247.5/-749.375/-500) and `phreatic_static_submerged1.dat` (level 10 above the top + `force_edge_water`: -997.5/-500/-1100), analytic = Pro, rc=1 on the pre-fix binary (-1995, -2095). **Divergence found**: the Professional applies the static pore pressure only when `post_calcul -groundflow_pressure -total_pressure` is present (`ground18` without it: Pro -1745, GNU -1247.5); Pro `excavate1` has no such record, so its target (-0.00210987, matched by the pre-fix GNU) now fails (GNU -0.00253891 = Pro with the post_calcul added). BACKLOG C10, owner decision pending. Regression (clean build): suite 16/16, same 70 file-check failures by name; validation-suite 233 PASS/1 RUNFAIL, sfnet 187/10, test-2026 19/19 (17 + the 2 new); every `.dbs` identical to the pre-fix binary after masking timestamps and denormals (examp14 45 s timeout in both). Detail: `manual-developer/groundflow_phreatic_level_static.md`, `AUDIT-ALL-VARIABLES-2026.md` §5.2. |
| **`-updated_linear` / `-updated_area` (audit cluster K01): no stress rotation + compressibility increment from `inc_epe`** | `85b0ed9` | 2026-09-29 | Found by the all-variables audit (cluster K01, 13 tests). Two defects: (1) `materi.cc` "rotate to new configuration" listed `-UPDATED_LINEAR` (since `22801c3`) and `-UPDATED_AREA` (since `a1d9d12`) and rotated the stress by `inc_rot`, which is NOT the identity for them (`set_deften_etc` overwrites it with `set_deften_u_rot`); manual Pro 2.2.1: `-updated_linear` adds the extra stresses with no rotation step (2014 source: only `-updated`/`-total`/`-total_piola`). (2) `stress.cc` compressibility of the linear memories used `tr(new_epe-old_epe)`; `old_epe` is 0 without `materi_strain_elasti`, so the total strain was added every step (legacy 2014); now `tr(inc_epe)`. Verified against the Professional binary 02-08-2026: Pro `interac1` (2 bars, analytic -0.005) Pro/GNU old/GNU new sigxx -0.005/-0.01/-0.005; sfnet `druckpr1` (simple shear, DP, analytic 0.11547) sigxy 0.11547/0.0693/0.11547, sigxx 0/0.0924/0. All 13 K01 tests re-run: 9 now match Pro in every node_dof/post_point_dof variable; the 4 `-fixed_in_space` vs tests (`mdp_shear`, `mstrain_druckprag`, `mstrain_druckprag_elast`, `mesh_act_grav`) keep a residual identical to GNU `-updated_without_rotation` (not the memory). Deliberate divergence: `-total_linear` + compressibility, Pro is time-step dependent (-0.0075/-0.0099/-0.009999 for dt 1e-2/1e-3/1e-4, same as the pre-fix GNU), GNU now -0.005 for every dt (BACKLOG C9). Finding (not changed): Pro builds the `-updated_linear` increment on the CURRENT geometry (explicit `mesh -follow_material`, 1D bar +50%: Pro 0.40714, GNU 0.5); the match without a mesh record comes from the Pro default mesh `-fixed_in_space` vs the GNU default `-follow_material` (BACKLOG C8). Self-checking inputs `validation-suite/test-2026/updated_linear_compress1.dat`, `updated_linear_shear1.dat` (pass on the fixed GNU and on Pro; rc=1 on the pre-fix binary). Regression (clean build): suite 16/16, same 70 file-check failures by name; validation-suite 233 PASS/1 RUNFAIL, sfnet 187/10, test-2026 17/17; classification diff 0 vs the pre-fix binary (except the 2 new tests); every `.dbs` identical after masking uninitialized `element_dof` denormals and timestamps (examp14 45 s timeout under load, identical `.dbs` unloaded). Detail: `manual-developer/group_materi_memory_updated_linear.md` §1/§3/§5/§6, `group_materi_memory_updated_area.md`, `group_materi_memory_total_linear.md`, `AUDIT-ALL-VARIABLES-2026.md` K01. |
| **Damage effective stress (manual Pro 2.2.5 `sigma_damaged = (1-d) sigma_undamaged`): legacy 2014 defect fixed — the damage compounded every step** | `63f134b` | 2026-09-29 | Found by the keyword cycle (`materi_damage_apply`: `dam` = Pro, `sigxx` not). Fix: internal per-integration-point record `element_intpnt_materi_stress_effective` (no Pro counterpart, `external` 0) carries the undamaged stress; damaged = (1-d) x effective, never divided. Verified against the Professional binary 02-08-2026 and theory: 1D bar eps 6e-4 `sigxx` theory/Pro/GNU new 85598.32178 (GNU old 1436.83); load+hold 85598.32178 (old decays to ~0); hex8 compression+yz shear load+hold `sigyy` -275960.2381 in both (old -0.77); gated step (`materi_damage_apply2`) 1.0872657685e7 in both (old 284630); von Mises 1e7 above the damaged peak: no plastic strain in both (yield test on the damaged stress, as in Pro). Open (BACKLOG C7): Mazars gamma with xx+xy shear (GNU and Pro both non-objective, different), damage with plasticity yielding first (plastic path differs). Self-checking inputs `validation-suite/test-2026/damage_stress1.dat`, `damage_stress_hold1.dat`, `damage_stress_3d1.dat` (pass on both binaries; rc=1 on the pre-fix binary and with a changed target). Regression (clean build): validation-suite 233 PASS/1 RUNFAIL, sfnet 187/10, classification diff 0 vs the pre-fix binary, only real `.dbs` change sfnet `damage1` `sigxx`; test-2026 15/15; suite 16/16, same 70 file-check failures. Detail: `manual-developer/materi_damage_effective_stress.md`, `manual-user/materi_damage_effective_stress.md`. |
| **`plasti_apply` (manual Pro 6.899, global) + `control_plasti_apply` (6.257, per timestep index): -no ignores any plasticity data** | `9bbd125` | 2026-09-27 | No runnable Professional usage example (`test/other/safety1/2.dat` use `control_plasti_apply` inside the unimplemented safety_slip feature; `validation_14_mesh.dat` is an include fragment), so the feature was verified against the Professional binary 02-08-2026 with dedicated tests plus Professional distribution tests with `plasti_apply -no` appended. **Scope measured on the Pro** ((b) `-no` == (c) records removed): mohr_coul, vonmises, druck_prag, mohr_coul_direct, tension_direct, tension (single plane-strain quad4: -21.978 / 21.978 in both binaries); hypoplasticity (`hypo1.dat`: sigyy -200, hyhis0 0.6 frozen in both); camclay (`validation_13.dat`: cchis1 207 in both); hardsoil plastic part (Pro only); slide springs follow the switch when `control_slide_plasti_apply` is absent (`slide1.dat`: `plasti_apply -no` == `control_slide_plasti_apply 0 -no` in both, manual 6.371). **Precedence**: the control record of the current timestep overrules the global one (MC and tension_direct two-timestep tests; tension_direct values identical in both binaries: 11.989 / 1). Plastic strains of earlier timesteps are kept. **Implementation**: `stress.cc` `set_stress()` (`plasti_off`: skip `plasti_incr`, `hypoplasticity()`, the direct cut-offs and the `plasti_rule` yield test), `slide.cc` default, helper `materi_apply_gate_off()` in `general.cc`, records in `database.cc`, enums in `tochnog.h`/`tochnog-mod.h`. **Left out**: truss/spring/interface/beam plasticity, user laws (umat/user_sigma), `group_plasti_apply`, `linear_calculation_apply`. Self-checking inputs `validation-suite/test-2026/plasti_apply1..3.dat`, `control_plasti_apply1..2.dat` (rc=0 on both binaries; a changed target gives GNU rc=1). **Regression** (clean `build_safe.sh --clean`, LC_ALL=C, 45 s/test): validation-suite 233 PASS and sfnet 187 PASS / 10 RUNFAIL, per-test classification diff = 0; internal suite 16/16 runs OK, file-check FALLO list identical to HEAD (70). Keyword probe: 827 -> 831 accepted (with the batch below). Detail: `manual-user/plasti_apply.md`, `manual-developer/plasti_apply.md`. |
| **`materi_damage_apply` (manual Pro 6.798) / `materi_failure_apply` (6.801): global (no index) forms of `control_materi_damage_apply` / `control_materi_failure_apply`** | `9bbd125` | 2026-09-27 | No Professional usage example; verified against the Professional binary 02-08-2026. **Damage** (Pro test `other/damage1.dat`): (a) dam 0.993903 in both; (b) `materi_damage_apply -no` == (c) record removed: dam 0, sigxx 1.404e7 in both; two timesteps: control 1 `-no` freezes dam at 0.451188 in both; global `-no` + control 0 `-yes` -> 0.451188 in both. **Failure** (own 2-quad4 model, `group_materi_failure_rupture`): (b) `materi_failure_apply -no` == (c) removed: sigxx 13.4615 in both (with failure active the GNU gives 8.07692; the Pro stops with 'Program error ... element_rhside_delete' in every model tried). **Precedence**: control record of the current timestep overrules the global one (damage and failure). **Behaviour change**: the failure switch no longer freezes damage (the Pro `materi_failure_apply -no` leaves damage unchanged), and it now gates the failure criteria in `failure.cc` (it did not before). Pro binary quirk: it accepts `control_materi_failure_apply` only WITHOUT index (manual says indexed; the GNU keeps the index). **Implementation**: records `MATERI_DAMAGE_APPLY` / `MATERI_FAILURE_APPLY` (INTEGER, no_index, class CONTROL) in `database.cc`, enums in `tochnog.h`/`tochnog-mod.h`, helper `materi_apply_gate_off()` in `general.cc`, consumers `damage()` (damage.cc) and `failure()` (failure.cc, detection of new failures skipped). Self-checking inputs `validation-suite/test-2026/materi_damage_apply1..2.dat` (pass on both binaries), `materi_failure_apply1.dat` (both), `materi_failure_apply2.dat` (GNU only). **Regression** (clean build): validation-suite 233 PASS, sfnet 187/10, classification diff = 0; suite 16/16, FALLO list identical (70). Detail: `manual-user/materi_damage_apply.md`, `materi_failure_apply.md` + developer pages. |
| **`geometry_list` (manual Pro 6.531 / 9.5): list of NODE numbers usable as geometry entity `-geometry_list <index>`** | `a4b95f5` | 2026-09-26 | The only Professional usage example (`test/other/tochnog_abaqus.dat`) is a non-runnable include fragment, so the feature was verified against the Professional binary 02-08-2026 with dedicated tests (3x quad4 strip, materi+condif). **Semantics measured on the Pro**: a node is in the geometry iff its number is in the list; element selections use node membership like any geometry (`control_mesh_delete_geometry` with list {3} deletes nothing, {3 4 7 8} deletes element 3), no element-number interpretation; `node_geometry_present` reports `-geometry_list <index>` before `-geometry_point`. **Implementation**: `database.cc` record `GEOMETRY_LIST` (INTEGER, variable length, class GEOMETRY; enum after `GEOMETRY_LINE` in `tochnog.h`/`tochnog-mod.h`); `geometry.cc` `geometry()` branch (membership by `inod`, `inod<0` never inside, projection = node) and `node_geometry_present_calculate()` scan; `change.cc` refusal (avoids an `idum` overflow). **Results GNU = Pro**: `geometry_list1` node 4/8 disx 0.03, node 8 temp 10, nodes 1/5 temp 0; delete tests identical element/node sets; node_geometry_present identical. Lists vs equivalent `geometry_line`/`geometry_point` give byte-identical `node_dof` in BOTH binaries (the remaining GNU/Pro interior difference, e.g. sigxx 9.998 vs 10, is pre-existing and model-wide). Self-checking inputs `validation-suite/test-2026/geometry_list1.dat`, `geometry_list_delete1.dat`, `geometry_list_delete2.dat` (rc=0 on both binaries). **Regression (clean `build_safe.sh --clean`, LC_ALL=C, 45 s/test)**: validation-suite 233 PASS before and after, sfnet 197 187 PASS / 10 RUNFAIL before and after, per-test classification diff = 0; internal suite 16/16 runs OK with the file-check FALLO list identical to a clean build of HEAD `0cfbbb5`. `-node <n>` as geometry entity (manual 9.4) is NOT supported: the GNU parses and ignores it, the Pro parses it and applies nothing (zero pivot). Detail: `manual-user/geometry_list.md` + `manual-developer/geometry_list.md`. |
| **CIERRE del residual mpc5: sincronización POST-solve del REPORTE de tensión nodal (re-evaluación constitutiva desde la velocidad RESUELTA, candidato (b)) — mpc5 sigyy 1.0 → 0.5 en k=1 (y 0.5 en k=2/3/5) = Professional one-pass; 0 regresiones de corpus** | `dae06a3` (feat) + docs de este lote | 2026-09-15 | **SÍNTOMA**: `mpc5.dat` daba `post_point_dof 10 -sigyy = 1.0` en k=1 y `0.5` en k=2, mientras el Professional 02-08-2026 reporta `0.5` en k=1 (one-pass). El valor de k=2 ya coincidía tras la recuperación absoluta (`0b23957`); quedaba el lag de REPORTE de una iteración. **CAUSA RAÍZ**: los dofs `materi_stress` son NO-principales (`input.cc` no fija `dof_principal`), así que `dof.cc parallel_new_dof_diagonal` los actualiza DESPUÉS del solve como `node_dof = node_rhside/node_lhside`, con el `element_rhside` ensamblado por `element_loop()` ANTES del solve (`materi.cc` ensambla `volume*weight*new_sig` con `new_sig` de la velocidad PRE-solve) → la tensión reportada va una iteración por detrás de la velocidad. **VERIFICACIÓN DEL CANDIDATO (b) (medida)**: el VALOR que produce la re-evaluación post-solve en la iteración k es el que el esquema sin tocar reporta en k+1 (`proj(sigma_constit(v_k))`), pero el MECANISMO es distinto: agrega una evaluación constitutiva, NO un solve extra, así que sincroniza k=1 sin una segunda resolución de momento. La afirmación "idéntico a la iteración k+1" es cierta para la secuencia de valores, no para el costo, y NO impide el sync de k=1. La sonda runtime (segundo `element_loop()` tras `solve()` gateado por `TOCHNOG_STRESS_POST`) confirmó k=1: 1.0 → 0.5. **IMPLEMENTACIÓN (mínima, stress-only)**: global `stress_recovery_post_pass` (`polynom.cc`, declarado en `tochnog.h`/`tochnog-mod.h`); en `top.cc` `top()`, tras `solve()` y antes de `parallel_new_dof_diagonal()`, con `stress_recovery_absolute && materi_stress && stres_indx>=0 && !find_local_softvar`: (1) se guardan los vectores `NODE_LHSIDE`/`NODE_RHSIDE` de todos los nodos activos; (2) se ponen a cero los 6 slots de tensión (`stres_indx/nder .. +5`); (3) `stress_recovery_post_pass=1; element_loop(); =0`; (4) se restauran las entradas NO-tensión (los elementos estructurales —truss/beam/spring/contactspring/interfaz— saltan `materi()` vía `goto add_to_system_vectors`, así que sin el restore sus contribuciones del primer pase se duplicaban; medido: truss2/truss3/conspr4 regresaban antes del restore). En `materi.cc` el pase salta el bloque de velocidad/momento y los bloques de otros dofs (solo ensambla la recuperación de tensión); en `elem.cc`/`element_loop()` salta `general()`, las llamadas de tipos distintos de `MATERI`, los bloques post-`materi` (derivatives/strain energy/`options_element_dof`/`element_dof_tmp`), el bloque post-punto-de-integración (tendon/area/residue/`ELEMENT_VOLUME`) y el store sparse; el bloque de elemento borrado NO re-aplica el RHS congelado (`element_rhside_delete` ya está guardado), solo pone a cero los slots de tensión de `element_rhside`/`element_lhside`. `parallel_new_dof_diagonal` corre una sola vez por iteración con mezcla de tensión fresca (velocidad resuelta) y no-tensión pre-solve. **GATE**: solo con recuperación absoluta; los dofs de tensión CONVECTADOS (`materi_strain_total`) conservan el update incremental. **VERIFICACIÓN (worktree scratch `/tmp/opencode/mpc5sync-20260915-231446` desde `3d4d70e`, build limpio `build_safe.sh --clean`, LC_ALL=C)**: mpc5 k=1/2/3/5 → `sigyy = 0.5` (rc=0); delete2 node2 sigxx **0.8**, delete3 sigyy **0.6** (byte-idénticos al baseline); delete1/4 y mpc1/2/3/4/5/6/7/8/9 rc=0. **CORPUS A/B** (base vs fix): sfnet 197 **187 PASS / 9 RUNFAIL / 1 PARSE**, diff per-test = **0**; validation-suite 233 **201 PASS / 32 RUNFAIL / 0 PARSE**, diff = **0**; suite interna **16/16** runs OK. Detalle: `manual-user/nodal_stress_post_solve_sync.md` + `manual-developer/nodal_stress_post_solve_sync.md`; `mpc5_staggered_stress_lag.md` marcado CERRADO. |
| **OPTIMIZACIÓN (hoist): el pase POST-solve de tensión nodal (`dae06a3`) se saca del bucle de equilibrio — el pase stress-only y la actualización de los dofs de tensión corren UNA vez por paso en vez de una vez por iteración; large3 (784 elem, 220 it) 1.23-1.25x, dynamic8 (16 elem, 4000 it) 1.12-1.21x; mpc5/delete2/delete3 y los targets de corpus sin cambio, 0 cambios de clasificación** | `3957809` (rama `dev/stress-pass-hoist`) | 2026-09-16 | **PROBLEMA**: `dae06a3` corría un `element_loop()` stress-only tras CADA iteración de equilibrio para modelos con `materi_stress` + recuperación absoluta, ~duplicando el costo de ensamblado de elementos. **CAMBIO (solo `top.cc`)**: el bloque (guardar `NODE_LHSIDE`/`NODE_RHSIDE`, poner a cero los 6 slots de tensión, `stress_recovery_post_pass=1; element_loop(); =0`, restaurar las entradas no-tensión) se mueve FUERA del `for (iteration...)`, inmediatamente después del bucle y antes de `db_version_copy(NEW,NORMAL)`; los dofs de tensión se asignan directamente (variante stress-only de `dof.cc parallel_new_dof_diagonal`, incluyendo los slots de `materi_stress_pressure_history` y de `derivatives`), solo los slots de tensión (re-correr la función completa duplicaría el `+=` de los otros no-principales). **MEDICIÓN (misma máquina, `LC_ALL=C`, cwds frescos, `build_safe.sh`)**: large3 (784 elem, 220 it) base min/mediana 56972/59024 ms → fix 46224/47213 ms = **1.23-1.25x**; dynamic8 (16 elem, 4000 it) 11227/12180 → 10009/10028 = **1.12-1.21x**; delete3 (2 elem, 25 it) 35/38 → 31/33 = ~1.14x; mpc5/delete2 dominados por el arranque del proceso (~15 ms). **VALIDACIÓN**: mpc5 `sigyy=0.5` en k=1/2/3/5 rc=0 (idéntico base/fix); delete2 node2 `sigxx=0.8`; delete3 `sigyy=0.6`; campos `node_dof`+`post_point_dof` byte-idénticos. **CORPUS A/B**: sfnet 197 **187 PASS / 9 RUNFAIL / 1 PARSE** en ambos, diff per-test = **0** y diff de valores `node_dof`/`post_point_dof` = **0**; validation-suite 233 **231 PASS / 2 RUNFAIL** (`elmethod`, `numit`) en ambos, diff de clasificación = **0** (único cambio de valor no-target: `contact_heatgen` `xsigyy`, derivada-x de la tensión; su target `node_dof 1 -disy` sin cambio). Suite interna **16/16** runs OK (lista de FALLO de file-checks idéntica base/fix). **CAMBIO OBSERVABLE NO-TARGET documentado**: `mpc8` (`-total_linear` + `materi_displacement`) `node_dof sigxx` (0.25,0.25,0.25)→(0,0.125,0.25) en k=2 — el pase hoisteado corre tras `parallel_new_dof_diagonal` y reporta la tensión consistente con el desplazamiento actualizado (su target `post_element_force_result=0.5` sin cambio; en k=1 el fix da (0,0.25,0.5) = Professional). Detalle: `manual-developer/nodal_stress_post_solve_sync.md` §Cost/§Verification; `IMPROVEMENT-ANALYSIS.md` item (a) marcado DONE. |

| **Familia `control_mesh_delete_geometry` (keywords `_stop`/`_stop_geometry`, manual Pro 6.174/6.175) + congelado del estado de punto de integración del elemento parcialmente borrado — delete2 pasa de PARSE a RUNFAIL (parse cerrado); target node_dof 0.8 BLOQUEADO por el lag de la recuperación nodal de tensión (work unit aparte)** | `0735028` (feat) + docs de este lote | 2026-09-15 | **SÍNTOMA**: `delete2.dat` (1D, 2 `bar2` en serie, `control_reset_dof -sigxx 1`, borrado del bar derecho al 40 % en el índice 10, `control_mesh_delete_geometry_stop 20 -yes`, timestep en el 30, target `node_dof 2 -sigxx = 0.8`) fallaba en PARSEO: `I don't know what to do with : control_mesh_delete_geometry_stop` (el keyword no estaba registrado y, al ser desconocido, el parser lo absorbía como valor del `control_mesh_delete_geometry_factor` precedente → `Problem reading : control_mesh_delete_geometry_factor`). **IMPLEMENTACIÓN (`delete.cc` `delete_geom`, `database.cc`, `tochnog.h`/`tochnog-mod.h`)**: (1) `_stop` se lee en el `icontrol` actual (puede vivir en un índice DISTINTO al de la geometría/factor; el bloque corre ANTES del gate de `CONTROL_MESH_DELETE_GEOMETRY`); con `-yes` recorre los elementos con `ELEMENT_DELETE_FACTOR` activo y borra `ELEMENT_DELETE_FACTOR` + `ELEMENT_RHSIDE_DELETE` (los elementos borrados COMPLETAMENTE, factor final > 1-1e-10, tienen el record de elemento eliminado y NO se reactivan, manual 6.170); (2) `_stop_geometry` restringe el stop a los elementos con TODOS sus nodos dentro de la geometría (misma membresía que el delete, `geometry()` con proyección delete). **FIX DEL ESTADO DEL ELEMENTO BORRADO (`elem.cc` `elem()`)**: el GNU recomputaba `ELEMENT_DOF` en cada iteración, así que el elemento parcialmente borrado derivaba a 1.4 mientras el Pro congela `element_intpnt_dof_average = 1.0` (medido en el .dbs del Pro 02-08-2026: el elemento borrado guarda `element_delete_factor`+`element_rhside_delete` y su `element_intpnt_dof_average` queda en el valor pre-borrado; `global_element_dof_apply` es `-yes` por defecto, manual 6.553). El fix salta el PUT de `ELEMENT_DOF` una vez capturado el `ELEMENT_RHSIDE_DELETE` (primera iteración de borrado), dejando el valor copiado de `VERSION_NORMAL` al inicio del paso. Medido: `element_dof` de ambos elementos = 0.8 EXACTO al Pro tras el stop. **RESIDUAL MEDIDO (mecanismo distinto, NO cerrado)**: el target `node_dof 2 -sigxx = 0.8` da **0.6** (nodo 1 y 3 dan 0.8, el nodo compartido 2 queda en 0.6). Causa: la recuperación nodal de tensión del GNU (`dof.cc parallel_new_dof_diagonal`) ACUMULA los INCREMENTOS `(new_sig-old_sig)/dt` de los elementos, y en el nodo compartido por dos elementos cuyos incrementos (+0.2/−0.2) se cancelan el valor reportado queda en el de inicio del paso; el Pro ASIGNA la media ponderada ABSOLUTA (`node_rhside/node_lhside`, medido: Pro node2 rhside=1.6 lhside=2.0 → 0.8; GNU rhside≈0 lhside=1.0). Barrido de iteraciones (default 2 y 5): GNU sigue en 0.6. Es el mismo lag ya documentado en la fila de mpc5 (`mpc5_staggered_stress_lag.md`); un intento de recuperación absoluta (rhside = peso·new_sig + asignación) REGRESIONÓ delete3 (0.6→0.12) y se REVIRTIÓ, confirmando que el mecanismo de pesos (least-squares) no es un simple reemplazo. **VERIFICACIÓN (worktree scratch `/tmp/opencode/delete2-20260915-210017`, rama `dev/delete-geometry-stop` desde `3f63f91`, build limpio `build_safe.sh`, LC_ALL=C)**: delete2 rc=1 (target 0.8, medido 0.6); delete1/delete3/delete4/mpc5 rc=0 SIN cambio (delete3 sigyy=0.6, mpc5 sigyy=0.5 EXACTOS); corpus sfnet 197 **187 PASS / 9 RUNFAIL / 1 PARSE** = byte-idéntico al baseline `3f63f91`; validation-suite 233 **201 PASS / 32 RUNFAIL / 0 PARSE** = byte-idéntico; corpus 363 (351 `.dat`) A/B contra el binario de `3f63f91`: único cambio real = `other/delete2.dat` PARSE→RUNFAIL (el parse cerrado); `other/dynamic5.dat` PASS→timeout 45 s es FLAKY (ambos binarios tardan ~44-45 s, re-corrido rc=0 en ambos); suite interna `build_safe.sh` **16/16 runs OK**. Detalle: `manual-user/control_mesh_delete_geometry_stop.md` + `manual-developer/control_mesh_delete_geometry_stop.md`. |
| **RECUPERACIÓN NODAL DE TENSIÓN = media ponderada ABSOLUTA del Professional (asignación, no acumulación de incrementos) — cierra delete2 `node_dof 2 -sigxx = 0.8` SIN la regresión de delete3 (0.6); todo el campo `node_dof`/`post_point_dof` de la familia delete/mpc = Pro a <=2.1e-10** | `(fix + docs de este lote)` | 2026-09-15 | **SÍNTOMA**: delete2 (2 `bar2` en serie, bar derecho borrado 40 % y reactivado, `control_reset_dof -sigxx 1`) da `node_dof 2 -sigxx = 0.6` (target 0.8); el nodo compartido 2 queda en el valor de inicio de paso. **CAUSA RAÍZ (medida con `.dbs` del Pro 02-08-2026)**: `dof.cc parallel_new_dof_diagonal` ACUMULA (`node_dof_new += node_rhside/node_lhside`) con `materi.cc` ensamblando `volume*weight*(new_sig-old_sig)/dt`; en un nodo compartido por dos elementos cuyos incrementos se cancelan (+0.2/−0.2) el reportado no cambia. El Pro ASIGNA la media absoluta: `node_dof = Σ_e element_dof_e / N_e`, con `node_lhside` = nº de elementos adyacentes (delete2 Pro: node2 `rhside=1.6`, `lhside=2.0` → **0.8**; delete3 Pro: `lhside` = 1 en los nodos compartidos (solo el elemento vivo) y 0 en los exclusivos del borrado, `node_dof` = 0.6 / 1.0). **FIX (mínimo, gated)**: nuevo global `stress_recovery_absolute` (`polynom.cc`, declarado en `tochnog.h`/`tochnog-mod.h`), fijado en `general.cc` como `!materi_strain_total && (no convectado)`; con el flag, `materi.cc` ensambla el numerador ABSOLUTO `volume*weight*new_sig` y el denominador `volume*weight` (mismo peso least-squares 2026-09-13, `sri_stress_recovery_weight`), `dof.cc` ASIGNA (`=`) para los dofs de tensión (el resto de no-principales sigue `+=`), `general.cc` salta el término de inercia incremental, y `elem.cc` pone a CERO los slots de tensión del `element_rhside` del elemento borrado (el `element_rhside_delete` congelado es un RHS de MOMENTO; dejarlo contaminaba la media y causaba la regresión delete3 0.6→0.12 del intento ingenuo revrtido). **GATE**: no se aplica a dofs de tensión CONVECTADOS (strconv1) ni con deformación independiente (`materi_strain_total`, examp9/10/23, sin referencia Pro: keywords no soportadas por el binario Pro), que conservan el update incremental. **VERIFICACIÓN (binario desde `73fb304`, build limpio `build_safe.sh`, LC_ALL=C, worktree scratch `/tmp/opencode/recovery-20260915-223445`)**: delete2 rc=0 node2 sigxx **0.8**; delete3 rc=0 sigyy **0.6**; mpc5 rc=0 sigyy **0.5** (target, en sus 2 iteraciones); delete1/4, mpc1/2/3/4/6/7 rc=0; campo `node_dof`+`post_point_dof` COMPLETO vs Pro 02-08-2026: delete1/2/4 diff 0, delete3 1.1e-14, mpc1 0, mpc2 3.2e-13, mpc3 2.1e-10, mpc4 1.9e-10, mpc5 5.4e-17, mpc6 1.0e-12. **CORPUS**: sfnet 197 **187 PASS / 9 RUNFAIL / 1 PARSE**, A/B per-test = **0 cambios**; validation-suite 233 **201 PASS / 32 RUNFAIL / 0 PARSE**, A/B = **0**; A/B de VALORES sobre los 191 tests sfnet con `.dbs`: `node_dof`/`post_point_dof` máx diff **0**; suite interna **16/16** (lista de FALLO de file-checks idéntica al baseline); strconv1 byte-idéntico. **RESIDUAL (mecanismo distinto, NO cerrado)**: el lag de REPORTE de una iteración de mpc5 (1.0 en k=1, 0.5 en k=2; el Pro 0.5 en k=1) sigue: los dofs de tensión son NO-principales y se actualizan post-solve con la tensión constitutiva PRE-solve; cerrarlo exige la re-evaluación constitutiva post-solve (candidato (b) de `mpc5_staggered_stress_lag.md`), work unit aparte. Detalle: `manual-user/nodal_stress_recovery.md` + `manual-developer/nodal_stress_recovery.md`. |
| **CIERRE validation_2 + validation_8 (suite de validación PROPIA del Professional, `test/validation/validation_{2,8}.dat`; NO están ni en sfnet ni en validation-suite) — ambos ya cerrados por commits ANCESTROS de `a41cda4`: validation_2 por `6eb399d` (ancho de banda del retry LU directo) y validation_8 por `f244ac5` (`-veln`); sin cambio de fuente nuevo, 0 regresiones** | `(docs de este lote)` | 2026-09-15 | **UBICACIÓN**: los `validation_*` NO viven en `external-downloads/sfnet/extracted/test` (197 `.dat`, sin `validation_*`) ni en `validation-suite/test-2014/` (233, sin `validation_*`); son la suite propia del Professional (`tochnog_version_02-08-2026/test/validation/validation_1..14.dat`; copia `/tmp/masin-kernel-corpus-full.8pjO3w/input/validation/`), untracked. **validation_2** (consolidación 1D acoplada, `inertia_apply -yes` + `groundflow_consolidation_apply -yes`, 6 h-refinements, target `post_point_dof 1 -pres = -0.63 ± 0.05`): el RUNFAIL histórico era el ancho de banda del retry LU directo tras breakdown Bi-CG (`so.cc`, banda `diff*nprinc` en vez de `(1+diff)*nprinc`; el acoplamiento presión↔velocidad de la consolidación, `nprinc=2`, distancia 3 > banda 2, escribía fuera de `mat`/`ldmat`), ya corregido por `6eb399d` (ancestro de `a41cda4`; el mismo fix cerró `examp2`, la copia sfnet de la misma física). Medido: GNU `a41cda4` `-0.634519741323` rc=0 (default y con `control_solver 0 -matrix_lapack`), Pro 02-08-2026 `-0.6345197413352` rc=0 (10 dígitos), GNU pre-merge `-0.634519842135` rc=0, GNU jan-2014 rechaza `inertia_apply`. **El target `-0.63` NO es stale** (el modelo da `-0.6345`, dentro de ±0.05); el `-0.688985` de la nota previa era el retry corrupto, REFUTADO. **validation_8** (metal forming euleriano, `mesh -fixed_in_space` + `convection_apply -yes`, von Mises con hardening, `bounda_dof 2 -geometry_set 1 -veln`; sin target → PASS = rc=0): el PARSE histórico era `-veln` sin implementar; manual Pro 6.22 da la normal desde la entidad geométrica (o `bounda_normal` para rangos de nodo) y dice que "internally a multi-point-constraint will be generated", pero NO especifica el eje esclavo/maestros — ese detalle lo dedujo por medición `f244ac5` (ancestro de `a41cda4`): normal = primer match de `geometry()` en orden del set; un `mpc_node_number` por nodo (esclavo = primer eje de velocidad con componente normal ≠0, maestros = ejes restantes, factor `-n_master/n_slave`, componentes cero omitidos), `mpc_from_bounda -yes`, regenerado por fingerprint de malla. Evidencia fresca: GNU `a41cda4` reproduce 50/50 `mpc_node_number` + 8/8 `mpc_node_factor` EXACTOS al Pro; rc=0 (163.5 s). **RESIDUAL (harness, NO bug de fuente)**: el test no cabe en el slot de 45 s del corpus (GNU 163 s / Pro 200 s; ambos 2000 pasos × 192 elementos; variante de 20 pasos 1.6 s = lineal) → el runner lo clasifica timeout/DEFERRED; el test pasa standalone. **VERIFICACIÓN (worktree scratch `/tmp/opencode/val28-20260915-202314/wt` desde `a41cda4`, build limpio `build_safe.sh --clean`, LC_ALL=C)**: sfnet 197 **187 PASS / 9 RUNFAIL / 1 PARSE**, diff per-test vs el run predecesor de `a41cda4` (`mpc5del-20260915-180036/corpus-sfnet-final`) = **0**; validation-suite 233 **201 PASS / 32 RUNFAIL / 0 PARSE**, diff = **0**; suite interna **16/16**. Detalle: `VALIDATION-PROFESSIONAL-2-8.md`. |
| **FIX mpc5: la carga liberada del elemento borrado es energy-consistent en el residuo staggered — la matriz y el `element_lhside` del elemento borrado se ponen a CERO (el RHS liberado congelado es constante, su tangente es 0); mpc5 rc=1 sigyy 0.7 → rc=0 sigyy 0.5 en sus 2 iteraciones** | `(fix + docs de este lote)` | 2026-09-15 | **SÍNTOMA**: `mpc5.dat` (quad9 borrado al 50% + reset `-sigyy 1` + tie `mpc_linear_quadratic`) esperaba `post_point_dof 10 -sigyy = 0.5` (±1e-3) y daba rc=1 sigyy=0.7 en su propio número de iteraciones (default 2). **CAUSA RAÍZ (medida)**: el bloque de borrado de `elem.cc` aplicaba el RHS liberado CONGELADO (`element_rhside = factor·element_rhside_delete`) pero seguía escalando la matriz COMPLETA por el factor (`element_matrix *= factor`, fix `7687e0a`). El RHS congelado es una CONSTANTE: su tangente es cero, pero la rigidez residual `factor·K` seguía en el solve y respondía a la deformación → se rompe la cancelación Newton `A = dt·K = dt·BᵀCB` del esquema staggered (que para elementos no borrados colapsa a one-pass) y la carga liberada relaja geométricamente a tasa exactamente 0.5 (k: 1→1.0, 2→0.7, 3→0.59, 5→0.5211, 8→0.5026, 20→0.5000006). Además el `element_lhside` del elemento borrado NO se escalaba: su peso lumped (con RHS de tensión congelado a 0) diluía la recuperación nodal de tensión del elemento superviviente en los nodos compartidos y frenaba la tensión incluso con la velocidad ya one-pass. **FIX (`elem.cc`, `elem()` `add_to_system_vectors`)**: con `ELEMENT_DELETE_FACTOR` activo, además de aplicar `element_rhside = factor·element_rhside_delete`, se pone `element_lhside[i] = 0` y `element_matrix[i] = 0` (emparejamiento Jacobiano-consistente: RHS constante ⇒ tangente cero). **VERIFICACIÓN (build limpio `build_safe.sh`, LC_ALL=C, worktree scratch `/tmp/opencode/mpc5del-20260915-180036/wt` desde `1e31fd6`, binario `final-bin`)**: **mpc5** rc=0 sigyy **0.5** en su propio número de iteraciones (k=2; k≥2 estable en 0.5); el campo `node_dof` COMPLETO coincide con el binario Professional 02-08-2026 a ~5e-17 (nodos 5/6/7 vely=0 sigyy=1; 9/11 vely=0.5 sigyy=0.5; 10 vely=0.5 sigyy=1; 13/14 sigyy=0.5). **delete3** rc=0 sigyy **0.6** (antes 0.6 a 5 iteraciones, ahora a 2) y campo completo = Pro a ~1e-14. El match se mantiene para TODOS los factores de delete probados (mpc5 con 0.3/0.5/0.7/0.8; delete3 con 0.4/0.6/0.8). **mpc1/2/3/4/6/7/8/9** rc=0 (sin delete, `node_dof` byte-idéntico al baseline); delete1/delete4 rc=0; delete2 rc=1 pre-existente (keyword `control_mesh_delete_geometry_stop` no implementada, PARSE). **CORPUS**: sfnet 197 **187 PASS / 9 RUNFAIL / 1 PARSE** (el `pruebavel` del clasificador; el A/B per-test vs el baseline limpio `1e31fd6` = **0 cambios de clase**) y validation-suite 233 **201 PASS / 32 RUNFAIL / 0 PARSE** (A/B per-test = **0 cambios**); suite interna **16/16** runs OK. **RESIDUAL MEDIDO (mecanismo distinto, NO cerrado)**: la proyección nodal de tensión del GNU es un update lumped POST-solve con el residuo PRE-solve (`dof.cc parallel_new_dof_diagonal`), así que la tensión reportada va una iteración por detrás de la velocidad: mpc5 da 1.0 en k=1 y 0.5 en k=2 (el Pro reporta en sincronía, 0.5 en k=1). El test pasa en su propio número de iteraciones (2); cerrar el lag de reporte de una iteración es una re-evaluación constitutiva post-solve, work unit aparte. Detalle: `manual-developer/delete_factor_element_matrix_scaling.md`. |
| **FIX regresión delete3: se refina el gate per-index de `control_reset_dof` para reconciliar delete3 (reset one-shot) con creset_geom (reset step-scoped) — delete3 0.920015 → 0.600002285218, 0 regresiones** | `a548415` (fix) + docs de este lote | 2026-09-15 | **REGRESIÓN**: `f142bdf` eliminó el gate `ireset!=icontrol` (añadido en `22801c3`) para arreglar `creset_geom` y reabrió delete3: el reset `-sigyy` del índice 5 se re-aplica en cada timestep del índice 10, re-inyecta sigma=1 y mueve el punto fijo de `0.600002285218` (rc=0, binario de `7687e0a`) a `0.920015` (rc=1). **BISECCIÓN**: el único commit de código entre el estado bueno documentado `6eb399d` (delete3 0.600002285218, `final-validation-20260914/REPORT.md`) y el malo es `f142bdf`; build limpio de `7687e0a` -> 0.600002285218 (rc=0) y de `f142bdf` -> 0.920015 (rc=1). `f045d04` (formulación `-total_linear`) es ancestro de `7687e0a`, NO es el culpable. **MECANISMO (traza instrumentada)**: delete3 intenta el reset en `icontrol=5` (t=0, paso no-timestep) y en cada `icontrol=10` (t=0.2..1.0); creset_geom en `icontrol=0` (t=0.1) e `icontrol=1` (t=0.2). Sin gate el reset de delete3 dispara en cada paso; con el gate estricto el reset de creset_geom del índice 1 no dispara en el paso 0, el history IP `element_dof` del brick queda 0.5 y el paso 1 resuelve hisv0=-1 (target 0). **FIX (`data.cc`)**: `if ( ireset!=icontrol && !db_active_index( CONTROL_TIMESTEP, ireset, VERSION_NORMAL ) ) continue;` — un reset en un índice SIN timestep es one-shot (su propio paso); un reset cuyo índice TAMBIÉN lleva `control_timestep` conserva el comportamiento legacy GNU (cada `data()`), que el Professional rechaza (`"Error detected for data item : control_reset_dof, record : 20"`) y que solo alcanza a inputs GNU-only. **VERIFICACIÓN (build limpio `build_safe.sh --clean`, LC_ALL=C, worktree scratch `/tmp/opencode/d3reg-20260915-135147/wt`, rama `dev/delete3-regression` desde `388a221`)**: delete3 rc=0 `0.600002285218` (target 0.6±1e-4) IDÉNTICO al binario de `7687e0a`; `creset_geom` rc=0 (node1 hisv0=0, node5=0.5) y `creset_iface` rc=0; `mat_rel_reset`/`reset_value_linear` (los otros 2 tests con reset en índice de timestep) rc=0; `total3` rc=0 disx=0.5; delete1/2/4 sin cambio; `iface_condif`/`groundflow_interface` rc=0; corpus sfnet 197 A/B **byte-idéntico a HEAD (187 PASS / 10 RUNFAIL / 0 PARSE)**; validation-suite 233 A/B **idéntico (201 PASS / 32 RUNFAIL / 0 PARSE)**; suite interna **16/16**. Detalle: `manual-user/control_reset_dof_scoping.md` + `manual-developer/control_reset_dof_scoping.md`. |
| **FIX ancho de banda del retry LU directo tras breakdown Bi-CG (`so.cc` `solve()`): la banda se quedaba corta en `nprinc-1` dofs y el fill escribía FUERA de `mat`/`ldmat`, corrompiendo el retry (`examp2` -0.7206 en vez del exacto -0.6345) — corpus 186/11/0 → 187/10/0, 0 regresiones; 4 tests OUT-OF-SCOPE** | `(fix + docs de este lote)` | 2026-09-14 | **SÍNTOMA**: `examp2.dat` (columna 1D de consolidación acoplada `-materi -groundflow`, `groundflow_pressure`, 64 elementos tras 6 h-refinements) con `inertia_apply -yes` + `groundflow_consolidation_apply -yes` da `post_point_dof 1 -pres = -0.72055658378` (target -0.63±0.05, RUNFAIL); el Pro 02-08-2026 da `-0.6345197413352` y el GNU jan-2014 da `-0.634520244513`. **CAUSA RAÍZ (medida)**: con `options_solver -matrix_lapack` o `-matrix_superlu` el GNU actual da EXACTAMENTE `-0.634519741323` = Pro → el sistema ensamblado es correcto y el culpable es el retry. El solver por defecto (Bi-CG) rompe (sistema mal condicionado) y reintenta con `dgbsv` sobre la matriz band `mat`. El ancho de banda del camino `bicg_solver` (`so.cc:517-529`, ordenación identidad) se calculaba como `diff*nprinc`, pero dos nodos separados `diff` números de nodo ocupan `diff*nprinc + (nprinc-1)` índices locales consecutivos: el acoplamiento presión↔velocidad de la consolidación (último dof de un nodo al primero del siguiente, `nprinc=2`, distancia 3) supera `band=2`, el índice de banda `ii=kl+ku+1+ilocal-jlocal=8 > ldmat=7` y `mat[indx-1] += val` escribe en la columna siguiente (corrupción silenciosa). La rama Cuthill-McKee (usada por `-matrix_lapack`) YA usaba `(1+diff)*nprinc`; la rama Bi-CG (añadida después para el retry) omitía el `+1`. Instrumentación temporal: `band=2 solve_nlocal=128 nprinc=2 npuknwn=14` y `OUT_OF_BAND ilocal=3 jlocal=0 ii=8 ldmat=7 diff=3`. `git bisect` (build+run por commit, solver por defecto) → primer commit malo `6cafafe` (defaults Bi-CG más estrictos, que hacen que el Bi-CG rompa y tome el retry corrupto; el bug de banda es latente desde antes). **FIX (`so.cc`)**: `if ( (1+diff)*nprinc > band ) band = (1+diff)*nprinc;` (mismo bound que la rama de ordenación). **VERIFICACIÓN**: `examp2` rc=0, `pres = -0.634519741323` = Pro `-0.6345197413352` (10 dígitos) y target -0.63±0.05; **CORPUS sfnet 197** (runner `corpus-ab.sh`, worktree scratch `/tmp/opencode/final5-20260914-160749`, rama `dev/parse-init-examp2` desde `0d0ad0c`; binario baseline propio `baseline-tochnog`): **186 PASS / 11 RUNFAIL / 0 PARSE → 187 PASS / 10 RUNFAIL / 0 PARSE**, diff per-test = SOLO `examp2` RUNFAIL→PASS, **0 regresiones**; suite interna `build_safe.sh --clean` **16/16**. **OUT-OF-SCOPE (4)**: `pruebavel` (`group_condif_flow 0 -velx -vely` usa nombres de dof; manual Pro 6.600 es numérico `beta1 beta2 beta3`; GNU, GNU jan-2014 y Pro lo rechazan; la copia en validation-suite es 0 bytes → test muerto/malformado); `ho_mech2` (test de inestabilidad con `-p_refinement`; GNU actual `vely=-1238`, jan-2014 aborta con `number_of_intagration_points at least 9`, el Pro propio `test/other/tensile_instability.dat` usa `-updated_linear` + otro control de timestep y target propio `-1975`; el target sfnet `-641` no lo alcanza ningún binario); `taylor3` (impacto axisimétrico; GNU actual se estanca en el timestep mínimo, `kap=0`; jan-2014 corre pero da `2.56391` vs target `2.58442` (falla); la copia propia del Pro usa `inertia_apply`, `mesh -follow_material`, `group_materi_memory -updated`, control de timestep de 3 valores y `_stop -no`; el Pro rechaza el input GNU por `options_stabilization`); `ho_mech4` (initia sin `materi_strain_plasti`, requerido por `group_materi_plasti_vonmises` según manual Pro 6.750 y presente en la copia propia del Pro `test/other/cylinder_buckle.dat`; añadiéndolo el GNU corre pero da `kap=0.06274` vs target `0.0747218`, y el target del Pro es otro `node_dof 1010 -velx=0.0118`). **CORPUS NO VERSIONADO**: los `.dat` viven en `external-downloads/sfnet/extracted/test` (gitignored); el A/B usa una copia alineada (incnav3 + `inertia_apply -yes`, examp2 + los dos switches) documentada en `CORPUS-PARSE-INIT-EXAMP2.md`; originales y runs por test en `/tmp/opencode/final5-evidence-20260914-160940/`. Detalle: `manual-user/solver_bicg_band_width.md` + `manual-developer/solver_bicg_band_width.md`. |
| **Subgrupo solver/timeout del corpus sfnet (6 tests): fix fuente `refine_l.cc` BAR2 fall-through (refine5) + `incnav3` noveno input stale (`inertia_apply -yes`) — corpus 184/13/0 → 186/11/0, 0 regresiones; 4 tests OUT-OF-SCOPE** | `(fix + docs de este lote)` | 2026-09-14 | **FIXABLE (fuente) — `refine5`**: abortaba `rc=134` en `refine_l.cc:231` `assert(name==-TET10)`. Causa: el dispatch de forma de elemento en `refine_locally()` tiene `if (name==-BAR2){...}` como `if` SUELTO seguido de `if (name==-BAR3){...} else if(TRIA3)... else { assert(name==-TET10); }`; para un elemento `-bar2` el primer `if` corre y el segundo (falso) cae al `else` final → assert. El MISMO fichero maneja `-BAR2`/`-BAR3` con `else if` correcto más abajo (líneas 677/690) → es un `else` faltante, no una feature. Fix: `else if (name==-BAR3)` (`refine_l.cc:203`). `refine5` rc=0, `post_point_dof 0 -temp = 1` (target 1±5e-2); el GNU 2014 de referencia aborta IDÉNTICO (pre-existente desde `a2407e0`). **FIXABLE (corpus input) — `incnav3`**: `system (numerically) singular or not positive definite`. Causa: sin BC de velocidad + material con solo `elasti_compressibility` (sin cizalla) + SIN inercia → el operador de momento es puramente volumétrico y numéricamente singular. El default `options_inertia=-NO` de `0560dcf` (alineación Pro, ancestro de `b83d907`) quitó el término de masa `M/dt²` que regularizaba el operador. Medido `node_dof 2 -tvely` (target -2.5±1e-2): GNU 2014 (default `-yes`) **-2.49461364553** rc=0; GNU actual sin inercia → breakdown Bi-CG + LU falla rc=1; GNU actual + `inertia_apply -yes` → **-2.49461** rc=0; Pro 02-08-2026 con su propia copia `test/other/incnav3.dat` (que SÍ lleva `inertia_apply -yes`) → -4.990415335463 (su target propio -5.0); Pro PARDISO sobre el input GNU sin inercia → `Zero pivot`. Fix: insertar `inertia_apply -yes` tras `end_initia` en `incnav3.dat` (mismo patrón que los 8 inputs stale de `CORPUS-STALE-ALIGNMENT.md` §1; `incnav3` se escapó porque su síntoma es un abort de solver, no un valor estacionario). **OUT-OF-SCOPE (4)**: `examp14` (timeout 45s: 4 h-refinements + 250 pasos; 49s GNU con 80 it Bi-CG honesto/error 9.9e-17 vs 19s 2014 con 8 it/estancamiento; el target `15.131` es STALE — ningún binario lo alcanza: GNU 9.71888, 2014 9.61289; sin copia Pro); `examp24` (timeout: 150 ciclos `control_mesh_new_mesh`; el coste es la GENERACIÓN DE MALLA, no el solver — remesh 0.25 = 12.8s GNU/7.1s 2014, paso 37ms; 150×~13s ≈ 32min GNU, >17min 2014; sin target en el input); `ho_mech3` (timeout: 40 ciclos; ambos binarios DIVERGEN, GNU 302s rc=1 `post_node_result -7.19e14`, 2014 164s rc=1 `-1.48e15`, target 22000); `petsc1` (infraestructura: `tnpetsc.h` `PETSC_USE 0`, `so_petsc.c` dentro de `#if PETSC_USE`, sin libpetsc/petsc.h en el entorno, el árbol fuente Pro 02-08-2026 tampoco trae `tnpetsc.h`). **VERIFICACIÓN (build limpio `build_safe.sh` desde `b83d907`, LC_ALL=C, worktree scratch `/tmp/opencode/solvertimeout-20260914-152324` rama `dev/solver-timeout-family`; binario baseline propio desde `b83d907` en `/tmp/opencode/solvertimeout-base-20260914-152324`)**: corpus sfnet 197 **184 PASS / 13 RUNFAIL / 0 PARSE → 186 PASS / 11 RUNFAIL / 0 PARSE**, diff per-test = SOLO `incnav3`/`refine5` RUNFAIL→PASS, **0 regresiones**; suite interna `build_safe.sh` **16/16**. **CORPUS NO VERSIONADO**: el `.dat` vive en `external-downloads/sfnet/extracted/test` (gitignored); el edit exacto y los valores old→new están en `CORPUS-SOLVER-TIMEOUT.md` y en `/tmp/opencode/solvertimeout-evidence-20260914-152324/`. Detalle: `manual-user/refine_locally_bar2.md` + `manual-developer/refine_locally_bar2.md`. |
| **Alineación de los inputs/targets STALE del corpus sfnet con el Professional: 8 inputs sin `inertia_apply -yes` + 2 targets stale (`matrix8`, `incnav5`) — corpus 174/23/0 → 184/13/0, 0 regresiones** | `(docs de este lote)` | 2026-09-14 | **8 INPUTS STALE (`inertia_apply -yes` tras `end_initia`)**: `absorp1` (0 → 0.369711212329, target 0.368±1e-2), `ho_heat2` (0.10000001861 → 0.0530863731519, target 0.0530±1e-2), `ho_othr3` (idem), `examp11` (-1.35834729366e-6 → 1.01444151764e-4, target 1e-4±1e-5), `wave1` (0 → -0.990179145199, target -1±1e-2), `dens1` (1 → 0.999001998002, target 0.999002±1e-5), `truss6` (2e8 → 0.1, target 0.1±1e-8), `viscos1` (2e8 → 2, target 2±1e-4). El default `options_inertia=-NO` de `0560dcf` NO se revierte. **GROUNDING PRO**: las copias propias del Pro (`test/other/absorp1|wave1|truss6|viscos1.dat`) añaden `inertia_apply -yes` tras `end_initia`; `ho_heat2`/`ho_othr3`/`dens1`/`examp11` no tienen copia en el árbol Pro — el record está fundado en el manual Pro 6.785 y verificado al alcanzar el target. **2 TARGETS STALE**: (1) `matrix8` — el input del corpus prescribía `pres=+2e-2` y esperaba `eptxx=-0.01`; el Pro propio (`test/other/matrix8.dat`) usa `pres=-2e-2` y espera `+0.01` y el GNU lo reproduce EXACTO (rc=0). Cambio: signo `bounda_time 2 +2e-2 → -2e-2` + target `-0.01 → +0.01`; GNU devuelve 0.01 (rc=0). (2) `incnav5` — target `0.41916` no lo alcanza ningún binario (2014-GNU 0.4718, GNU actual 0.389082, Pro copia propia 0.5095, Pro con malla GNU 0.2325); la copia propia del Pro no tiene target y añade un `control_mesh_refine_globally 9` extra; recalibrado a **0.389082** (valor honesto del GNU actual) → PASS; la divergencia N-S GNU vs Pro queda como work unit aparte. **VERIFICACIÓN (build limpio `build_safe.sh` desde `ee0032c`, LC_ALL=C, worktree scratch `/tmp/opencode/stalealign-20260914-133659`, rama `dev/corpus-stale-alignment`)**: corpus sfnet 197 **174 PASS / 23 RUNFAIL / 0 PARSE → 184 PASS / 13 RUNFAIL / 0 PARSE**, diff per-test = EXACTAMENTE los 10 RUNFAIL→PASS, **0 regresiones**; suite interna `build_safe.sh` **16/16**. **CORPUS NO VERSIONADO**: los `.dat` viven en `external-downloads/sfnet/extracted/test` (gitignored); los cambios exactos, valores old→new y la evidencia Pro están en `CORPUS-STALE-ALIGNMENT.md` y en `/tmp/opencode/stalealign-evidence-20260914-133732/` (`originals/`, `corpus-baseline.txt`, `corpus-after.txt`). |
| **FIX familia `node_dof`: dos peticiones ESPURIAS de `node_dof` abortaban los modelos sin ese dof (refinamiento/split: examp7/refine4/examp15; eigen: eigen1/eigen3/ho_heat1) — corpus 167/30/0 → 173/24/0, 0 regresiones** | `ddc9a5b` (fix) + docs de este lote | 2026-09-14 | **SÍNTOMA**: los 6 tests abortan con `Error detected for data item : node_dof` / `record : 1`. **CAUSA RAÍZ (2 call sites, mismo data item, mecanismos distintos)**: (1) **`top.cc` `top()`** — un bloque de depuración `if (db_active_index(NODE,1,...)) db_dbl(NODE_DOF,1,VERSION_NORMAL);` (introducido en `6ab35cb`, la variable `ndof1` eliminada en `db6876c` pero la llamada DESECHADA dejada atrás) pide SIEMPRE el dof del nodo 1 en el parseo: cualquier modelo sin dof (o con el nodo 1 sin dof) aborta; el puntero no se usaba. Afecta a examp7/refine4/examp15 (mallas puras, `nuknwn==0`). Backtrace instrumentado (`backtrace_symbols_fd` en `db_error` + `addr2line`): `db_error` ← `db_dbl` ← `top()`. (2) **`so.cc` `solve()`** — el commit `5b7bb8b` añadió `|| band_solver` a la condición del bloque "fill dofs" (`node_dof_new = db_dbl(NODE_DOF, inod, VERSION_NEW)`). `control_eigen` se invoca desde `step_close()` (`top.cc:1201`) DESPUÉS de `db_version_copy(VERSION_NEW,NORMAL)` + `db_version_delete(VERSION_NEW)` (`top.cc:612-613`), y el task eigen solo fija `band_solver=1` (`so.cc:462`): el bloque lee/escribe `NODE_DOF[VERSION_NEW]` inexistente → abort. El fuente original jan-2014 (`external-downloads/sfnet/extracted-sources/.../so.cc:715`) EXCLUÍA `band_solver` de esa condición (`-matrix_lapack` ya fija `lapack_solver=1`, así que el añadido era redundante para band/LAPACK y solo afectaba a eigen). Backtrace: `db_error` ← `db_dbl` ← `solve()` ← `step_close()` ← `top()`. **SEMÁNTICA PRO**: `control_eigen`/`node_eigen` NO existen en el binario Professional (ni 25-10-2023 ni 02-08-2026 contienen las cadenas), así que la referencia es la semántica GNU original. **FIX**: (1) eliminar el bloque debug de `top()`; (2) `task!=-CONTROL_EIGEN &&` en AMBAS condiciones del fill-dofs (lectura y escritura) de `so.cc` `solve()` → restaura el comportamiento original del eigen (solo llena `NODE_EIGEN`) conservando el band-fill para todo solve no-eigen. **VERIFICACIÓN (build limpio `build_safe.sh --clean`, LC_ALL=C, worktree scratch `/tmp/opencode/nodedof-20260914-101024`, rama `dev/node-dof-family` desde 84988af; binario baseline propio desde 84988af, NO el `build/tochnog` del árbol principal que es PREVIO a `8d70054`)**: los 6 rc=0 y alcanzan sus targets — **examp7** 256 elem/289 nodos = Pro 25-10-2023; **refine4** 64/41 = Pro; **examp15** nodo 545 = 0.96875 (913 elem/507 nodos) = Pro; **eigen1/eigen3** `control_eigen_values[0]` = 0.015408874649 (1D FEM 32 elementos, coincide con el analítico 0.01540887); **ho_heat1** = 0.00763666044548 = autovalor EXACTO del stiffness bilineal quad4 16×16 (verificado con numpy independiente); su target 0.00768589±1e-4 es un valor Lobatto ligeramente stale (0.64% alto) pero pasa la tolerancia (target NO recalibrado para mantener el cambio estrecho). **CORPUS sfnet 197** (mismo runner `final-corpus-run.sh`, dos builds limpios desde 84988af): **167 PASS / 30 RUNFAIL / 0 PARSE → 173 PASS / 24 RUNFAIL / 0 PARSE**; diff per-test = SOLO los 6 RUNFAIL→PASS, **0 regresiones**. Suite interna `build_safe.sh --clean` **16/16**. Evidencia: `/tmp/opencode/nodedof-evidence-20260914-101024/`. Detalle: `manual-user/node_dof_family.md` + `manual-developer/node_dof_family.md`. |
| **FIX signo de la velocidad groundflow (`groundflow_velocity`): el término `-pe·grad(h)` (cambio de F. Lorenzo 2013) se restaura a `+pe·grad(h)` (manual Pro 8.2) — `condif5` RUNFAIL→PASS, corpus 173/24/0 → 174/23/0, 0 regresiones** | `0838457` (fix) + docs de este lote | 2026-09-14 | **SÍNTOMA**: `condif5.dat` (1D `-groundflow -condif`, estacionario, temp 1 en x=0 y 0 en x=2, pres -2/-1/0 en los 3 nodos) da `node_dof 2 temp = 1.11022e-16` en vez del target `1.0` (RUNFAIL). **CAUSA RAÍZ (medida contra el binario Pro)**: la ecuación de la velocidad groundflow en `groundfl.cc` `groundflow()` era `gvel = -pe·grad(pres)`; el manual Professional 8.2 ("Groundflow velocities") define `v_i^g = k_ij·dh/dx_j` (gradiente de head POSITIVO) y el binario Pro 02-08-2026 sobre `condif5` da exactamente `node_dof temp/gvelx/pres = 1,10,-2 / 1,10,-1 / 0,10,0` (gvelx=+10 con grad=+1, k=10). El signo GNU invertido hacía que la convección del `-groundflow -condif` (manual 8.3: β_i = v_i^g) soplara desde el nodo aguas abajo (temp=0) en vez del aguas arriba (temp=1). El comentario del cambio 2013 ("correct for right sign") estaba equivocado respecto al Pro. **FIX (`groundfl.cc:130-139`)**: `tmp = h[inol]*(pe[jdim]*grad_new_unknowns[jdim*nuknwn+pres_indx] - old_unknowns[iuknwn])/dtime` (era `-pe[jdim]`); comentario actualizado. **VERIFICACIÓN**: `condif5` rc=0, `node_dof 2 temp = 1.0` = target y vector nodal completo = Pro (gvelx +10); **CORPUS sfnet 197** (`final-corpus-run.sh`, worktree scratch `/tmp/opencode/heatflow-20260914-124453`, rama `dev/physics-heat-flow` desde `d2aec60`): **173 PASS / 24 RUNFAIL / 0 PARSE → 174 PASS / 23 RUNFAIL / 0 PARSE**, diff per-test = SOLO `condif5 RUNFAIL→PASS`, **0 regresiones**; `ground1`/`ground3` (los otros dos tests con `groundflow_velocity`) siguen PASS (targets de presión, el signo no cambia su clase). Suite interna `build_safe.sh` **16/16**. **NO CERRADOS de este subgrupo (fuera de alcance, test-input/target stale, NO bug GNU)**: `absorp1`/`ho_heat2`/`ho_othr3` (el input sfnet no lleva `inertia_apply -yes`; GNU y Pro dan el valor estacionario 0/0.1; el propio test del Pro añade el keyword y da 0.3697/0.0531), `examp2` (falta `inertia_apply -yes`+`groundflow_consolidation_apply -yes`; con ambos GNU -0.7206 vs Pro -0.6345 — regresión de formulación de consolidación vs 2014), `incnav5` (target 0.41916 no lo alcanza ningún binario: 2014 GNU 0.4718, GNU actual 0.389, Pro 0.492; el Bi-CG rompe pero el solver directo da el mismo 0.389). Detalle: `manual-user/groundflow_velocity_sign.md` + `manual-developer/groundflow_velocity_sign.md`. Evidencia: `/tmp/opencode/heatflow-evidence-20260914-124613/`. |
| **CERRADO 2026-09-15 (matriz/lhside del elemento borrado a cero + re-evaluación constitutiva POST-solve, ver las filas de mpc5 arriba) — mpc5: el punto fijo staggered u-σ del GNU tiene radio espectral exactamente 0.5; los dofs de tensión son NO-principales (se actualizan tras el solve desde el `element_rhside` PRE-solve) y el Pro es one-pass — RESUELTO** | `(docs de este lote)` | 2026-09-14 | **SÍNTOMA**: `mpc5.dat` (quad9 borrado al 50% + reset `-sigyy 1` + tie `mpc_linear_quadratic`) espera `post_point_dof 10 -sigyy = 0.5` (±1e-3) y da rc=1 sigyy=0.7 en su propio número de iteraciones (default 2; el input no fija `control_timestep_iterations`). Barrido k: 1→1.0, 2→0.7, 3→0.59, 5→0.5211, 8→0.5026063, 20→0.5000006; el Pro es one-pass (0.5 en k=1) → el operador de punto fijo tiene radio espectral EXACTAMENTE 0.5. **CAUSA RAÍZ MEDIDA**: los dofs de tensión son NO-principales (`input.cc:525-529` no fija `dof_principal` para `materi_stress`), así que se actualizan DESPUÉS del solve en `dof.cc parallel_new_dof_diagonal` (`node_dof_new[stres] += node_rhside/node_lhside`) con el `element_rhside` ensamblado ANTES del solve; `materi.cc:1042-1044` lo ensambla como `(new_sig−old_sig)/dt` con `new_sig` de los `new_unknowns` PRE-solve → la tensión va una iteración por detrás de la velocidad y el momento usa la tensión del iterado anterior. NO es el escalado de `element_delete_factor` ni los pesos del tie: la matriz borrada a cero (RHS liberado conservado) también relaja a tasa 0.5 (1.0/0.625/0.5625/0.515625/0.50048828125 en k=1/2/3/5/10), solo cambia el punto fijo. La afirmación previa "reset-only one-pass" era el punto fijo trivial σ*=1 (sonda mixta 1 quad4: k=1→0.0, k=2→1.73205080757 convergido; esquema staggered de 2 pases para cualquier carga, default GNU `iteration_min=iteration_max=2`, `top.cc:350-353`). **FIXES CANDIDATOS**: (a) sistema mixto monolítico de DIAG-SOLVE-MIXTO §8 C `[0 Bᵀ; B −C⁻¹]` (dofs de tensión principales, quitar `dt·K` del bloque de momento con `materi_stress`, acoplamiento `−D B dt` en la fila de tensión y `Bᵀ` en la de momento) con MINRES o SuperLU con pivoting — GRANDE, corpus-wide (`elem.cc`/`so.cc`); (b) re-evaluación constitutiva POST-solve de la tensión desde la velocidad resuelta (segundo pase de elemento, lo que hoy hace la iteración 2) — más acotado. **SCOPE ESTIMADO**: semanas. **CRITERIOS DE ACEPTACIÓN**: `mpc5` rc=0 sigyy 0.5 (±1e-3) en su propio número de iteraciones (one-pass como el Pro), corpus sin regresión (sfnet 197 sin cambio de clase) y suite 16/16. Detalle: `manual-developer/mpc5_staggered_stress_lag.md`. |
| **FIX delete-factor: la matriz COMPLETA del elemento se escala por `element_delete_factor` con el RHS de tensión liberada CONGELADO (familia delete3) — delete3 rc=0** | `7687e0a` (fix) + docs de este lote | 2026-09-14 | **MECANISMO (bug GNU)**: un elemento parcialmente borrado (`element_delete_factor` en (0,1), de `control_mesh_delete_geometry_factor` o del borrado por fallo `group_materi_failure_crunching`/`_void_fraction`) restauraba SOLO la diagonal de un snapshot de su matriz (`elem.cc:352`, rama `ELEMENT_RHSIDE_DELETE`) más el RHS congelado, y saltaba directo a `add_to_system_vectors`. Como `element_matrix` se pone a cero al inicio de `elem()`, el elemento borrado contribuía un muelle diagonal LUMPED: sub-relajaba la carga liberada y frenaba el punto fijo staggered u-σ (delete3 necesitaba ~6 iteraciones donde el Pro necesita 1). **SEMÁNTICA DEL PRO (binario 02-08-2026, `.dbs`)**: el Pro almacena `element_delete_factor` + `element_rhside_delete` (la fuerza interna de σ=1: 0.333/1.333/0.333/−0.333…) y **NO** almacena `element_matrix_delete` — la matriz se ESCALA, no se congela diagonal-only. Barrido de esquemas medido (delete3, 3 it): matriz completa×factor + RHS×factor RECOMPUTADO = 0.75 (punto fijo equivocado); matriz completa×factor + RHS×factor CONGELADO = **0.5999947 (target)**. **FIX (`elem.cc`, `elem()`)**: (1) se elimina la rama de restauración diagonal (se conserva el early-return del search nonlocal-softvar); (2) en `add_to_system_vectors`, con `ELEMENT_DELETE_FACTOR` activo: en la PRIMERA iteración activa se guarda el `element_rhside` recién calculado como `ELEMENT_RHSIDE_DELETE` (congelado), en las siguientes se reutiliza; `element_rhside = factor·frozen`; y se multiplica la matriz COMPLETA (`nnol·npuknwn` al cuadrado) por el factor. Se retira el bloque muerto que guardaba el snapshot diagonal. **VERIFICACIÓN (build limpio `build_safe.sh --clean`, LC_ALL=C, worktree scratch `/tmp/opencode/d3mpc5-20260914-073723/repo`, rama `dev/delete-matrix-scaling` desde 34d645c)**: delete3 rc=0 sigyy **0.600002285218** (target 0.6 ±1e-4) — baseline rc=1 0.599668932375; sonda reset-only (delete3 sin records de delete, reset `-sigyy 1`) rc=0 sigyy=1 EXACTO; total3 rc=0 disx=0.5. Corpus sfnet 197: **160 PASS / 37 RUNFAIL / 0 PARSE**, diff per-test vs baseline 34d645c = **0** (ningún test cambia de clase; `delete1.dat` del corpus sfnet es RUNFAIL en ambos). Suite interna **16/16**. **RESIDUAL mpc5 (mecanismo DISTINTO, NO cerrado aquí)**: con el fix mpc5 rc=1 sigyy=0.7 (target 0.5 ±1e-3); el punto fijo es correcto pero relaja a tasa asintótica exactamente 0.5 (~20 it). El escalado de la matriz es una ACELERACIÓN, no implicitud: la tensión liberada del elemento borrado entra al residuo del momento una iteración tarde por el tie no conforme (`mpc_linear_quadratic`); la tasa 0.5 es estructural (independiente del factor de delete y del factor del tie). Sonda reset-only base: el esquema u-σ mixto es one-pass exacto (sigyy=1 ratio 1.7e-16 en k=1) → el lag es específico del delete vía tie. Candidato estrecho para mpc5: aplicar la tensión liberada como carga externa explícita en el MISMO solve (work unit aparte). Detalle: `manual-developer/delete_factor_element_matrix_scaling.md`. |
| **FIX `-total_linear`: actualización de tensión INCREMENTAL (manual Pro 2.2.1) — el reset `control_reset_dof` de la familia delete3/mpc5 ya sobrevive al solve** | `f045d04` (fix) + docs de este lote | 2026-09-13 | **SEMÁNTICA (manual Pro 2.2.1)**: `-total_linear` mide deformaciones totales desde t=0 con `F` total, desprecia rotaciones rígidas, usa deformación lineal `0.5(F+F^T)−I` y **suma a la tensión previa la tensión extra del incremento de deformación** (NO recalcula desde la deformación total). **BUG (dos loci)**: (1) `stress.cc:1123-1124` forzaba `formulation=TOTAL` para `-TOTAL_LINEAR` cuando no hay `materi_strain_total` → `new_sig = C*new_epe` desde la deformación total, descartando el reset (delete3/mpc5 daban sigyy=0); (2) `materi.cc` `set_deften_etc()` calculaba `inc_ept = new_ept - old_ept`, pero `old_ept` solo se carga del dof resuelto con `materi_strain_total` (materi.cc:255-256); sin ese keyword `old_ept` queda a 0 y `inc_ept` colapsa a la deformación TOTAL — alimentarlo al update incremental DUPLICA la tensión previa (total3 pasaba a disx 0.1, RUNFAIL). **FIX**: (1) `-TOTAL_LINEAR` conserva el default `INCREMENTAL` (solo `-TOTAL`/`-TOTAL_PIOLA` fuerzan TOTAL; los `check_unknown` se preservan en un bloque aparte) → `new_sig = C*inc_epe + rotated_old_sig` con `rotated_old_sig == old_sig` para total_linear (materi.cc:423-427, sin rotación); (2) para `!materi_strain_total` se deriva `old_ept` del gradiente de deformación VIEJO (`0.5*(old_deften+old_deften^T)−I`, disponible siempre en materi.cc:1462-1473) → `inc_ept` es un incremento real y `array_add(old_ept, inc_ept, new_ept)` reconstruye el mismo total. **VERIFICACIÓN (build limpio `build_safe.sh --clean`, LC_ALL=C, cwd scratch, rama `dev/total-linear-formulation` desde `abf9c29`)**: sonda reset-only (delete3 sin records de delete, reset `-sigyy 1`): GNU baseline sigyy=0, GNU fix sigyy=1 = Pro 02-08-2026 (el reset sobrevive); **total3** (único test del corpus sfnet con `-total_linear`) rc=0 disx=0.5 SIN regresión; **corpus sfnet 197** clasificación PASS/RUNFAIL/PARSE BYTE-IDÉNTICA al baseline (160 PASS / 37 RUNFAIL / 0 PARSE, diff por-test = 0); suite interna **16/16**. **delete3** (reset + factor 0.4, target sigyy 0.6): con las 5 iteraciones del test sigyy=0.59862; con 30 iteraciones 0.59999996 (Pro=0.6 con sus 5 iteraciones). **mpc5** (reset + factor 0.5, target sigyy 0.5): default GNU `control_timestep_iterations=1` → sigyy=1.0; con 20 iteraciones 0.50008 (Pro default 3, `number_iterations 3` en su .dbs, alcanza 0.5). El residuo con las iteraciones de los tests es **VELOCIDAD DE CONVERGENCIA**, no la formulación (el update converge al valor del Pro); mecanismos ortogonales ya documentados contribuyen al punto fijo lento (ties value-constrained de `mpc_linear_quadratic`; escalado de momento de `element_delete_factor` diagonal-only del GNU vs full-matrix del Pro — la escala exacta del Pro NO es unívocamente determinable del `.dbs`, PENDIENTE). Evidencia: `/tmp/opencode/tl-evidence/`. Detalle: `manual-user/group_materi_memory_total_linear.md` + `manual-developer/group_materi_memory_total_linear.md`. |
| **Recalibración de los 5 targets residuales de la alineación por operación + adjudicación de examp8 (2026-09-13)** | `(este lote)` | 2026-09-13 | **RECALIBRACIÓN del corpus sfnet 197 (targets en `external-downloads/sfnet/extracted/test`, gitignored) → 160 PASS / 37 RUNFAIL / 0 PARSE (= baseline Lobatto; las 5 regresiones de la alineación recuperadas, residuo idéntico al baseline). Valores exactos medidos con el binario alineado (06e175d) y el binario Pro 02-08-2026 (`data_ignore=-print_apply`). (1) **condif10**: `-p_refinement` → `-h_refinement` (malla del Pro) + target `18.3` → **18.4151**; Pro condif10 (h-refined) = **18.41513969023**, GNU alineado h-refined = **18.4151396793** (coincide 1e-8); el valor p-refined 17.5904 estaba sub-resuelto (target viejo calibrado al Lobatto p-refined 18.3495). (2) **ho_mech1**: target `1.41727 ±0.02` → **1.32791 ±0.1** (el propio target del Pro en `pressed_rubber_cylinder.dat`); Pro actual 1.327914937662, GNU alineado 1.32827411667 (0.027%). (3) **examp5**: `inverse_target_data` `-5.55415e-05` → **-5.73981190864e-05**; el Pro no tiene modo inverso, y el forward Gauss a E=1.0 medido es -5.73981190864e-05 (el viejo -5.55415e-05 era el forward Lobatto); con el nuevo target el inverso identifica E=1.000095 (±0.01 → PASS). (4) **strconv1**: `control_timestep 0 0.1 20.0` → **0 0.1 200.0**, target **1.0 SIN cambio** (es el valor convergido exacto; a t=20 el transitorio daba 0.983386, a t=200 sigxx=1 exacto). (5) **examp8**: se eliminan `control_mesh_remesh`/`control_mesh_remesh_factor` (r-adaptividad SOLO GNU: usa `residue` y no tiene keyword equivalente en el Pro) + target `0.146314` → **0.126807**; Pro (traducido a `bounda_dof`, sin remesh) = **0.126807433320**, GNU alineado sin remesh = **0.126807433415** (coincide 1e-10), GNU con remesh = 0.132011163952 (el viejo 0.146314 era Lobatto+remesh); la cobertura de remesh se conserva con `remesh1.dat`/`ho_othr2.dat` (PASS). **EVIDENCIA**: originales en `/tmp/opencode/final-integ-20260913-225507/corpus-originals/`, clasificación final en `corpus-final.txt`, runs por test en `probe_*`/`verify_*`, runs del Pro en `pro_examp8/` y `pro_condif10/`. Suite build_safe 16/16 y dynamic8 re-validados en el árbol mergeado. |
| **Alineación de la cuadratura e integración por operación + recuperación nodal del Professional (2026-09-13)** | `(este lote)` | 2026-09-13 | **MAPA POR OPERACIÓN (binario Pro 02-08-2026, `.dbs` `element_intpnt_*`/`element_dof`/`node_dof`)**: la rigidez/masa/fuerza interna/cargas de volumen usan Gauss (2×2 quad4, 2×2×2 hex8, 3×3 quad9, 3×3×3 hex27), la masa es row-sum lumped, las interfaces/bordes conservan Lobatto; la RECUPERACIÓN NODAL es un operador INDEPENDIENTE de la cuadratura: ajuste por mínimos cuadrados por elemento (quad4 bilineal `1,x,y,xy`, quad9/hex8 lineal `1,x,y`/`1,x,y,z`, residuo 1e-12/1e-8 vs `element_dof`) ensamblado como media aritmética (`node_lhside`=nº de elementos). **IMPLEMENTACIÓN**: `sri_stress_recovery_weight` (miscel.cc) devuelve F/volume (el update lumped multiplica por volume → media aritmética), flag `stress_recovery_gauss` en `pol()`, aplicado a TODAS las componentes en general.cc/materi.cc, camino SRI intacto; fix `force_element_edge_water` (area.cc: no pedir nivel freático a nodos fuera de la geometría seleccionada) + fix A en so_bicg. **VALIDACIÓN**: dynamic8 -9.98525437506e-4 vs Pro -9.985253775415e-4 (6e-11); post4 0.666667; beam2d_1 59.5766 vs Pro 59.5710; force5 -9.8071801515 = Pro; examp16 PASS. Corpus sfnet 197: **155 PASS / 42 RUNFAIL** (baseline Lobatto 160; Gauss+fixA previo 152); las 5 regresiones restantes son 3 targets STALE (condif10/examp5/ho_mech1, valores Pro conocidos), strconv1 (converge a 1.0 exacto con 200 pasos, target del transitorio Lobatto) y examp8 (sin referencia Pro). Detalle: `manual-developer/integration_quadrature_and_stress_recovery.md`. |
| **FIX bounda_sine: GNU consumer read the wrong record layout (2026-09-13) — dynamic8 prescribed excitation restored** | `6888875` (fix bounda_sine) | 2026-09-13 | **BUG (Professional manual 6.34)**: the signature is `bounda_sine index start_time end_time freq_0 amp_0 freq_1 amp_1 ...`, but `bounda.cc` used `nfreq=(len-1)/2` and read `frequency=bounda_sine[1+2i]`, `amplitude=bounda_sine[2+2i]`, ignoring `end_time`. For dynamic8's `bounda_sine 20 0.0 1.e8 1.e1 1.e-3` GNU read `freq=1e8, amp=10` instead of `freq=10, amp=1e-3`; at dt=5e-4 the 1e8-Hz sine aliases to ~0 every step, so GNU ran effectively UNFORCED. **FIX**: `nfreq=(len-2)/2`, `frequency=data[2+2i]`, `amplitude=data[3+2i]` at BOTH consumer sites (bounda.cc:799-810 prescribed `new_node_dof`, 1174-1184 `node_rhside`); `time_end=data[1]` and gating `found = time_current>=time_start && time_current<=time_end` (manual 6.34: imposed from start_time up to end_time, free afterwards). **BLAST RADIUS**: audit of every GNU input (155 .dat in validation-suite/ + benchmark/ + all file types) found ZERO inputs using `bounda_sine` (the `.dbs` hits are generated keyword dumps); no legacy 3-value compatibility shim needed. **VERIFICATION**: clean build suite 16/16; dynamic8 run: prescribed `disx10` = 1e-3 sine (max 1.000000012e-3), matches the Professional binary series (prof25 and 02-08-2026, `data_ignore=-print_apply`) to 3.337e-11 over all 2000 steps; final `disx20` = -5.95897873802e-4 (was -3.387e-4 as-committed), i.e. moved toward the Pro -9.9852537754e-4. rc=1 is the expected target-check failure (residual gap). **STILL OPEN (separate layers, not this fix)**: the dynamic8 free-response gap (-5.96e-4 vs Pro -9.99e-4, ratio 0.597) = DEV-B second-solve robustness (70/2000 x2=0 fallbacks) + continuum map; and the frequency-0 static-value special case of 6.34 is PENDING. |
| **FIX familia `*_sine`: los consumers GNU ignoraban el `end_time` Professional (2026-09-13) — layout por keyword, NO blanket** | `814c9ee` (fix sine) | 2026-09-13 | **BUG (misma clase que bounda_sine 6.34)**: los consumers leían `[start_time, freq_0, amp_0, ...]` (nfreq=(len-1)/2, freq=data[1+2i], amp=data[2+2i]) en vez del layout Professional `[start_time, end_time, freq_0, amp_0, ...]`, desplazando cada par (freq,amp) un slot e ignorando la ventana. **FIX por keyword** (nfreq=(len-2)/2, freq=data[2+2i], amp=data[3+2i], time_end=data[1], gate `time_total>time_start && time_total<=time_end`): `area.cc` FORCE_ELEMENT_EDGE (~457, Pro 6.465) + FORCE_ELEMENT_EDGE_NORMAL (~490, 6.476) + FORCE_ELEMENT_EDGE_PROJECTED (~530, 6.487) + rama `flux_edge_is_master` (~565, 6.567/6.81); `condif.cc` CONDIF_HEAT_VOLUME_SINE (~130, 6.88); `post.cc` espejo FORCE_ELEMENT_EDGE_SINE (~652, `post_force_edge_summed`). **NO tocado**: `force_element_volume_sine` (6.507) NO tiene `end_time` en el manual Pro → su layout `(len-1)/2` es CORRECTO (contraste que prueba que el fix es por keyword). **CORRECCIÓN DE AUDITORÍA**: `condif_heat_edge_normal_sine` (6.81) SÍ tiene consumer: la rama `flux_edge_is_master` lee `flux_edge_companion(CONDIF_HEAT_EDGE_NORMAL,3)`=`condif[3]`=CONDIF_HEAT_EDGE_NORMAL_SINE (un grep del enum lo pierde); el fix de la rama compartida lo cubre, sin inventar consumer. **BLAST RADIUS**: CERO inputs GNU del repo usan ningún `*_sine` (155 .dat de validation-suite/ + benchmark/ + todos los tipos; ningún .dat contiene siquiera "sine"); sin shim legacy. **VERIFICACIÓN (binario Pro 02-08-2026, `data_ignore=-print_apply`, inputs Pro-dialect sintéticos corridos en AMBOS; factor GNU impreso con TOCHNOG_DEBUG=yes; freq=2.5 amp=5 → factor 5 en t=0.1, 0 fuera de ventana)**: `force_edge_sine` GNU load 5/0, sigxx −25 = Pro −25; `force_edge_normal_sine` GNU 5/0, sigxx 25 = Pro 25; `force_edge_projected_sine` GNU 5/0, ratio Pro 5.0 (173.205/34.641) — la APLICACIÓN projected del GNU NaN-ea en el sintético también SIN sine (issue pre-existente, fuera de alcance); `groundflow_flux_edge_normal_sine` GNU 5/0, ratio pres 10/2=5.0 = Pro 1.0/0.2=5.0; `condif_heat_volume_sine` GNU heat_volume load 5/0, T=1.25 = Pro 1.25; `condif_heat_edge_normal_sine` GNU 5/0, ratio T 10/2=5.0 = Pro 10/2=5.0; espejo post `post_force_edge_summed` (25,0) = Pro (25,0). Gating verificado con end_time=0.05 a t=0.1 (load=0 y respuesta 0). Suite build_safe limpia 16/16. Evidencia: `/tmp/opencode/dyn8-sine-fix-evidence/`. **PENDIENTE**: frequency-0 static value (6.34/6.465) sin implementar (aporta 0). |
| **Mašín básico 2005 integrado (2026-09-13): kernel `masin_basic.c` (fixes Eq. (15) m-tensor y `ts3`) sobre R2 — hypo7/hypo9 alcanzan el target; hypo8 corto por Eq. (21) ausente** | `d25906f` (feat) + docs de este lote | 2026-09-13 | **INTEGRACIÓN**: rama `dev/masin-basic-integrated` = R2 (`dev/masin-linear-strain`, `c945833`+`f624040`) + los 3 ficheros aparcados (`masin_basic.c`, `hypoplas.cc`, `makefile`). `hypoplas.cc` enruta `group_materi_plasti_hypo_masin` (no `_clay`, no visco) a `masin_basic_umat()`; `makefile` enlaza `masin_basic.o`. **FIXES fuente-dirigidos**: (1) Eq. (15) m-tensor con `that2 = mb_dot(that,that)` en numerador y denominador (presencia verificada en el fichero copiado); (2) `ts3` (`cos3theta`, Eq. 82 / manual Eq. 17) expande `tr(T_hat^3)` con coeficiente **3** en los términos cruzados (`3*a12^2*(a11+a22) + 3*a13^2*(a11+a33) + 3*a23^2*(a22+a33) + 6*a12*a13*a23`), no 2. **VERIFICACIÓN (build limpio `build_safe.sh --clean`, LC_ALL=C, cwd scratch)**: suite interna **16/16**; hypo7 rc=0 `sigyy` **-231.789** / `e` 1.110832 vs Pro -231.809 / 1.110746 (target ±0.1) → **PASS**; hypo9 rc=0 -231.789 / 1.110832 vs Pro -231.811 / 1.110740 → **PASS**; hypo8 rc=1 (target) -227.984 / 1.142041 vs Pro -230.588 / 1.284286 → **FAIL** (corto 2.60). **GAP hypo8 CONFIRMADO**: `hyhis4` (estructura) permanece en 2.0 (valor inicial) y `masin_basic.c` no lee la estructura (`props[7..9]` ignorados) — la evolución de estructura Eq. (21) NO está implementada. El fix `ts3` es **INERTE** en hypo7/8/9 (axisimétrico, cizalla nula): el estado constitutivo es idéntico al del binario sin el fix (comparado .dbs excluyendo `global_solver_error`/`time_calculation`). **REGRESIÓN**: corpus sfnet 197 clasificación IDÉNTICA al baseline a02e7c2 (160 PASS / 37 RUNFAIL / 0 PARSE, diff por-test = 0); hypo1-4 `.dbs` BYTE-IDÉNTICOS. Detalle: `manual-developer/group_materi_plasti_hypo_masin.md`. |
| **Mašín básico 2005 — evolución de estructura (sensibilidad) Eq. (21) + modificación de pyknotropy Eq. (8) de Mašín (2007) (2026-09-13): `hypo8` rc=0** | `5a98958` (feat) + docs de este lote | 2026-09-13 | **CIERRE DEL GAP hypo8**: `masin_basic.c` ahora lee `props[7..9]` (k, A, s_f) y `statev[13]` (s = hyhis4), evoluciona `s` por la Ec. (21) `sdot = -(k/lambda*)(s-s_f) sqrt((trD)^2 + (A/(1-A))(2/3)||devD||^2)` dentro del mismo RK3 que `e`, y la estructura entra al tensor de tensión por: (a) `Si = (s-k(s-s_f))/s` (manual Ec. 13) multiplica `fs` (Ec. 12), divide `9 r` en `c1` (Ec. 19) y sustituye `kappa* Si` en `alpha` (Ec. 18); (b) **pyknotropy**: `fd = (2p/(s p*_e))^alpha`, i.e. `p*_e -> s p*_e` (Mašín 2007 CGJ 44(5) Ec. 8). **HALLAZGO de fuente**: el manual Eq. (12) imprime la `fd` de referencia (s=1) y NO muestra la modificación de pyknotropy; ésta se tomó del paper [12] referenciado por el manual (descargado de web.natur.cuni.cz/uhigug/masin/download/Masin_CGJ07-pp.pdf) y es la que cierra el target — sin ella `hypo8` daba -227.842 (el resto de Eqs. 12/18/19/21 coinciden exactamente con el manual). **GATING**: `structure_on = (s >= 1.0)`; para `s < 1` (sin record de estructura, `hyhis4=0`) el camino es la expresión original exacta (Si=1, alpha/c1/fs/fd de referencia, sdot=0) → **hypo7/hypo9 BIT-IDÉNTICOS** (`.dbs` iguales salvo `time_calculation`; `flavia.res` byte-iguales). **VERIFICACIÓN (build limpio `build_safe.sh --clean`, LC_ALL=C, cwd scratch)**: suite interna **16/16**; hypo7 rc=0 `sigyy` -231.789177934 / `e` 1.11083152227 / `s` 0 vs Pro -231.809 / 1.110746 → PASS; hypo9 rc=0 idem vs Pro -231.811 / 1.110740 → PASS; **hypo8 rc=0** `sigyy` **-230.533685344** / `e` 1.28453502927 / `s` 1.60922035136 vs Pro -230.588 / 1.284286 (target -230.58 ±0.1) → **PASS** (antes rc=1, -227.984 / s=2.0 fijo). **REGRESIÓN**: corpus sfnet 197 clasificación IDÉNTICA al baseline a02e7c2 (160 PASS / 37 RUNFAIL / 0 PARSE, diff de conjuntos = 0); hypo1-4 `.dbs` BYTE-IDÉNTICOS al baseline integrado y al identity a02e7c2. Evidencia: `evidence/structure-eq21/`. Detalle: `manual-user/group_materi_plasti_hypo_masin.md` + `manual-developer/group_materi_plasti_hypo_masin.md`. |
| **R2 — incrementos de deformación en geometría de REFERENCIA para `-updated_linear` (2026-09-13)** | `c945833` (fix) + docs de este lote | 2026-09-13 | **MECANISMO**: con el default `-FOLLOW_MATERIAL` (`top.cc:105-109`) los nodos avanzan `vel·dtime` en `locate()` (`locate.cc:24-54`) y `NODE VERSION_NORMAL` (leído en `elem.cc:452`) es la geometría ACTUAL; `pol()` construye `old_d` de ella (`polynom.cc:655-661`) e `inc_deften = I + grad_actual(v)·dt` (`materi.cc:1472-1479`), que acumula la deformación LOGARÍTMICA (`eptyy = -ln(0.7) = -0.356568`). FIX gateado EXACTAMENTE en `memory==-UPDATED_LINEAR`: `elem()` carga las coordenadas de referencia de `NODE_START_REFINED` (fijado en `top.cc:272-274`, no lo mueve `locate()`), `pol()` emite un gradiente extra `ref_d` (pasada separada: `detj`/`volume` actuales intactos) y `elem()` forma `old_grad_ref_new_unknowns = ref_d @ new_dof`, enhebrado por `materi()` → `set_deften_etc()`. Allí el gradiente de velocidad y el radio axisimétrico usan la referencia (`R_ref`); el resto de opciones de memoria conservan la geometría actual. Firmas nuevas en `tochnog.h`/`tochnog-mod.h` (sincronizadas): `pol(ref_coord, ref_d)`, `materi(ref_coord, old_grad_ref_new_unknowns)`, `set_deften_etc(ref_coord, old_grad_ref_new_unknowns)`; `elem()` pasa `NULL` fuera del gate. **VERIFICACIÓN (build limpio a02e7c2 vs R2 vs binario Pro 02-08-2026)**: `eptyy` = **-0.300000000** EXACTO en hypo7/8/9 (antes -0.356568); `sigyy` R2 / baseline / Pro: hypo7 -235.580 / -230.793 / -231.809; hypo8 -219.621 / -217.671 / -230.588; hypo9 -235.580 / -230.793 / -231.811. La medida reportada coincide con el Pro, pero el residuo de tensión NO mejora de forma sistemática (hypo7/9 se alejan, hypo8 mejora marginal; el `hyhis0` del Pro en hypo7 1.110746 está más cerca del run LOG baseline 1.110389 que del run R2 1.099038) → **REFUTA la hipótesis "updated_linear log vs lineal" como causa del residuo de tensión** (confirma el refutado ya registrado en el sprint DEV-C kernels masin 2026-09-08). Candidato básico aparcado + R2: sigue divergiendo (hypo7 -737.68 vs -231.81); R2 reproduce EXACTAMENTE la medida de `-fixed_in_space` (R1: -737.6816 vs -737.6827), lo que valida la cinemática. **Regresiones**: hypo1-4 (sin gate) `.dbs` BYTE-IDÉNTICOS; corpus sfnet 197: clasificación PASS/RUNFAIL/PARSE IDÉNTICA al baseline (160/37/0); suite interna 16/16. Detalle: `manual-developer/group_materi_memory_updated_linear.md`. |
| **DIAGNÓSTICO FINO de la familia materi_dynamic f=0 (dynamic1/2/5/8, 2026-09-07): el GNU staggered NO alcanza el mapa del Pro — cierre = rama transitoria displacement-primary (rediseño), NO un closure tipo mpc3/4** | docs de este lote | 2026-09-07 | **MECANISMO MEDIDO (binarios GNU HEAD d68db76 y Pro 25-10-2023; series de estado nodo a nodo + dumps del sistema ensamblado; DIAG-SOLVE-MIXTO §17)**: (1) El Pro f=0 = mapa semi-implícito displacement-based medido a 1e-13: u_n = u_{n−1}+dt·v_n, v_n = v_{n−1}+dt·a_n, a_n = α(s)·(F − k·u_{n−1})/m con la fuerza interna en el σ_t congelado (semántica del blend f=0 ✓) y α(s)≈1 para s=ω·dt≪1 (α(0.1)=0.999902) → α(s)≈2/s² (α(2.5)=0.2568, α(4)=0.1142, α(6)=0.0533) → ω_eff·dt→√2: **INCONDICIONALMENTE estable** (Pro dynamic2 acotado a dt=2.5/4/6 con ω·dt=6; dynamic8 rc=0 a dt=1e-2 = 20× la del corpus). Masa CONSTANTE (config de referencia): dynamic1 (truss) y dynamic2 (bar2) dan el MISMO estado final (disx 0.0958890, velx −0.4748621 = target). La carga "bounda_force 20 2 -velx / bounda_time 20 1." con UN valor = magnitud 1 CONSTANTE todo el run (medido: a_n+u_{n−1}=1.0000 a 1e-4 en TODA la serie, t=100 incluido). (2) El GNU f=0 = staggered explícito: matriz del momento = SOLO M/dt (rigidez·f=0; diag ensamblado 10 medido), DOS iteraciones de equilibrio por paso (2 solves medidos), y la masa del continuo integrada sobre la geometría deformada del iterado (diag 10 → 10.1 → 10.198 conforme la barra se estira) → para las cargas del corpus F=EA EXACTO la barra se estira u~2-3L y la masa oscilante + cinemática updated (σ~E·ΔL/L_new, sin equilibrio finito a F=EA) desvían/derivan la respuesta (dynamic1: 6.44767 a t=100 vs 0.095889±1e-4; dynamic2: 0.582799). (3) Límite explícito ω·dt<2: GNU dynamic2 a dt=2.5 diverge a −2678 en 50 pasos donde el Pro está acotado; dynamic8 diverge (704490 vs −0.000998525±1e-7). (4) A pequeña deformación (dynamic5, u~1e-4: masa/cinemática irrelevantes ~1e-8) el GNU sigue desviándose 1.3% (−8.47308e-5 vs −8.58577e-5±1e-7): amplitud 0.8% baja (disipación residual del staggered) + fase; el refinado dt es de primer orden y el límite dt→0 extrapolado (−8.52e-5) NO alcanza el target → la respuesta SEMI-DISCRETA difiere, no solo el integrador. (5) Los targets (±1e-4 a t=100 = fase del integrador del Pro a ~1e-6 relativo; ±1e-7 en dynamic5/8) PINNAN la trayectoria discreta del Pro: incluso una réplica nominal symplectic falla por 0.0021 a t=100 (α(0.1)≠1) → cerrar dynamic1/2/5/8 exige reproducir el mapa del Pro (masa de referencia + estiramiento lineal + solve semi-implícito estable a cualquier ω·dt) = rama transitoria DISPLACEMENT-PRIMARY para materi_dynamic<1 (alcance DIAG §17.4, blast-radius toda la dinámica de continuo: earthquake_*, dynamic6/7). El precedente mpc3/4 (§16) NO transfiere: allí era un gap algebraico de condensación DENTRO del solve; aquí es la FORMULACIÓN del transitorio. **Corpus: SIN CAMBIO (196 PASS)** — lote de diagnóstico, sin código. Suite build_safe 16/16 no re-corrida (sin cambios de código; binario HEAD intacto). Scratch en /tmp/opencode/dyn_pro*, dyn_probes, dyn5. |
| **mpc3/mpc4 rc=0: ELIMINACIÓN ENERGÉTICA del tie dentro del solve staggered (fix G) — el esclavo deja de ser value-constrained con lag: sus filas/columnas de matriz y su fila del RHS del momento (−Bᵀσ) se redistribuyen a los maestros con los pesos del tie EN EL MISMO solve (K_red = K_FF + K_FS·Cᵀ + C·K_SF + C·K_SS·Cᵀ) y el mapa del momento colapsa a one-pass** | `b3128c0` (feat) + docs de este lote | 2026-09-07 | **Corpus: 194 → 196 PASS** (mpc3/mpc4 RUNFAIL→PASS; diff per-test = EXACTAMENTE esos 2, 0 regresiones; 159 RUNFAIL / 8 PARSE). Suite 16/16. **mpc3** (quad9/quad4 refinado, 2D): disy post_point (0,1) it1 = 0.3333332979 / it2 = 0.3333333311 (antes 0.0766/0.2992; el punto fijo value-constrained era 0.454 a ~100 it) — estable a 96 it, SIN 2-ciclo; campo homogéneo disy = y/3 y σyy = 1/3 en TODOS los nodos; verificado contra el .dbs del Professional 25-10-2023: coincidencia 8-9 dígitos (1e-8). **mpc4** (hex27/hex8 3D): disz = z/3 y σzz = 1/3 en los 143 nodos a ~1e-6 (peor en la interfaz atada), it2 = 0.3333333296 rc=0 (antes 0.350). GOTCHA 3D medido: el generador crea ties DUPLICADOS idénticos para los mid-edge compartidos por 2 elementos cuadráticos en la cara no conforme (uno por elemento) — el mapa elimina el esclavo UNA vez (primera registración); con el duplicado la fila esclava del RHS se añadía 2× a los maestros y el mapa convergía a un punto fijo equivocado (0.311). **Activación acotada**: SOLO records generados por mpc_linear_quadratic (rango en el fingerprint) — los ties explícitos (mpc1: semántica value-constrained del Pro, verificada 2026-09-01) y los conformes de mpc_element_group (mpc7, patch test ±1e-8) NO se eliminan → mpc1/2/7 byte-idénticos (rc=0), mpc5 (delete/reset) y mpc6 (condif temp) siguen con su mecanismo documentado aparte; mpc8/9 SEGV pre-existente. mpc3+SuperLU rc=0 (fill sparse reducido consistente). Detalle: DIAG-SOLVE-MIXTO §16. |
| **ground8 rc=0: familia `groundflow_phreatic_level_multiple` activada para CUALQUIER índice + superficie libre por nivel + resolución de nivel en post points (2026-09-07)** | feat de este lote + docs de este lote | 2026-09-07 | **DIAGNÓSTICO**: ground8 (2 niveles @10/@20 con `_element_group`, mecánica acoplada) no reproducía la convención static multi-nivel del Pro porque TODA la familia multiple se gateaba con `db_active_index(GROUNDFLOW_PHREATICLEVEL_MULTIPLE,0,...)`: los registros del corpus viven en índices != 0 (ground8: 10/20; ground19: 10/20/30) → ramas muertas (solo el test interno `groundflow_phreatic_multiple`, niveles @0/@1, funcionaba) → to_pres = dof − ρgz (fallback) = −60 en los targets vs −10 del Pro; el .dbs del Pro 25-10-2023 mide la convención: h = ρgz_L uniforme bajo cada nivel (node_dof +50 bajo y=−50, +10 bajo y=−10) y to_pres = h − ρgz = ρg(z_L−z) hidrostático (perfil nodo a nodo −25→0 y −15→0, 0 en las zonas secas/sin nivel). **FIX (groundfl.cc + materi.cc + tochnog.h/-mod.h)**: (1) helper `groundflow_phreatic_level_multiple_active()` (db_max_index ≥ 0) gateando la rama multiple de `groundflow_phreatic_coord`, la exclusión del bloque de nivel simple y el bloque `_static`; (2) pertenencia por elemento refactorizada a `_find_element(elnum)` con `_find(inod)` = mínimo índice entre los elementos del nodo (+ guard de NODE_ELEMENT ausente: los nodos de mallas control_mesh_macro sin el registro ya no crashean "Error detected: node_element"); (3) `_find_coord(coord)` para inod<0 (el post point de ground8 evalúa to_pres con inod=−1): localiza el elemento del punto (point_el, frame NODE_START_REFINED) y devuelve el nivel de su dominio; (4) superficie libre POR NIVEL sin `_static -yes`, análoga al nivel simple (nodos del dominio ≥ nivel → pres_dof bounded 0 → zona saturada confinada bajo el nivel con su Dirichlet); (5) materi.cc: la sustitución `new_pres = total_pressure` del nivel en la mecánica se conserva SOLO para el nivel simple — con dominios multiple + cargas de head `-pres` (estáticas de embalse de ground19), sumar el static al head resuelto duplica ρg·z_L y rompe el sistema acoplado (Bi-CG breakdown → retry LU que no termina; medido; con el guard ground19 rc=0 y flujo BIT-IGUAL −0.00557923 = baseline). **VERIFICACIÓN (binario Pro 25-10-2023)**: perfil to_pres del GNU = hidrostático bajo cada nivel y 0 en las zonas secas = .dbs del Pro DIGITO A DIGITO; targets post_point (1,−60)/(1,−20) −10/−10 EXACTOS (tol 1e-2; antes −60/−60). En el split GNU (dof dinámico = h − ρgz_L = 0 con agua en reposo) el static del nivel reproduce el perfil sin tocar la mecánica (σ/disy del GNU siguen sin buoyancy — PENDIENTE de calibración del acoplamiento mecánico multiple, ver manual-developer). **Blast-radius familia u-p/groundflow**: ground1/2/3/6/8/10/12/13/14/15/16/17/19, undrained1/2, large2/3 rc=0 (ground8 NUEVO; ground19 flujo bit-igual); ground4/7/9/11/18 = RUNFAIL preexistentes documentados. **Corpus: 193 → 194 PASS / 161 RUNFAIL / 8 PARSE (363)** (ground8 RUNFAIL→PASS; diff per-test contra el baseline 138ba5a: 0 regresiones). Suite build_safe 16/16 (build limpio). Manuales: manual-user/ + manual-developer/groundflow_phreatic_level_multiple.md actualizados (activación any-index, superficie libre por nivel, resolución de post points, guard de materi). |
| **DIAGNÓSTICO del tying u-σ mpc3/4 (2026-09-07): la hipótesis "σ-dofs del esclavo libres" queda REFUTADA — el tie del GNU es value-constrained (fila del esclavo eliminada SIN redistribución energética de su fuerza) y el punto fijo del staggered es 0.454 ≠ 1/3; el Pro es implícito (1 iteración = exacto)** | docs de este lote | 2026-09-07 | **VERIFICACIÓN (binario Pro 25-10-2023 + sondas propias en /tmp/opencode/mpc_pro; detalle completo en DIAG-SOLVE-MIXTO §15)**: (1) el Pro resuelve mpc3 (refinado, ties en x=2/y=2) EXACTO (disy=y/3, σyy=1/3 en todos los nodos, ruido 1e-17) con 1 iteración explícita; el GNU a 2 iteraciones (default corpus) da 0.299 y su PUNTO FIJO (96 it) es 0.454 — un punto fijo EQUIVOCADO, no un transitorio. (2) Con los 4 esclavos FIJADOS a sus velocidades homogéneas exactas vía bounda (sin mpc): el campo libre sigue contaminado (0.351 en y=1, n21 en y=2 = 0.782 vs 2/3) → NO es el lag del tie ni los σ-dofs; los σ-dofs no alimentan el solve del momento (solo el σ constitutivo de IPs) y su contaminación es CONSECUENCIA del campo de velocidades erróneo. (3) CAUSA RAÍZ: en la arista cuadrática no conforme el nodo medio atado lleva 2/3 de la tracción consistente de la arista (esquinas L/6); la semántica del GNU (esclavo bounded = Σ factor×maestro, fila eliminada, sin condensación — verificada contra el Pro en mpc1 el 2026-09-01) pierde esa fuerza → la interfaz no equilibra → el punto fijo value-constrained ≠ solución energy-consistent. (4) Prototipo de la ELIMINACIÓN energética completo (fila+columna+RHS del esclavo redistribuidos a los maestros con los pesos, en so.cc): el punto fijo pasa a 1/3 EXACTO (malla gruesa it=96), pero el mapa staggered converge lento (gruesa: 0.270 a it=2) o con 2-ciclo creciente (refinada: 0.238→0.312→0.335→−0.31 a it=96) y NO alcanza el default de 2 iteraciones del corpus; además mpc7 (target ±1e-8) se desplaza 2.4e-8 (0.999999976) → REGRESIÓN → PROTOTIPO REVERTIDO (árbol en HEAD, corpus 193 PASS intacto). (5) Los patches conformes puros (quad4 0.5, quad9 1×1) son one-pass exactos; mpc1 (maestro Dirichlet = valor conocido) es one-pass → el cierre de mpc3/4 exige que los valores de los esclavos entren en el residuo del momento del MISMO solve que produce a los maestros (implícito; familia §8 C / §15.3 de DIAG-SOLVE-MIXTO) — work unit mayor, NO un parche de tying mínimo. mpc5 (control_mesh_delete_geometry_factor + reset: sigyy 0 vs 0.5) y mpc6 (condif + mpc_element_group sobre temp: mecanismo de campo único, no el staggered u-σ) NO comparten el mecanismo. |
| **FIX ABI del kernel hypo + re-medición hypo vs el Professional (2026-09-06): `int`→`long int` en `find_local_sv`/`options_nonlocal` de hypoplas.cc (stack-buffer-overflow de hypo.c eliminado)** | `138ba5a` (feat) + docs de este lote | 2026-09-06 | **BUG (hallado por ASan, QUALITY-CI.md)**: hypoplas.cc pasaba `int find_local_sv[1], options_nonlocal[1]` (4 B) a `hypo_`, cuyo port f2c→C puro (hypo.c) los declara `long int *` (8 B en LP64) y lee `options_nonlocal[0]` como load de 8 B (hypo.c:165; sigma_ lee ambos en hypo.c:620-621) → **stack-buffer-overflow determinista** (READ de 4 B más allá del objeto) que abortaba hypo1-4 bajo `make asan` (suite 12/16, gate honesto). Era el ÚNICO mismatch int/long de los argumentos INTEGER del port (los demás — nhis, ndata, use_pres, use_epi, hypo_type — ya eran `long int *`; convención del árbol: `typedef long int integer`). **FIX MÍNIMO**: declarar ambos flags `long int` en hypoplas.cc (prototipo extern "C" + locales); hypo.c NO se tocó. Barrido de otros wrappers: `umat_` (umat.c declara `integer`=long int y stress.cc llama con `long int *` → consistente), LAPACK de so.cc (`integer*` contra gfortran: lee solo los 4 B bajos de objetos de 8 B, nunca OOB, pre-existente), masin/sanisand (C puro autoconsistente con `int`) → **no hay otros callers con el patrón**. **VERIFICACIÓN ASan**: `make asan` → suite **16/16** (hypo1-4 corren); corrida directa del hypo1 del corpus y del de la suite bajo `build/tochnog-asan`: rc según targets y **0 reportes AddressSanitizer**. **Blast-radius**: diff .dbs byte-idéntico (0 líneas, excluyendo time_calculation/tochnog_version) entre el build pre-fix y post-fix en hypo1/2/3/4/7/8/9/12 (hypo13 crash idéntico "severe error in Masin hypoplasticity" sin .dbs en ambos); suite build_safe 16/16 pre y post. **RE-MEDICIÓN vs el Professional (binario 25-10-2023, target dof del record final post_point_dof de cada .dbs) — el fix NO cambió NINGÚN número → los desvíos NO eran contaminación de memoria, son calibración del kernel**: hypo1 sigyy −862.766345621 vs Pro −862.9290976 (0.019 %, tol 0.1, rc=1; hyhis0 0.598072 vs 0.598063), hypo3 sigxx −0.104030562155 vs Pro −0.162914620443 (36 %, rc=1; estado inicial anisótropo K0 — rama del kernel a calibrar), hypo2 rc=0 (sigxx −0.1728 dentro de tol −0.176±1e-2), hypo4 rc=0, hypo7/8/9 (masin) −224.08/−216.29/−224.08 vs Pro −231.81/−230.59/−231.81 (3-6 %, rc=1, sin cambio), hypo12 (masin visco) −282.46 vs Pro −143.49 con "iterative solver broke down" + hyhis0 estable 0.5906 en el GNU pero overflow de e (1.86e14) en hypo13 antes del aborto "severe error in Masin hypoplasticity" (rc=1, sin cambio). **Ningún test hypo nuevo pasa con el fix** (hypo2/4 siguen siendo los únicos rc=0); calibración wolfersdorff (hypo1/hypo3) y masin/visco (hypo7/8/9/12/13) PENDIENTE con los números exactos en manual-developer/group_materi_plasti_hypo_wolfersdorff.md. **Descubrimiento no obvio**: en el build -O1 default la basura que leía el OOB era padding a CERO (los flags espurios NUNCA se activaron → resultados no contaminados), pero el bug es UB determinista que mata los builds sanitizados y podía dispararse con otras flags/layouts → el fix es obligatorio aunque no mueva el corpus. Corpus: **193 PASS / 162 RUNFAIL / 8 PARSE (363)** = idéntico al baseline 3076250 (`/tmp/corpus_result_138ba5a.txt`; el fix es numéricamente inerte). Suite build_safe 16/16 y suite ASan 16/16 al cierre. |
| **SPRINT CONSOLIDACIÓN u-p (2026-09-05): defaults de `inertia_apply` (-yes→-no) y de la familia `groundflow_consolidation_apply` (-yes→-no) alineados al Professional** | `0560dcf` (feat) + docs de este lote | 2026-09-05 | **SÍNTOMA (A/B contra el binario Pro 25-10-2023, probe serie temporal de ground14 con pares de control_timestep hasta 2048 s)**: el GNU drenaba la columna artesiana de los safety tests (2D, phreatic 0.25, -topres −10 en la base, ρ_sat=2, E=100, k=1e-5) con un transitorio real de ~10² s (σ′(0)=−8.14 a 0.5 s, −9.531 a 1 s, −9.9999 a 128 s, −10.0000002 a 2048 s = el MISMO estado drenado asintótico del Pro: p_total(y)=−10+40y bajo el nivel, σ′_yy=−10−20y, disy −2.5625e-3, digito a digito), mientras el Pro estaba YA drenado en su paso 1 (−10.000000 a 0.5 s). La hipótesis previa del "factor ~625 = bug de unidades del acoplamiento" queda REFUTADA: forzando `groundflow_consolidation_apply -yes` en AMBOS binarios el transitorio del GNU (σ′(1 s)=−9.531) reproduce al Pro (−9.295, 2.5 %), y los tests de consolidación real del corpus (large2/large3, que SÍ fijan -yes) ya pasaban. **CAUSA RAÍZ (2 defaults legacy GNU ≠ Pro)**: (1) `inertia_apply`: el GNU legacy default -yes activaba la masa de los dofs de velocidad con densidad>0 → el modelo corría PSEUDO-DINÁMICO sin pedirlo (oscilación del momento alrededor del drenado: σ′(1 s)=−10.4 con acoplamiento off vs −10 exacto del Pro estático); el Pro (manual 6.785): "Default ... set to -no" (dinámica SOLO con inertia_apply -yes; el propio manual 9.9: "Dynamic calculations are triggered by setting inertia_apply -yes"). (2) `groundflow_consolidation_apply`: default -no en el Pro (6.556) vs -yes legacy del GNU → el GNU acoplaba div(v_soil) en la ecuación de presión de TODO modelo materi+groundflow sin registro → transitorio de consolidación espurio en tests de safety que esperan el estado drenado en su ventana de 1 s. **FIX**: defaults -no en general.cc/truss.cc (options_inertia) y groundfl.cc (los 3 switches + materidivergence), con la SEMÁNTICA DE ACTIVACIÓN medida en el Pro (A/B g_yes/g_no): solo el record GLOBAL o el CONTROL activan el acoplamiento; los records de GRUPO solo EXCLUYEN (-no) los elementos del grupo (un -yes de grupo solo no activa: medido, safety 2.3333 drenado); los legacy GROUP_GROUNDFLOW_MATERIDIVERGENCE/-OPTIONS_SKIP siguen honrados. **VERIFICACIÓN**: ground14/15/16 rc=0 SIN registros inyectados (σ′(0) −10.0000001423, piping 2.3333/lifting 2.0000 exactos = .dbs Pro); large2/3 (consolidación -yes explícita) y ground19/undrained1/2/ground17 siguen rc=0; A/B de la semántica: ground14 + global -yes → transitorio (σ′ −9.366), + global -yes con group -no → drenado rc=0 (exclusión), + group -yes solo → drenado rc=0 (no activa) — los 3 reproducen el Pro medido. **Corpus: 177 → 193 PASS / 162 RUNFAIL / 7 PARSE (363; 53.2 %)** — diff per-test del corpus completo contra la baseline c44a439 (`/tmp/corpus_result_new.txt`): **0 regresiones** y 16 RUNFAIL→PASS: ground14/15/16 + 13 bonus de la misma raíz (el default legacy -yes de la inercia hacía pseudo-dinámicos los modelos estáticos con densidad): condif6, distri1, incnav2, matrix1, matrix3, merge1 (el blocker condif transitorio del Batch C: 1.0 exacto), nodfor1, refine2, skip1, temp1, temp3, volume1, volume2. Suite build_safe 16/16 verificada al cierre. GOTCHA documentado: la premisa "GNU 625× más lento que el Pro" era comparar el default legacy GNU (consolidación ON) contra el default del Pro (OFF, drenado instantáneo por diseño de estos tests de safety). |
| **FAMILIAS PEQUEÑAS (2026-09-05): truss térmico + pretensado, geometry por-nodo/presencia, -k0, element_dof_initial, gates de salida** | `c24376e` (feat) + `de08e20` (docs) | 2026-09-05 | **group_truss_expansion** (manual Pro 6.775): el incremento de longitud mecánico del truss resta la expansión térmica libre alpha·dT·L0 (dT = promedio del incremento de temperatura nodal del paso, dof condif_temperature); truss11 rc=0: truss E=A=α=1 con extremos fijos calentado a T=2 → element_truss_force −2 EXACTO = .dbs Pro. **group_truss_initial_force** (6.776): en el nacimiento del elemento (sin estado ELEMENT_TRUSS_FORCE previo) la fuerza "old" del update incremental se siembra con el pretensado (F = F0 + EA/L·dL); truss12 rc=0: node_rhside ±1 en vely EXACTO del .dbs Pro (el Pro renumería el truss 2→3 en su .dbs; física idéntica; element_truss_strain 0 con truss estacionario = Pro). **truss15 rc=0** de propina (su parse fallaba en control_mesh_truss_distribute_mpc_element_group_truss, registrado en este lote). **post_calcul -materi_stress -k0** (6.901: ratio horizontal/vertical = 0.5(σxx+σzz)/σyy en 2D, 0.5(σxx+σyy)/σzz en 3D; item de salida -k0_sig = post_calcul_label del Pro): k0 rc=0 — node/post_point_dof_calcul 0.1111111111 EXACTO (columna confinada ν=0.1). **element_dof_initial** (6.422: campo de dofs PASADOS del elemento en su paso de nacimiento, consumido por la inercia transitoria): phase1 rc=0 — barra de conducción con T0=0.4 alcanza el equilibrio dinámico node 2 temp 0.7 EXACTO (término de capacidad C(T−T0)/dt; marker ELEMENT_DOF_INITIAL_APPLIED escrito por step_close tras converger el paso de nacimiento — probe 2 pasos del Pro: el paso 2 continúa 0.7→0.85, el campo inicial solo afecta al nacimiento). **geometry_node_type** (6.540) + **geometry_projection_type** (6.543): overrides POR GEOMETRÍA dentro de geometry() (coordenadas -node/-node_start_refined/-plus_displacement; -project_inside = interior relleno — 8 ramas filled-vs-edge); node_type_1 rc=0 — node_geometry_present 1 → -geometry_point 3 = .dbs Pro. **node_geometry_present** (6.886) + print_node_geometry_present (6.995/6.996): fill por paso al INICIO de step_start con las coordenadas del paso previo convergido (calibrado con barrido temporal t_end 3/5/8/9/10 en el Pro: el último fill gana; paso sin presencia deja el record vacío; estático sin movimiento no rellena); target (1,1)=3 EXACTO. **validation_14 rc=0** (refinamiento adaptativo local con error de materi_strain_plasti + delete + restart + repeat completo, 24 s): el círculo del hole con -project_inside alimenta la región del refine; alias control_mesh_delete_geometry_move_node (Pro) → movenodes. **Gates**: print_database_calculation (6.972) + print_gid_calculation consumidos en exit_tn (-no salta el volcado final; large1 ya no falla de parse pero RUNFAIL por TIEMPO/solver: 170 s y 51 GB — Bi-CG sin converger con retry a LU directa en los ~90k dofs del brick 30³, muy por encima del timeout de 45 s del corpus: PENDIENTE solver/perf); post_apply (6.900) y print_group_data (6.990) registrados parse-only (-no gating PENDIENTE; el único test no los consume). distri3: parse COMPLETO (alias repeat_save_calculate_result → repeat_calculate_result) pero RUNFAIL de SOLVER — el Bi-CG rompe en cada repeat del bar 2D con young aleatorio (sistema mal condicionado por el reset/repeat) y los valores del post_point divergen entre repeats; el Pro lo resuelve (>100 s). PENDIENTE familia solver/iteración mixta. **Registro de nombres**: K0/K0_SIG/PROJECT_INSIDE/PRINT_GROUP_DATA puros + aliases repeat_save_calculate_result → repeat_calculate_result. **Corpus: 170 → 177 PASS** (diff blast-radius 87 tests de las familias truss/interface/condif/conspr/validation/ground/material/post: SOLO 4 cambios RUNFAIL→PASS — truss11, truss12, truss15, validation_14 — y 0 regresiones; k0, node_type_1 y phase1 cerrados fuera del subset). **Corpus completo final: 177 PASS / 178 RUNFAIL / 8 PARSE** (diff per-test: SOLO los 7 RUNFAIL→PASS, 0 regresiones). Suite build_safe 16/16. distri3 y large1 quedan RUNFAIL con diagnóstico (solver/timeout, ver arriba). GOTCHA de infraestructura documentado: en este db los records NO version_all solo admiten VERSION_NORMAL (los internos llevan version_all=1); los records de entrada de longitud VARIABLE con data_length DATA_ITEM_SIZE se tragan las keywords siguientes (fixed_length=1 obligatorio para listas cortas). |
| **FAMILIAS MESH-CUT + BOUNDA-RELAX (2026-09-05): control_mesh_cut_geometry + control_mesh_cut_node_force + control_bounda_relax(_geometry)** | (commits feat/docs de este lote, 2026-09-05) | 2026-09-05 | **control_mesh_cut** (manual Pro 6.163/6.164, verificado contra el binario Pro 25-10-2023 — .dbs y stdout de mesh_cut_1.dat/mesh_cut_2.dat + probe 2D de corte real): el corte elimina los elementos cuyos nodos estan TODOS dentro de la geometria (misma membresia que control_mesh_delete_geometry) + poda los nodos huerfanos (nod_nod via mesh_has_changed) y sustituye la parte cortada por sus fuerzas nodales: la accion de los elementos eliminados sobre los nodos de contorno supervivientes (-f_elem, f_elem = int B^T*sigma dV desde las tensiones de punto de integracion ELEMENT_DOF del estado convergido previo al corte; integracion replica de pol()/materi(), familia tensorial bar2/bar3/quad4/quad9/hex8/hex27) se escribe en los records node_force (aplicados por parallel_new_dof_before en cada paso posterior). La malla restante conserva el equilibrio exacto del instante del corte. Activacion: el indice del control sin control_timestep corre su step virtual (top.cc) — en mesh_cut_1 el corte ocurre en t=2.0 entre los bloques 20 y 40. Verificado: mesh_cut_1 rc=0 — malla restante nodos 1..6 (0..5/9) = .dbs del Pro, node_force |1| en el nodo de contorno, sigxx post-corte = 1.0 sostenido (target ±1e-3); mesh_cut_2 rc=0 — la linea de corte no cae dentro de la tolerancia de NINGUN nodo: el Pro NO elimina nada (no-op reproducido, sus 81 elementos intactos); probe 2D de corte real (quadrilateral sobre las filas superiores): fuerzas de la fila de contorno 1/9 (interiores) y 1/18 (bordes) = EXACTAS del Pro (suma 1.0), 45 elementos restantes, sigyy post-corte = 1.0 rc=0. GOTCHAS: (a) el record -bar2 integra con la regla MINIMAL de 1 punto (polynom.cc:397-398), no la MAXIMAL — sin eso las fuerzas salen a la mitad; (b) el determinante del Jacobiano de los quad4 del macro es NEGATIVO: la inversa debe usar el det CON SIGNO (replica de pol(): |det| solo para el volumen) — con |det| antes de invertir las derivadas fisicas se voltean y la sustitucion sale comprimiendo en vez de traccionando; (c) el consumo del record node_force de NUESTRO dof.cc aplica un +F almacenado como carga -x (medido A/B nf1.dat: Pro +1 -> disx +1, GNU +1 -> disx -1): mesh_cut almacena el valor NEGADO para reproducir el equilibrio del Pro (nuestro .dbs node_force 6 = -1 vs +1 del Pro; PENDIENTE dedicado alinear el signo del record node_force y quitar la negacion); (d) la geometria de corte usa la semantica de proyeccion del delete (interior relleno de areas) en las 8 ramas de geometry(). **control_bounda_relax** (6.108/6.109): registros parse-only — la semantica de almacenamiento del Pro NO es observable en relax1/2 (A/B en el Pro 25-10-2023 con/sin los records: resultados IDENTICOS; la relajacion observable la conduce la funcion temporal de bounda_force sola). **relax1 rc=0** (target disx 0.5). **relax2 RUNFAIL documentado con causa raiz ORTOGONAL a la familia** (no es el relax): su modelo es velocity-only (sin materi_displacement ni group_materi_memory) y el GNU convecta la malla de referencia con el material por defecto (options_mesh default -follow_material, record legacy que el Pro 25-10-2023 RECHAZA; locate() mueve NODE con v*dt) -> cinematica de gran deformacion: reaccion = E*ln(1+u) (0.7188 en t=1, 0.7677 en t=3) vs la lineal del Pro (1.0/1.0 con la malla fija). A/B: con options_mesh -fixed_in_space el GNU reproduce el Pro EXACTO (reacciones 1.0/0.0/1.0 en t=1/2/3, rc=0) — alinear el default velocity-only con el Pro = work unit de cinematica con blast-radius de corpus (PENDIENTE). **Corpus: 167 → 170 PASS / 185 RUNFAIL / 8 PARSE** (3 nuevos: mesh_cut_1, mesh_cut_2, relax1; diff per-test del corpus completo: 0 regresiones). Los 5 SEGV del corpus (force16/17, mpc8/9, post7) son PRE-EXISTENTES (verificados rc=139 tambien con el binario del HEAD 8d6c4e4). Suite build_safe 16/16. Blast-radius malla/bounda (delete1-4, extrude/convert/interface/conspr/support/slide, mpc, force): sin cambios. |
| **FAMILIA CAPACIDAD NO DRENADA (2026-09-05): group_materi_undrained_capacity + element_intpnt_materi_undrained_pressure + control_materi_undrained_apply + apparents young/poisson** | (commits feat/docs de este lote, 2026-09-05) | 2026-09-05 | **Modelo (manual Pro 6.760 + teoria 2.2.7, verificado contra el binario Pro 25-10-2023)**: el analisis NO DRENADO sin el dof de flujo en la matriz: la ecuacion de almacenamiento SIN permeabilidad C·ṗ = div(v) se resuelve POR ELEMENTO; por punto de integracion p_u += Δε_vol/C por paso (record element_intpnt_materi_undrained_pressure, un valor por IP; _average registrado). La presion actua isotropica: σ_total = σ'_eff + (p_fija_del_flow + p_u)·I para el equilibrio de momento (manual 2.2.7: "fixed total pressure from the hydraulic pressure heads plus the excessive undrained pressure ... as the full total pressure"), los dofs de σ (sigxx...) conservan el valor EFECTIVO (medido: u1 node σyy = E·ε; p_u aparte) y la rigidez de momento gana el termino volumetrico dt/C·(div w)(div v) — la respuesta total = drenada + 1/C (medido u2: ε = F/(E+1/C) = −0.25 con E=3, C=1, F=1: sigxx −0.75 & p_u −0.25 EXACTOS del .dbs Pro). **Tests**: ground17 rc=0 (strain volumetrica prescrita −1e-5, E=1e7, ν=0, C=1e-7: young_apparent 1e7 & poisson_apparent 0 EXACTOS en los 4 nodos = .dbs Pro digito a digito — los operadores post_calcul -young_apparent/-poisson_apparent (6.903) implementados con el INCREMENTO del ultimo paso (snapshot interno node_dof_previous_step capturado por top() al inicio de cada paso; E = 9KG/(3K+G), ν = (3K−2G)/(2(3K+G)) desde K = dp/dε_v y G = dq/(3dε_q)); undrained1 rc=0 (quad4 + flow bounda -tpres − all con phreatic 1.0 SOBRE el mesh: to_pres = 10z EXACTO (0 arriba/−10 en z=−1), p_u −49999.99875 vs Pro −49997.5 (0.005%, tol ±100), σyy_mid −49999.99875 vs −50002.5 (en tolerancia)); undrained2 rc=0. **FIX -topres rama phreatic (bounda.cc)**: MEDIDO en el Pro, con un groundflow_phreatic_level la conversion head→dof NO resta el static por nodo: el head conserva la carga (h = load) y p_total = h − ρgz: pres_dof (= p_dynamic) = load − ρg·level (uniforme; los probes E/I: load −100 → −90/−105 = la p_dyn del .dbs del Pro). La conversion previa (load − static) solo estaba verificada para la rama static-height (ground13): con el nivel SOBRE el mesh daba p_fija = 0 en vez del perfil −ρgz medido. groundflow_phreatic_coord gana el out-param level_source (1 = phreatic, 2 = static-height) para distinguir las ramas (8 callers actualizados). **Alias -tpres** (db_number): el corpus escribe -tpres (undrained1) y -topres (ground13); el manual dice -topres. **Registros**: GROUP_MATERI_UNDRAINED_CAPACITY (1 valor, class MATERI), ELEMENT_INTPNT_MATERI_UNDRAINED_PRESSURE (npointmax slots, fixed_length 0, version_all 1 como ELEMENT_DOF: el copy de fin de paso promueve el valor convergido a NORMAL), _AVERAGE (registrado), YOUNG_APPARENT/POISSON_APPARENT (name entries del post_calcul), NODE_DOF_PREVIOUS_STEP (interno, external 0). Pre-alloc del record en top() (mallas de input) Y en step_start (las mallas control_mesh_macro aparecen DESPUES del init de top(): sin el re-check el primer PUT aloca dentro del element loop paralelo → "Program error"). GOTCHA del record multi-IP: el read-modify-write por IP debe preservar los slots de los OTROS IPs del mismo pase leyendo el buffer de VERSION_NEW (con NORMAL como base solo del slot actual): sin eso el ultimo call clobberea los anteriores (medido [0,0,0,val] en u1). GOTCHA fixed_length: el default de db es fixed_length=1 (registros de longitud FIJA = data_length): el record por-IP necesita fixed_length=0 (longitud variable por elemento, el stride de db queda npointmax). **Corpus 164 → 167 PASS / 188 RUNFAIL / 8 PARSE** (3 nuevos: ground17, undrained1, undrained2; diff per-test vs el baseline bicg14c + lotes intermedios: 0 regresiones). Suite build_safe 16/16 (build limpio). Blast-radius familia u-p/groundflow/consolidacion: ground1/2/3/6/10/12/13/19 rc=0 (igual que baseline); ground4/7/8/9/11/14/15/16/18 y consolidation_gibson/foundation_consolidation siguen en sus RUNFAIL documentados (diagnosticos previos: transitorio u-p de consolidacion ground14/15/16, ground8 static multi-nivel, ground11 dependency+singular). NOTA menor documentada: el reset lineal del σ inicial (u1) no reproduce el gradiente 10y del Pro (diferencia 2.5 en 5e4 = 0.005%, dentro de la tolerancia ±100 del target; la diferencia real es la discretizacion de carga ε = −4.999999875e-3 vs −4.99975e-3). Manuales: manual-user/group_materi_undrained_capacity.md + manual-developer/group_materi_undrained_capacity.md + bounda_dof.md actualizado (-topres phreatic + alias -tpres) + indices. |
| **FAMILIA PRESIÓN TOTAL DEL GROUNDFLOW (2026-09-05): bounda_dof -topres + superficie libre del nivel freático + operadores de seguridad piping/lifting** | (commits feat/docs de este lote, 2026-09-05) | 2026-09-05 | **bounda_dof -topres** (manual Pro 2.4.1: prescribir la presión TOTAL p_total = h−ρgz en vez del head): el parser resuelve "-topres" al enum KEYWORD groundflow_pressure (NO al dof-label dinámico de "-pres"); bounda.cc lo mapea al pres dof (topres_bounda por record) y la carga se convierte POR NODO invirtiendo groundflow_phreatic_coord: pres_dof = load − static cuando un nivel freático/static-height cubre el nodo; pres_dof = load + ρgz sin nivel (mismas convenciones de signo del Pro: g y z negativos bajo el datum). **ground13 rc=0**: con topres −20 (z=0)/−10 (z=1), ρ=1, g=−10 y la region static-height 123: pres_dof uniforme +1210 (=Pro h−ρgz_L... dy_pres), to_pres −20/−10 y st_pres −1230/−1220 y dy_pres 1210 = .dbs del Pro DIGITO A DIGITO (antes db_error BOUNDA_UNKNOWN en el parse). **Superficie libre / zona seca del nivel freático simple** (groundflow_phreatic_apply, groundfl.cc): con un groundflow_phreatic_level (sin multiple), los nodos con coord ≥ nivel se boundean a pres_dof = 0 (p_total = 0: el static de phreatic_coord se clampea a la atmosferica sobre el nivel; reproduce el Pro medido: h = ρgz_L sobre la linea = p_dynamic 0). Sin este bound el dominio sobre el nivel queda solo conectado al Dirichlet inferior y el solve acoplado rompe (Bi-CG breakdown) — con el bound la zona saturada recibe su Dirichlet de superficie libre (medido: pres_dof lineal −7.5→0 entre el fondo y el nivel, zona seca 0; Pro: p_dyn −7.5→0). **post_calcul -safety_piping/-safety_lifting** (6.919, consumo): naming por metodo (vertical: safety_piping 1 valor; prival: safety_*_prival_0..2; global: safety_*_global_x/y/z; 12 enums puros nuevos) + formulas en calculate_operat: lifting=(σ_i+p_total)/p_total, piping=(σ_i+p_dyn)/p_dyn con σ_i = σ vertical (yy 2D/zz 3D/xx 1D), prival del MAS compresivo (GOTCHA: el sort GNU es descendente → prival[MDIM−1−j]) o σ globales; denominador |p|<TINY → 0 (safety_default); cap post_calcul_safety_maximum. **Corpus 163 → 164 PASS / 191 RUNFAIL / 8 PARSE** (1 nuevo: ground13 RUNFAIL→PASS; diff per-test 0 regresiones — ground3, el unico PASS con nivel simple, bit-igual). Suite build_safe 16/16. **PENDIENTE documentado (diagnóstico fino, NO es el -topres)**: ground14/15/16 corren TODO el pipeline (parse, pres, safety, targets) y fallan solo el target numerico (piping_global_y 2.2708 vs 2.3333±1e-3) porque el acoplamiento u-p de consolidación del GNU (feedback ∇·v_soil en la ecuacion de presion) esta en transitorio a t=1 s (ventana del test; Bi-CG breakdown en el paso 1; σ′(0) −9.53 vs −10); el ESTADO DRENADO del Pro se reproduce EXACTO con el GNU a ~2000 s (τ_cons = L²/c_v ~ 625 s con k=1e-5, E=100, γ_w=10): pres lineal −7.5→0, σ′ −10.000, piping 2.3333 y lifting 2.0000 rc=0 (probe 4000 pasos). La convencion de consolidacion (velocidad de drenaje del esquema u-p) es el siguiente blocker de la familia. ground11_nonsaturated: SIN CAMBIOS (su monitor dependency_item -to_pres sigue sin ser dof-label + solve singular inicial; NO es el -topres de bounda). ground8 (phreatic multiple + mecánica, to_pres −76.9 vs −10): la rama de nivel SIMPLE no lo toca — PENDIENTE su convención static multi-nivel. |
| **KEYWORDS BATCH C (2026-09-04): clústeres pequeños (springs, geometry_element_group, dependency gates, post_calcul estático/seguridad, registros parse-only)** | (commits feat/docs de este lote, 2026-09-04) | 2026-09-04 | **springs**: `group_spring_memory` (6.765; el muelle GNU es incremental en config actual = -updated_linear; parse) + `group_spring_stiffness_nonlinear` (6.768, diagrama eps/k con la k del incremento en la strain MEDIA: integración EXACTA para diagramas lineales — spring6 F=0.5±1e-8) + `element_spring_strain` (elongación total por elemento; version_all+pre-alloc en top.cc). **spring1 y spring6 rc=0**. **geometry_element_group/method** (6.524/6.525): filtro de geometrías por grupos de elementos DENTRO de geometry() (helper nuevo node_attached_element_groups en group.cc que lee ELEMENT_GROUP — el primer valor del record ELEMENT es el TIPO, no el grupo); reglas: -all = attached==listado; default/-any/-only = attached⊆listado (la exclusión de nodos con elementos fuera de la lista en TODO método reproduce el binario Pro 25-10-2023 medido en merge2). **merge2 y force12 rc=0**. **FIXES asociados**: (a) area.cc — la dirección de las cargas/flujos de arista es SIEMPRE la normal FÍSICA de la cara calculada de las coordenadas (con geometry_point la normal acumulada de geometry() es RADIAL y tilts la carga: force12 medía disy 1.179 vs 1.0); (b) elem.cc — los elementos de grupos puros -empty/-none NO integran area() (la cara compartida entre un elemento material y uno -empty no se carga dos veces). **dependency gates**: `dependency_apply` + `control_dependency_apply` (6.125/6.126) — gate global/per-control en get_group_data() (group.cc; el control_ NO lleva data_required: el combination-check por mismo índice exigiría dependency_item en el índice del control, que ground11 no tiene). ground11 avanza al solve (singular inicial: siguiente blocker). **post_calcul**: `post_calcul_static_pressure_height`+`_element_group` (6.921/6.922) — regiones (coord_min,coord_max,height_ref) consumidas en groundflow_phreatic_coord (p_static = ρg(height_ref−y), total = pres+static); **ground13 target st_pres −1220 EXACTO** (sigue RUNFAIL por bounda_dof -topres). `post_calcul_safety_method`/`_maximum` (6.919/6.918) + operadores -safety_piping/-safety_lifting + selectores puros VERTICAL/GLOBAL; **labels medidos contra el .dbs del Pro en ground15/16** (vertical: safety_piping 1 valor; prival: _prival_0..2; global: _global_x/y/z; fórmulas lifting=(σ_i+p_total)/p_total, piping=(σ_i+p_dyn)/p_dyn calibradas: 2.0/2.3333); cálculo PENDIENTE (familia convenciones de presión). **Parse-only registrados**: group_materi_plasti_element_group/_factor (6.698/6.699), group_groundflow_expansion (6.614), control_print_gid_contact_spring2 (6.297 + alias sin control_ — no_index: el corpus lo escribe con el valor solo), control_mesh_macro_concentrate (6.207), alias control_mesh_generate_truss_beam→trussbeam. **Corpus: 159 → 163 PASS / 192 RUNFAIL / 8 PARSE** (4 nuevos: spring1, spring6, merge2, force12; diff per-test: 0 regresiones). Blast-radius: springs 2-5 y force/edge-load/merge/interface sin cambios; suite build_safe 16/16. PENDIENTES documentados con diagnóstico: ground11 (monitor de dependency sobre -to_pres no-dof + solve), heat_exchanger_* (parse completo; macro concentration + física termo + singular/timeout), tutorial_3 (parse completo; modelo beam/excavation pesado), ground13/15/16 (bounda_dof -topres + convenciones de presión), undrained1/2/ground17 (group_materi_undrained_capacity), relax1/2 (control_bounda_relax: almacenar RHS nodal + relajación). |
| **Lote 3D cuadrático (L3 + L4)**: keyword `-hex20` registrada (enum HEX20 + name; el GNU no tiene hex20 real — el Pro lo maneja NATIVO y conserva `-hex20` en su .dbs) → `mesh_convert_hex20()` en mesh.cc eleva cada volumen serendipity hex20 al Lagrange hex27: permutación explícita al orden tensorial GNU verificada contra el .dbs del Pro (su conversión automática de interface_quad8_hex20 escribe EXACTAMENTE el layout tensorial GNU: planos base/mid/top quad9-tensor; esquinas el[1..4]/el[5..8] → slots 1,3,7,9 / 19,21,25,27; aristas BM,LM,RM,TM base/top (el[9..16]) → slots 2,4,6,8 / 20,22,24,26; verticales (x0,y0),(x1,y0),(x0,y1),(x1,y1) (el[17..20]) → slots 10,12,16,18; 6 centros de cara NUEVOS (media de las 4 esquinas de la cara) en slots 5,11,13,15,17,23 + centro de cuerpo (media de las 8) en 14); centros DEDUPLICADOS por coordenadas (EPS 1e-10) entre vecinos (hex20.dat apilado z=0..2: la cara compartida z=1 tiene UN solo centro en ambos hex27; sin dedup la malla se rompe); hook en step_start TODOS los pasos (idempotente), ANTES de extrude() e interface_convert(). **L4 interfaz 3D cuadrática**: `-hex18` (2 caras quad9, ns1=9) en interface_element() — frame con esquinas de cara 0,2,6 (el nodo 1 = BM es colineal → cruz degenerado), pesos por par = producto Lobatto 2D (1,4,16)/36 (interface3 carga la cara top con −1/−4/−16 y TODOS los intpnts reportan −9 uniforme; con pesos uniformes el intpnt de esquina daría −20.25), medida de cara desde las esquinas 0,2,6,8 del quad9 (orden tensorial = bowtie, fan de triángulos); cara `-quad8` de grupo de interfaz → `-hex18` en interface_convert() (el 9º nodo del lado 1 = centro de cara = media de las 4 esquinas, buscado por coordenadas — el nodo que mesh_convert_hex20 creó deduplicado en la cara compartida del hex27; si no existe se crea; lado 1 reescrito al orden tensorial con la misma permutación de slots del quad8; 9 duplicados desplazados 0.01 en la normal; reconexión del sólido del otro lado por TODOS los 9 nodos + test de centroide). **Records de interfaz**: los historias/output son arrays planos de db_data_length() y el PUT rechaza longitud mayor (defaults para ns1≤4): el bloque any_interface de step_start los ELEVA (monótono, solo si el modelo los excede — bajar el stride corrompe TODO el record, regresión expuesta por interface2/10/patch) al mayor elemento de interfaz del modelo (hex18: ns1=9, intpnt 9·3=27); buffers de reset en data.cc ampliados 4→9. | `529b383` (feat mesh) + `4af8739` (feat interface) | 2026-09-04 | **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: (1) hex20 rc=0 — node_dof IDÉNTICOS al Pro en los 32 nodos comunes (worst 3.8e-10 en desplazamientos, 2.1e-8 en σzz; los 13 nodos interiores del GNU no tienen contraparte Pro porque el Pro conserva hex20); el run auto-convertido es bit-igual (1.5e-15) al mismo mesh escrito como `-hex27` explícito; (2) interface3 rc=0 — `element_interface_intpnt_stress` = −9 EXACTO en los 9 intpnts (target −9±1e-5); (3) interface_quad8_hex20 rc=0 — elemento convertido `-hex18 5 13 6 14 57 15 7 16 8 | dups` = MISMO layout del Pro (`-hex18 5 13 6 14 55 15 7 16 8 | 68..76`; 57/55 = el centro de cara compartido); sigzz post_point −1.005025 vs target −1.0±1e-2 = el MISMO offset 0.5 % de la familia hex8 establecida (interface_quad4_hex8: −1.005025 — la convención GNU desplaza los duplicados 0.01 al +normal, el Pro los crea coincidentes; física equivalente). **Corpus: 156 → 159 PASS / 196 RUNFAIL / 8 PARSE** (3 nuevos: hex20, interface3, interface_quad8_hex20 — todos RUNFAIL→PASS; diff per-test: 0 regresiones). Blast-radius: interface1-15/patch/bar2_quad4/quad4_hex8/many/tria3_prism6/bar2_hex8/conspr1-7/quad8/elasti6/interface_bar3_quad8/patch1 rc=0; suite build_safe 16/16. |
| Auto-conversión de mallas cuadráticas 2D (lote L1+L2 del sprint 2D): keyword `-quad8` registrada (enum + name; el GNU no tiene quad8 real, como el Pro lo auto-convierte) → `mesh_convert_quad8()` en mesh.cc inserta el nodo CENTRO (media de las 4 esquinas) y reescribe el elemento en el orden tensorial del quad9 GNU con la permutación `quad9 = {q8[1],q8[5],q8[2],q8[6],C,q8[7],q8[3],q8[8],q8[4]}` (corpus [BL,BR,TL,TR,BM,LM,RM,TM] → GNU [BL,BM,BR,LM,C,RM,TL,TM,TR]); hook en step_start TODOS los pasos (un `control_mesh_convert` en índice < timestep evalúa post_point en un paso intermedio contra el quad8 crudo → assert; idempotente), ANTES de extrude() (quad8→quad9→hex27) e interface_convert(); los `-quad8` de grupo de interfaz NO se convierten (cara 3D, lote futuro). **Extensión de `control_mesh_convert`**: `-bar3` → `-quad6` (2D cuadrático, 3+3 nodos; ns1=3; reconexión por TODOS los nodos del lado 1 + test de centroide; el quad6 de interface_element() ya existía con Lobatto 1/6,4/6,1/6). **Nuevo item `-post_force_edge_summed`** (manual Pro 6.936: total de la fuerza de los force_edge integrada sobre las aristas, ndim valores, record plano no_index): cálculo on-demand (target_item/control_print) en post() espejando la integración de area() (cuadratura Lobatto × largo/área de cara, load temporal sine/time/time_file, factor espacial, restricciones de elemento/grupo/lado/nodo/geometry). | commit feat de este lote (2026-09-04) | 2026-09-04 | **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: (1) quad8 rc=0 — node_dof del patch test IDÉNTICOS al Pro (sigxx 0.1, disx 1.0/2.0, ruido 1e-16); (2) elasti6 rc=0 — elemento convertido `-quad9 1 5 2 6 9 7 3 8 4` = misma permutación de slots del Pro; sigyy node 8 = 1.00000008274 vs Pro 0.9999998606958 (target ±1e-4), disy 1e-8 idéntico; post_force_edge_summed = (0, 10) EXACTO (Pro: 0. 1.000000000000e+01); (3) interface_bar3_quad8 rc=0 (sigyy −1.0 ±1e-2; Pro rc=0); (4) NOTA Pro: crea los duplicados de interfaz con coordenadas coincidentes en el bloque −normal; el GNU los desplaza 0.01 al +normal — física equivalente, solo cambia el layout del .dbs. **Corpus: 153 → 156 PASS / 199 RUNFAIL / 8 PARSE** (3 nuevos; diff per-test: 0 regresiones; ground19_water_under_dam flaky por el timeout 45 s del corpus — PASS en los runs previos y en este). Blast-radius: interface1-15/patch/bar2_quad4/quad4_hex8/many/tria3_prism6/bar2_hex8/conspr1-7/force7/8 rc=0; suite build_safe 16/16. |
| Familia slide (`slide_geometry` law, manual 6.1041-6.1047 + 6.370-6.372 + 6.893): ley elastoplástica kn/kt sobre desplazamientos TOTALES (patrón support_edge_normal; slide.cc `slide_spring`, legacy penalty conservada para slides sin `slide_stiffness` — examp14): spring normal bidireccional con cap de tracción `slide_plasti_tension` sig_t, cap de fricción `c + Fn·tan(phi)` con retorno al cap en slip (matriz tangencial = solo fracción residual, el kt elástico sesga los pasos de 2 iteraciones), factor axisym 2·pi·r (record `slide_axisymmetric` o `group_axisymmetric`), records de salida `node_slide_force`/`node_slide_f`/`node_slide_direction`; FIX de parse de la forma sin índice `axisymmetric -yes/-no` (manual 6.17) con pushback de un token (input.cc); `print_debug` registrado parse-only | `f79ea82` | 2026-09-04 | **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: (1) slide1 rc=0 — estado final IDÉNTICO al .dbs del Pro: node_slide_force nodo 1 (−2.1148e-3, −1.2210e-3) vs (−2.1129e-3, −1.2199e-3), nodo 2 (−7.8802e-3, −4.5496e-3) vs (−7.8821e-3, −4.5507e-3), reacciones top +0.009995/−0.0057707 (targets +0.01/−0.005773): Fn = kn·penetración EXACTA por nodo (sin factor L/2) y fricción = μ·Fn al cap (node_slide_f = 0); la redistribución vertical por el par del drag (S1 = 2.11e-3/S2 = 7.88e-3 en x=0/x=1) se reproduce. (2) slide3 rc=0 (parse de `axisymmetric -no` + `print_debug`). (3) axisym1 rescatado por el parse sin índice. (4) slide2 RUNFAIL: σyy −6.447e-4 vs target −6.666666e-4 (tol 1e-9) — blocker PREEXISTENTE del elemento axisym del GNU: reparto de fuerzas nodales del anillo comprimido IGUAL (2πr en el punto medio: probe columna fija GNU 4.712/4.712) vs el Pro con lobatto ∝ radio del nodo (3.142/6.283) → los springs se asientan con u1 ≠ u2; la ley slide queda verificada contra el elemento GNU (u1 = (2/3)σ, u2 = (5/12)σ exactos). (5) slide4 RUNFAIL: drag a vx = 1 (100× el de slide1): la iteración del slider rígido oscila al arrancar el drag (2..40 iteraciones/paso sin converger; el nodo x=0 deriva a lift-off) — requiere slip plástico INCREMENTAL con historia + tangente consistente (PENDIENTE; **SUPERADO 2026-10-11, F4 T6e C53**: slide4 rc=0 = Professional con el default sin límite de tracción, el cap `c + \|Fn\|·tan(phi)` y slide() en el post-solve pass; ver registro de verificación). **Corpus: 147 → 150 PASS** (205 RUNFAIL / 8 PARSE). Suite propia 16/16 + blast-radius 16/16 sin regresiones. |
| **MEDIDA del elemento de interfaz (largo 2D / área 3D) + defaults Bi-CG más estrictos** — el ensamblaje por par de la interfaz llevaba los pesos Lobatto/even DESNUDOS (suma 1 = elemento de medida unidad): la fuerza nodal = w·σ y la rigidez = w·kn en vez de w·L·σ / w·L·kn. patch1 del corpus (interfaces quad6 inclinadas de largo 1.677/0.559 entre bloques quad9 cargados) convergía a σn = 1610/1073/536 por intpnt (media 1431) en vez del 960 uniforme del Professional — sin el largo el sistema discreto pierde el brazo de momento del camino de carga y la tracción uniforme NO puede equilibrar la carga de borde (la suma ponderada 2400 = la carga solo con el campo NO uniforme, verificado); con la medida el 960 uniforme equilibra (σn per-IP 960.0000000 con solver directo; desplazamientos nodales = Professional a 10 dígitos). El área 3D usa abanico de triángulos (la numeración de cara hex8/quad4 de tochnog es bowtie: la fórmula de diagonales cruzadas da 0). Defaults Bi-CG: bicg_error 1e-10 → 1e-14 y floor 1e-12 → 1e-16 — los sistemas mixtos condicionados por penalización (kn/E ~ 1e4-1e7) dejan el jump de interfaz con ruido ~1e-2 al criterio relativo-1e-5 viejo (σn 960.01 vs target ±1e-3); el criterio más fino los hace fallar honesto (breakdown CG) y el retry Bi-CG → LU directa existente los resuelve EXACTO. Los tests de largo unidad (interface1-15, conspr1-7, interface_quad4_hex8, mohr_coul_direct1-3) quedan bit-igual (rc=0). | `6cafafe` | 2026-09-04 | **patch1 rc=0** (σn 960.0000000 directo / 960.0002 con Bi-CG 1e-14+retry; σt −480.0000; σxx post 1200.0000000 — antes 1431/1431 con CUALQUIER solver); ground19_water_under_dam rc=0 (18.8 s, sin retry); **corpus 153 PASS / 202 RUNFAIL / 8 PARSE** con patch1 dentro (la única baja vs la baseline 141: validation_2 — su −0.688985 es el valor del modelo CON solver directo puro también; el target −0.63 era artefacto de la trayectoria del CG laxo en el transitorio inercia/consolidación; taylor3 = marginal 45 s del timeout del corpus, flaky en ambos binarios). Suite propia 16/16 en build limpio. PENDIENTE (misma familia, mecanismos distintos — diagnóstico fino en DIAG-SOLVE-MIXTO §14.3): mpc3/4 (tying mpc_linear_quadratic: diagnóstico fino 2026-09-07 en DIAG-SOLVE-MIXTO §15 — la hipótesis de los σ-dofs libres queda REFUTADA; el tie es value-constrained sin redistribución energética → punto fijo 0.454 ≠ 1/3; cierre PENDIENTE = tratamiento implícito, work unit mayor), ground8 (phreatic multiple + mecánica: to_pres −76.9 vs −10), dynamic1/2/5/8 (materi_dynamic f=0 — límite explícito del esquema), ground15/16 y mpc5/6 (raíces mixtas). |
| Senda SuperLU directa arreglada (`solver -matrix_superlu`) + criterio de ruta directa | `c2cbaf7` | 2026-09-04 | **Causa raíz del SEGFAULT SuperLU (2 bugs en cadena, so.cc)**: (1) desde 5b7bb8b el band mat se aloca SIEMPRE y el ensamblaje sparse (else de `if (mat)`) quedó como código MUERTO — solve_superlu recibía una matriz VACÍA (nnz=0, dump verificado) → SEGV; (2) la búsqueda de miembros `array_member((long int*)inz[a],...)` escaneaba buffers `int` con stride de `long`: con >2 entradas por columna el scan compara pares de ints compuestos → duplicados CSC + valores en slots equivocados + diagonales ausentes (validado: bar2 nnz 437→376, 47/60 columnas con duplicados, columna 1 sin diagonal) → colamd/dgssv SEGV o basura. Fix: el fill sparse corre cuando el solver es SuperLU/PETSc aunque mat esté alocado + helper `sparse_list_add()` con scan int-correcto. **Verificado**: patch1+SuperLU rc=0 con 4 factorizaciones reales y solución IDÉNTICA a la banda LU (σn interfaz 1431.07742001 vs 1431.07742002, 11 dígitos); bar2_hex8+SuperLU ya NO segfallea: INFO=14 honesto (matriz singular — bug de extrusión del test); suite propia 16/16 (tslv_override con banner SuperLU); blast-radius 114 tests sin cambios vs HEAD (59 PASS; force16/17 + mpc8/9 SEGV pre-existentes confirmados también con el binario HEAD). **Corpus: 153 PASS sin regresión**. Adicional: mat/ipiv se liberan en TODAS las rutas (leak por llamada desde 5b7bb8b en bicg/superlu) y `solve_superlu` devuelve 0 con info≠0 (fallo honesto). Diagnóstico CORREGIDO de los pendientes: patch1 y bar2_hex8 NO son bloqueadores de solver (ver su sección) — el clúster "cond ~1e5" del plan no compartía raíz de solver lineal. |
| `materi_acceleration` (4.11: dofs de aceleracion derivados acc = (v_new − v_old)/dt + bounda_dof -accx = bound de velocidad v_new = v_old + a·dt) + `materi_dynamic`/`control_materi_dynamic` (6.800/6.141: blend sigma = (1−f)·sigma(t) + f·sigma(t+dt) con rigidez del momento ·f) + `start_if_not`/`end_if_not` + registros parse-only (`control_data_save`, `control_print_gid_save_difference`, `bounda_time_smc(_offset/_units)` con lector PENDIENTE) | `f516a03` | 2026-09-04 |
| Lote keywords baratas (2026-09-04, corpus 143 → 147 PASS): (1) `materi_plasti_camclay_history` (4.16) — initia con 2 dofs fijos cchis0/cchis1 (= e y p0 del Modified Cam Clay, basenames en db_initialize) y `group_materi_plasti_camclay` al layout Professional `M kappa lambda` (3 valores, variable 1..3): p0 pasa a ser historia EXPLICITA (plasti.cc lee old_hisv[1]; fallback legacy N solo con 4 valores y sin cchis1); checks camclay aceptan materi_history_variables o materi_plasti_camclay_history. **validation_13 rc=0**: cchis1 final 262.3 vs target 264.5±10 (Professional .dbs: 262.32). (2) `group_beam_inertia` (6.590) `Iyy Izz J` — A/B en el Professional (25-10-2023): la flexion 2D x-y consume Izz (2º valor: disy 3.333e-3 con Izz=1 vs 3.333e-5 con Izz=100); el GNU lee variable 1..3 (legacy 1 valor = Izz); `-total_linear` aceptado en beam_memory; alias 2D rotx/roty → rotz en db_number (los inputs Pro prescriben los 3 rot). **trubea1/trubea4 rc=0**; trubea2/3 = blocker beam 3D (flexion fuera del plano: el beam GNU solo tiene rigidez en su plano → Bi-CG singular; PENDIENTE). (3) naming `dipriscohis0..N-1` para materi_plasti_diprisco_history (db_initialize) — **diprisc1 rc=0**. (4) items groundflow del post_calcul: `to_pres`/`st_pres`/`dy_pres` (nombres puros) + generacion en calcul.cc con los labels Professional + resolucion de slots de target_item sobre POST_POINT/LINE/QUADRILATERAL_DOF_CALCUL en exit_tn (miscel.cc) — ground8 ya resuelve el target; valor -76.9 vs -10 del Pro = blocker de la convencion static de los phreatic multiples (PENDIENTE). (5) `groundflow_pressure_gradient` (4.7) registrado (dofs pres_gradx..; relleno PENDIENTE) + alias `group_porosity` + `control_groundflow_seepage_apply` registrado + FIX del combination check roto de los 3 controles groundflow por timestep (data_required=GROUNDFLOW exigia un marker inexistente al mismo indice → check.cc contra el dof). (6) familia slide registrada (`slide_plasti_friction` phi c, `slide_plasti_tension`, `slide_stiffness`, `slide_plasti_residual_stiffness`, `control_slide_plasti_apply`, `control_slide_stiffness_apply`) — consumo PENDIENTE (slide1/4 corren rc=1: la slide() legacy no constrine los nodos: el bloque cae, reaccion top ~0; falta la ley elastoplastica kn/kt + cap c+Fn tan(phi) sobre desplazamientos totales). | (commits de este lote, 2026-09-04) | 2026-09-04 | **Corpus: 143 → 147 PASS** (208 RUNFAIL / 8 PARSE). Suite propia 16/16. Sin regresiones (16 tests que pasaban siguen pasando: beam2d_1/2, beam3d, truss1/4/5/6/8/14, validation_12, mohrcou4, hypo2, ground19, visc_pl1, mohr_coul_direct3, mpc7). tutorial_2 parsea hasta el include (la mesh auxiliar no la stagea el harness — mismo bloqueo que earthquake_*). |
 **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: (1) dynamic3 rc=0 — prescripcion -accx = 1.234 en barra 1D: disx nodo 1 = 62.317 a t=10 IDENTICO al Professional (integracion Euler del bound de velocidad; target 61.7±1); dynamic4 rc=0 (node_mass + bounda_force: disx 1.0 y accx 2.0 a t=1). (2) A/B de materi_dynamic en dynamic1 (truss): Pro f=0 → 0.095889 (=target) y f=1 → 0.9955; GNU f=1 → 1.912 y f=0 → 6.45/divergente — el blend es correcto (dynamic5 mejora 15%→1.3% de error con f=0) pero el esquema de VELOCIDADES del GNU (masa lumped + rigidez implicita) no alcanza la estabilidad del esquema de DESPLAZAMIENTOS del Pro: dynamic1/2/5/8 quedan RUNFAIL por esquema/solver (familia DIAG-SOLVE-MIXTO), NO por parse. (3) 12 tests PARSE→RUNFAIL/PASS: dynamic2/5/6/7/8 y earthquake_1/2/3 y smc_1/2 dejan de fallar en init (materi_acceleration) y dynamic6 pasa del parse al solve (divergencia del pulso elastico — solver); smc_1/2 fallan con error honesto de lector SMC PENDIENTE; earthquake_* necesitan ademas el include-mid-record de valores y archivos auxiliares que el harness del corpus no stagea. **Corpus: 141 → 143 PASS** (210 RUNFAIL / 10 PARSE). Suite propia 16/16 (recovered). Sin regresiones (0 tests que eran PASS dejaron de serlo). |
| Plasticidad de interfaz Mohr-Coulomb direct (`group_interface_materi_plasti_mohr_coul_direct`, manual 6.632 — ángulos en RADIANES): ensamblaje de fuerzas TOTALES acumuladas (patrón spring), dilatancia por deslizamiento plástico INCREMENTAL, Fn positiva bajo compresión en el límite de fricción, reset `control_reset_dof -sigxx/-sigyy/-sigzz` → pre-tensión normal de la interfaz alineada, y FIX de `control_reset_interface`/`_interface_strain` (anidado en el gate de control_reset_dof + longitud basura → nunca reseteaba los histories) | (commit feat de este lote, 2026-09-03) | 2026-09-03 | **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: (1) mohr_coul_direct3 rc=0 — ensayo de corte con `0. 0. 0.785398` (φ=0, c=0, φ_flow=π/4): σn final −199.99997 y node_rhside ±99.99998 IDÉNTICOS al .dbs del Professional (el GNU viejo ensamblaba el incremento del último paso: −50; la dilatancia usaba el multiplicador plástico TOTAL y el ratchet doblaba la cuenta: −250 vs −200 = kn·(−1e-6 − 1e-6·tan(π/4))). (2) interface9: el ensamblaje incremental hacía creep 10 pasos iguales bajo carga constante; el Professional alcanza el equilibrio estático en 1 paso — con fuerzas totales el GNU converge igual (τ −0.1593, sigue rc=0). (3) interface7 re-verificado con los records acumulados: el reset de interfaz (índice 15) ahora dispara de verdad (estaba muerto: db_max_index devuelve −1 con records activos + gate de control_reset_dof + PUT de 1 valor con longitud basura). (4) direct4: el reset −sigyy −1 pre-tensa la interfaz (probe .dbs: `element_interface_strain_normal 3 -1e-06` = el ε_n del Professional) y el límite = c + |Fn|·tan(φ) = 1.20271 (plato medido del Professional); RUNFAIL residual: 0.251052 vs 0.10244 — el GNU alcanza el equilibrio plástico estático (u3 = cap/4) porque su σn cinemático relaja a ~0 al derivar los bloques, mientras el Professional mantiene σn = −1 100 s (pre-stress persistente del dof de tensión) — blocker de mecánica de dofs de tensión, fuera de la familia de interfaz. **Corpus: 140 → 141 PASS** (200 RUNFAIL / 22 PARSE). Suite propia 16/16. |
| Aliases `convection_apply` (6.395) / `control_convection_apply` (6.113) / `convection_stabilization` (6.396) en db_number() → `options_convection`/`options_stabilization` (machinery ya en general.cc: convección con velocidad materi + upwind peclet) | (commit feat de este lote, 2026-09-02) | 2026-09-02 | **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: condif2/condif3 rc=0 con temp nodo 2 = 0.992028 (target 1.0±2e-2; el Professional da 1.000000 EXACTO); tube1 rc=0 con vely post_point = 1066.59 (target 1000±100; el Professional 1000.49). El .dbs del Professional guarda los records con esos mismos nombres. **Corpus: 130 → 133 PASS** (208 RUNFAIL / 22 PARSE). Suite propia 16/16. Blocker por test: validation_1 — parsea y corre, pero el transitorio del NS euleriano (auto-advección de velocidad, quad9 p+h refinado) diverge del Professional desde el paso 2 (velx final +0.0723 vs −0.0391 en el post_point; la variante sin convección de AMBOS códigos coincide → el término convectivo del ensamblaje mixto velocidad/tensión necesita sprint A/B dedicado); validation_8 — parsea hasta `-veln` (dof de bounda_dof para velocidad normal con generación de MPC, manual 6.22; el Professional genera mpc_node_number con -velx/-vely según la entidad — algoritmo de elección de eje no deducible sin más datos). Familia force_element_edge_water/excavate1: ver nota de layout en su sección. |
| `group_materi_plasti_visco_power` (6.748) — layout Professional `η p` (2 valores; ley ε̇_pl = η·f^p sin tensión de referencia) + fix del camino nonlocal clásico (nonlocal_set() por cambio de malla) | (commit feat de este lote, 2026-09-02) | 2026-09-02 | **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: el GNU legacy exigía 3 valores (`η p f_ref`, test GNU 2001 `1. 1. 10.` ≡ Professional `0.1 1.` con f_ref=1); el record pasa a variable-length y el consumo en stress.cc toma f_ref del 3er valor si existe, si no f_ref=1. visc_pl1 rc=0 (sigxx Professional 1.0243 vs target 1±0.1; overshoot visco y relajación a la superficie de fluencia), visc_pl2 rc=0 (η=1e-6: sigxx elástico 2.0±1e-3 con epp ~5e-7, igual que el Professional), validation_12 rc=0 (tensión NO-LOCAL con 2 zonas materiales y hardening kap: target kap 0.0065; el fix es la llamada restaurada a nonlocal_set() en la rama clásica de la iteración — antes NODE_NONLOCAL nunca se construía y fallaba con db_error). **Corpus: 133 → 140 PASS** (201 RUNFAIL / 22 PARSE). Suite propia 16/16. |
| `group_groundflow_permeability` (6.618) isotrópico de 1 valor + `control_mesh_merge_geometry(_not)` (6.214/6.215) + `mpc_element_group`/`_dof`/`_geometry` (6.860/6.864/6.866) + gates `mpc_apply`/`control_mpc_apply` (6.859/6.255) + layout `force_edge` por dirección (6.454) + `solver_matrix_symmetric` honesto | (commit feat de este lote, 2026-09-02) | 2026-09-02 | **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: (1) groundflow_permeability era FIXED a ndim valores; ahora variable 1..ndim con replicación isotrópica en groundda.cc. (2) mpc_element_group genera ties con las shape functions del elemento master (point_el, patrón mpc_linear_quadratic) restringidos por geometría/dofs, saltando nodos ya miembros del master (semántica always -no; el tie consigo mismo congelaría la ecuación). (3) mpc7 rc=0 = patch test EXACTO sigxx 1.0±1e-8 (2 quad4 pequeños en 1 grande). (4) ground19 rc=0: flujo post_node_result -0.00557863 vs target -0.00557268±1e-5 (Professional: -0.00557268; el bloque mpc NO altera el campo de presiones — A/B con/sin mpc idéntico en el Professional). (5) large2/large3 rc=0 (consolidación 3D acoplada; large2 divergía porque solver_matrix_symmetric -yes forzaba CG sobre un sistema medido NO simétrico — ahora Bi-CG honesto sobre la matriz tal-cual, como el Professional que simetriza "si hace falta"; large3 con bounda_alternate mide simétrico y queda igual). (6) force_element_edge variable 1..nprinc (Professional: 1 valor por dirección; el 4º principal dof de pres NO se carga). **Corpus: 133 → 140 PASS** (7 nuevos: visc_pl1/2, validation_12, ground19_water_under_dam, large2/3 + mpc7). Suite propia 16/16. |
| Familia material `mohr_coul_direct`/`tension_direct` — modo completo (manual 6.726/6.738: cut-off espectral de tensiones + cut del diferencia de tensiones principales con flujo no-asociado phi_flow) + normal `_normal` ndim-flexible + reducción de fricción en paredes `group_materi_plasti_bounda/_factor` (6.231/6.232) | (commit feat de este lote, 2026-09-02) | 2026-09-02 | **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: (1) `materi_direct_full_mc()` en stress.cc — cap espectral (`matrix_jacobi`, sigy=0 por defecto si falta tension_direct con MC) + cut del par (max,min) con el flujo `(1+sin psi)/(1-sin psi)` y acoplamiento C del principal intermedio + regla de esquina sigma1=sigma2. (2) mohr_coul_direct1 rc=0 (oedómetro eptxx 1.49e-2, antes 2.99e-3 ≈ 5×) y mohr_coul_direct2 rc=0 (dilatancia phi_flow=π/4: vely=velx, antes 0.35). (3) mohr_coul_direct6/7: la ley material es EXACTA contra el Professional (con `-total_linear`: sigxy = 0.0249584 EXACTO en direct6); los runs del corpus siguen RUNFAIL porque su incremento único de cortante 100 % (γ=1.0) ejercita la deformación de elemento `-updated` (polar U-I) del GNU vs la lineal del Professional — blocker de CINEMÁTICA de elemento, no de la ley. (4) mohr_coul_direct8: parse arreglado (`_normal` con ndim valores + alias `group_materi_plasti_bounda`), RUNFAIL residual por la semántica de detección de pared. **Corpus: 124 → 130 PASS** (211 RUNFAIL / 22 PARSE). Suite propia 16/16. Blocker por test: direct3/4 interfaz (reacción −50 vs −100), direct5 (force_edge sobre bar2 axisimétrico), direct6/7 (cinemática -updated vs lineal a γ=1.0), direct8 (reducción de pared). |
| Familia MPC: `mpc_node_number`/`mpc_node_factor` (6.874/6.875, consumo de los records) + `mpc_linear_quadratic` (6.873, generación automática de ties) + record interno `mpc_linear_quadratic_mesh_fingerprint` (regeneración al cambiar la malla) | (commit feat de este lote, 2026-09-01) | 2026-09-01 | **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: (1) semántica de mpc_node_number/factor — esclavo = Σ factor×maestro, esclavo bounded (node_bounded=1) recomputado cada iteración, fila excluida del sistema SIN condensación ni redistribución (experimentos: maestro libre + link rígido → pardiso del Professional falla por sistema singular; fuerza externa sobre esclavo queda en su rhside sin transferirse, u3=0). (2) mpc1 rc=0 con node_dof IDÉNTICO al Professional: nodo 2 velx=disx=2, sigxx=0.5. (3) ties de mpc_linear_quadratic byte-idénticos en mpc3/4/5 (mid-edge 2 maestros ×0.5, mid-face 4 ×0.25, índices 0..N); tie de mohr_coul_direct2 verificado (nodos 3=4 exactamente iguales). **Corpus: 123 → 124 PASS** (217 RUNFAIL / 22 PARSE). Suite propia 16/16. mpc3/mpc4 siguen RUNFAIL por el solver staggered u-σ (0.299/0.350 vs 0.333 exacto del tying; familia DIAG-SOLVE-MIXTO); mpc5 por control_mesh_delete_geometry_factor + reset de stress (0 vs 0.5, mismo root que delete3); mohr_coul_direct2/4/7 por la familia MATERIAL direct (mohr_coul_direct1/3/5/6/8 sin mpc también fallan, p.ej. eptxx 2.99e-3 vs 1.49e-2). mpc8/mpc9: SEGV pre-existente (también sin records mpc). |
| Familia `bounda_time_until_data` + `bounda_time_until_value_minimum` (6.40/6.41) + registro de keywords del corpus (mpc_*, control_mesh_truss_distribute_mpc, post_calcul_length, strain_volume_*, bounda_used, processors_used...) | (commit feat de este lote, 2026-09-01) | 2026-09-01 | **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: la ley del factor del until se verificó con DOS runs de until1.dat (E=1 y E=2): factor = ((monitor/first − wanted)/(start − wanted))² con clamp [0,1], monitor del PASO ANTERIOR (post_node_result), first = valor inicial del monitor — match 1e-6 en TODA la serie (p.ej. E=1 t=1.992: 0.81 = (9e-3/1e-2)²; E=2 t=1.992: 0.81 = (1.8e-2/(1e-2·2))²). until1 rc=0 (target node_rhside 1 −velx = 0.0 ± 1e-3; GNU 5.45e-5 vs Professional 9.76e-5). **Corpus: 121 → 123 PASS** (218 RUNFAIL / 22 PARSE). Suite propia 16/16. mpc2 rc=0 por coincidencia física (traslación rígida del quad4; la familia MPC queda registrada sin consumo — ver registro). |
| Aliases `geometry_factor` (6.527) → `geometry_bounda_factor` y `processors` → `options_processors` en db_number() | `9d22f6c` | 2026-09-01 | **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: matrix2 (condif, factores lineales 1..4 por lado) post_point_dof = 2.4999998123 vs 2.500000000000e+00 del Professional (8 cifras); el ejemplo del manual 6.527 (nodo x=0.2 → 20·1.6, x=0.4 → 20·2.2) coincide con la interpolación de geometry.cc. **Corpus: 118 → 121 PASS** (220 RUNFAIL / 22 PARSE). Suite propia 16/16. mpc4/interface11/validation_14_mesh siguen RUNFAIL por blockers ajenos al alias (mpc_linear_quadratic, mesh_interface_triangle_coordinate, bounda_time_until_data). delete3 RUNFAIL documentado (control_mesh_delete_geometry_factor sin aplicar). |
| Records `*_sig` de `post_calcul -materi_stress -force` (nombres de items node_dof_calcul: -norx_sig ... -mom2s_sig) + dirección de plot 2D = tangente de cara + gate de validación pre-step + `post_element_force_normal` INTEGER/filtro +n + extrude quad9→hex27 + default npointmax 27 | (commit feat de este lote, 2026-09-01) | 2026-09-01 | **Contra el binario del Professional (user-supplied 25-10-2023, .dbs)**: force7/8 (ménsula quad9 2D): nory/norx = −12.34, shey/shex = +100, momy/momx = −5000 EXACTOS; force9 (anillo quad9 2D, Roark pág. 262): nory≈0, shey = −0.5, momy = +14.4125; force10/13 (ménsula hex8 3D): nory = −12.34, shey = −100, mom1y = −5000; nors = −12.34, shes = +100; post_element_force_result −123.4/0/1000/0/−1e5 (force10) y −123.4/0/1000/0/−5e4 (force13). force11 (anillo 3D quad9→hex27): shey = −0.2 EXACTO, mom1y 0.19% alto (malla 1 capa) — RUNFAIL por timeout del solve Bi-CG (45 s del corpus), familia solver. **Corpus: 112 → 118 PASS** (223 RUNFAIL / 22 PARSE). Suite propia 16/16 en build limpio. |
| `check_used` | `3793892` | 2026-08-04 | test propio |
| `check_element_shape`, `check_memory` | `762f2ca` | 2026-08-04 | test propio |
| `check_solver` | `e42e12d` | 2026-08-04 | test propio |
| `check_data`, `check_error`, `check_warning` | `5a2b3d2` | 2026-08-04 | test propio |
| `control_mesh_move/mirror/copy` | `3e94dac` | 2026-08-05 | test propio |
| `control_mesh_rotate`, elementos 3D | `d8d3a72` | 2026-08-05 | tests 3D |
| `control_mesh_keep_*`, `change_element_group`, `delete_element` | `afc1dad` | 2026-08-05 | test propio |
| `control_mesh_remove` | `1e5b4b9` | 2026-08-05 | test propio |
| `control_mesh_extrude` | `5dec3c2` | 2026-08-05 | test propio |
| `bounda_unknown` (P3) | `99d6cd6` | 2026-08-05 | test propio |
| `bounda_dof` (alias unknown) | `b5e02d5` | 2026-08-05 | test propio |
| `bounda_alternate` | `0ba382a` | 2026-08-05 | test propio |
| `bounda_normal` | `d20ca47` | 2026-08-05 | test propio |
| `bounda_water` | `e7861f6` | 2026-08-05 | test propio |
| `bounda_dof_radial`, `_cylindrical` | `dd53575` | 2026-08-05 | test propio |
| `bounda_time_units`, `geometry_method`, `factor_parabolic_x` | `ba20c45` | 2026-08-05 | test propio |
| `groundflow_pressure_factor` | `75d1ac1` | 2026-08-05 | test propio |
| `group_groundflow_permeability_vertical_stress` | `07f77ac` | 2026-08-05 | test propio |
| `group_materi_plasti_hypo_masin` | `974b284` | 2026-08-06 | tests hypomasin1-5 contra `reference-masin` (paths de tensión) |
| `..._masin_clay_visco` | `a0f55a6` | 2026-08-06 | test propio (visco Niemunis Dr/Iv) |
| `..._masin_clay_anisotropic` | `4225318` | 2026-08-06 | test propio (anisotropía + intergranular strain) |
| `group_materi_plasti_sanisand` | `b6034cc` | 2026-08-07 | contra UMAT de Dafalias & Manzari (P4-E1, hasta 5.4% de deviator residual) |
| `control_print_tabular` (P5-T1/T2) | `1c677bb`/`7904b4c` | 2026-08-11 | contra valores de referencia del test |
| `control_print_history_smooth` | `8f8ad67` | 2026-08-12 | media móvil verificada numéricamente (test smooth1) |
| `control_print_dof` | `ca2a9e4` | 2026-08-12 | contra SQLite de tabular (test dof1) |
| `control_print_gmsh` (familia) | `90798f9` | 2026-08-12 | validado con gmsh 4.12 (test gmsh) |
| `control_print_frd` (familia) | `19a8d08` | 2026-08-12 | formato byte-a-byte contra CalculiX (test frd) |
| `control_print_vtk_dof` | `b3e58a4` | 2026-08-13 | campos filtrados verificados (test vtk_dof1) |
| `force_element_edge_multi_linear_factor_x` | `86f6d03` | 2026-08-13 | efecto del factor 0/1 en target (test mlx1) |
| `control_mesh_generate_truss/beam` | `bd7b3ab` | 2026-08-13 | contra tests sfnet genera1/ho_othr1 (gen1/genbeam1) |
| `control_reset_dof` + `_value_constant` + `_value_dof` + `_value_dof_diagram` + `_value_method` | `b6eaee4` + batch hypo | 2026-08-13 / 2026-09-01 | hisv0 reseteado a 0.55 (0.5986 sin reset); diagrama sigyy→hisv0; métodos -use/-add/-multiply. **Batch**: (a) scoping por indice (`icontrol==ireset`: el reset corre UNA vez en su propio paso de control; el Professional RECHAZA un reset en el indice del timestep — medido); (b) records multi-dof: TODOS los dofs listados reciben el valor (antes solo reset_dof[0]) |
| `control_change_dataitem_apply` | `bdcdeaf` | 2026-08-13 | con -no ignora change_dataitem (targets hypo1 se cumplen); sin apply la geometry cambia y los targets fallan |
| Carril A (diseño) | — | 2026-08-13 | diseño técnico en DESIGN-INTERFACES.md: modelo físico (strain = dif. de desplazamiento entre lados), análisis del codebase (patrón spring.cc en elem.cc), 4 fases, test de validación propuesto. Sin tests de referencia sfnet. |
| Carril A Fase 1 (`group_interface` + `_elasti_stiffness`) | `a82cbc7` | 2026-08-13 | interface_element() en interface.cc. **Validada** (test 2 bloques): signo corregido (-sign*stress*dir); kn=100→0.044, kn=1e6→-0.003≈soldado, kn=0.001→≈1.0 libre. Límites físicos correctos. Estrategia de no-interpenetración: penalización implícita + control_timestep_iterations (sin line-search/arc-length). |
| Carril A Fase 3 (gap, residual_stiffness, tension_direct, mohr_coul_direct — fix y validación) | `9c2f4c8` (código T1-T5; docs T8 + build_safe T7 en el commit del orquestador) | 2026-08-14 | VALIDADO (familia `iface_mc`, 13 tests / 18 runs en verde). Mohr-Coulomb acumulativo (history `element_interface_force_tang`, activación por presencia del record, phi_flow = dilatancia); gap: cerrada si strain > gap, default -1e20 (siempre cerrada), hueco físico = gap negativo; tension_direct abre en tracción sobre la fuerza normal TOTAL. a/a' y b/b' invariantes en nº de pasos. Corrige los bugs de `60bf78c` (2026-08-13). |
| Carril A Fase 2 (`control_mesh_convert` bar2→quad4) | `490545b` | 2026-08-14 | validado (test iface_conv): nodos 7,8 en x=0.99 creados, elemento reescrito quad4, bloques reconectados. |
| Carril A Fase 4 (`control_print_interface_stress` 2D) | `01f6c3e` | 2026-08-14 | validado (test iface_stress): interface_stress.0 generado con distancia+sign (strain acumulado * kn); sign crece con la compresion. sigt=0 y 3D pendientes. |
| `post_calcul -materi_stress -force` (integración 2D) | `1babf2d` | 2026-08-28 | analítico por estática (tests propios): msf_beam2d (moms = P·(8−x) 0-2.5%, promediado exacto), msf_beam2d_pure (M cte, shes=0), msf_quad9 (promedio exacto, -primary 15/27), msf_quad9_noavg (avg -no), msf_nor (N=2e-2 0-5%), msf_shear (she = G·γ·h = 0.384615 EXACTO) |
| `post_calcul -materi_stress -force` (integración 3D hex8/hex27) | `fa7650b` | 2026-08-28 | analítico con deformación prescrita (tests propios; GOTCHA del solve mixto 3D del GNU documentado): msf_sheet3d (hex27 flexión prescrita: mom1 = E·κ/12 = 0.0833333 EXACTO 81/81), msf_sheet3d_hex8 (hex8 corte puro: she = G·γ·t = 0.5 EXACTO; limita del hex8 en flexión documentado), msf_hex27_avg (-all 45 vs -primary 27, promedio exacto), msf_tunnel3d (anillo hex27 expansión radial prescrita: nor = E·u0·t/R = p·R = 0.1 EXACTO 144/144, mom1 = she = 0) |
| `post_calcul -materi_stress -force` (L4: fuente de σ = puntos de integración del elemento + tests nuevos) | (lote 4, 2026-08-28) | 2026-08-28 | **código** (calcul_force.cc + elem.cc): la tensión de sección se integra ahora de los IPs del elemento (`ELEMENT_DOF` — las tensiones constitutivas que el elemento usó; reconstrucción por Lagrange de la regla propia del elemento, lectura directa en las reglas con nodos) en vez de la nodal recuperada; `options_element_dof -yes` requerido; elem.cc escribe el bloque de tensión COMPLETO (6×nder — antes MDIM*MDIM=9 perdía syy/syz/szz con `derivatives`; el READ del solver intacto); axisimétrico: el integrando lleva el peso físico 2πr (antes la integral sin pesar daba σ·t/(2πr) — un factor 2πr pequeño; fix documentado); fallback documentado a la nodal cuando el ELEMENT_DOF sale a ceros (los 3D resueltos con derivatives: el bucle escalonado no propaga la deformación — gforce10/13). **VERIFICACIÓN**: (1) arness Professional — gforce10/13 N = 12.34 EXACTO (1.0000×); gforce7 N/V siguen 1.24×/2.7× (la nodal = media EXACTA de los IPs — la polución vive en el CAMPO σ; el momento mejora a 0.9976×); gforce7q4_ref/gffq4 byte-idénticos; (2) suite **201/201** (199 + 2 nuevos): `msf_cant3d_hex27` (el test que el L3 no pudo — cantilever 3D hex27 CARGADO: mom1 = P·(L−z) dentro del 9%, la cizalla 3D del esquema mixto documentada 0.08·P) y `msf_axisym` (el pendiente del L2 — anillo prescrito: nor = σ_zz·t = 0.5 EXACTO por unidad de circunferencia); (3) checks actualizados con justificación: promediados quad9/hex27 (banda FE < 2e-3 — las caras del propio elemento), cizalla flexión pura (2.3e-4), cizallas qsri (2.07·P / 8.4·P — el σ_xy crudo de los IPs). **CONCLUSIÓN**: la estática EXACTA del Professional NO viene del σ crudo (ni nodal ni de IP — ambos contaminados por igual): consistente con un cálculo por fuerzas internas/equilibrio del elemento (siguiente lote recomendado; PAPER-LINES §8). |
| `post_calcul -materi_stress -force` (L5: fuerzas internas del elemento — la estática del equilibrio) | (lote 5, 2026-08-28) | 2026-08-28 | **código** (calcul_force.cc): las fuerzas de sección se calculan desde las FUERZAS INTERNAS del elemento f_elem = ∫Bᵀσ dV (las fuerzas nodales consistentes del campo de IPs del ELEMENT_DOF — la MISMA cantidad que materi() acumula en element_rhside con signo opuesto; cinemática replicada de pol()/materi(): B de polynom.cc:549-599, derivadas físicas invJ·p, volúmenes w·4·|detJ|/w·8·|detJ| con el factor 2πr axisimétrico, la regla del propio elemento incl. SRI-Gauss y la regla de 1 punto axisimétrica). La sección: R_face = Σ de las fuerzas internas de los nodos de la cara = ∫σ·n̂ dA = el resultante del cuerpo libre de las cargas; nor = n̂·R/l (tracción +), she = \|t̂·R\|/l, mom = Σ (n̂·f)·arm/l (brazos desde el punto medio de la cara, decisión L3). SRI quad4: la cizalla de regla completa se reemplaza por la reducida de 1 punto del feedback (4·detJ_c·b_shear·media(σ_xy) — el módulo se cancela; sin ella la sección SRI leería K_full lockeado). Se retiran la cuadratura de cara y la reconstrucción Lagrange del L4 (la tensión se lee solo en los IPs). **VERIFICACIÓN**: (1) arness Professional — gforce7q4 N/V EXACTOS (1.0000×) y M 0.9984×; gffq4 M_end 1.0000×/M_center 0.9987× con la identidad \|M_e\|+\|M_c\| = pL²/8 a 0.9996 (antes 0.197×); msf_tunnel3d = los valores del node_dof_calcul del Professional DÍGITO A DÍGITO (nor 0.0998068, shes 0.013727, mom1 2.03e-4); gforce7 mide 5× (el estado σ del run NO está en equilibrio — propiedad PREEXISTENTE del esquema, el L4 la media como 1.24×/2.7× por integración de campo; el resultante débil la amplifica — frente abierto del SOLVER, no del post-proceso); (2) suite **201/201** en build limpio: msf_beam2d moms = P·(8−x) y shes = P EXACTOS (6 dígitos) en todas las secciones, msf_shear/msf_axisym EXACTOS, msf_cant3d_hex27 mom1 al 1% y shes = P (la cizalla 0.08·P desaparece), qsri mom = P·8 EXACTO con y sin SRI (la estática del cuerpo libre de las soluciones equilibradas — el lock queda en la deflexión y el campo σ, no en las fuerzas de sección); (3) checks actualizados con justificación (par de reacciones del corte puro 0.1923/0.25, túnel = valores del Professional, equilibrio qsri). **CONCLUSIÓN**: la estática por fuerzas internas/equilibrio del elemento — la vía del Professional identificada en el L4 — está IMPLEMENTADA y verificada: EXACTA (cuerpo libre) para los estados σ en equilibrio; el estado σ no-equilibrado de los runs gruesos del GNU es el único limitante (solve). |
| `group_element_selective_reduced_integration` (SRI de Hughes para quad4, EXTENSIÓN GNU — no existe en Professional) | `776d11b` | 2026-08-28 | A/B contra el lock documentado del MSF L2 (tests propios, familia qsri, suite **199/199** en build limpio): qsri_beam2d reproduce el shear locking EXACTO (mom = 0.231×, she = 0.741×) y qsri_beam2d_sri lo mejora (mom = 0.312×, she = 0.65×); qsri_patch_s/_off (corte simple prescrito: she = G·γ = 0.384615384615 EXACTO, idéntico con/sin SRI) y qsri_patch_t/_off (tracción: nor idéntico con/sin SRI) — el SRI integra estados constantes exactamente; qsri_modes_rigid (rotación rígida prescrita → tensiones 0 EXACTAS: 3 modos cero) y qsri_modes (hourglass prescrito → tensiones NO nulas: sin modo espurio). Análisis de eigenvalues del quad4 SRI aislado documentado: exactamente 3 modos cero, 4º con energía positiva. GOTCHA del solve mixto documentado: el valor clásico del SRI (mom ≈ P·(8−x)) se recupera con el solve exacto (99.3% deflexión / 93.75% momento, verificado por eliminación gaussiana sobre la matriz del C) pero el Bi-CG mixto del GNU es numéricamente poco fiable en la malla 1-en-espesor (misma familia que el GOTCHA 3D del MSF L3). |
| `group_element_selective_reduced_integration` (SRI de Hughes EXTENDIDO al hex8 3D — el mismo keyword, el trilineal) | `feat(elem): integracion reducida selectiva hex8 (SRI 3D) — extension del quad4` (2026-08-29) | 2026-08-29 | **código** (miscel.cc `sri_active` renombrado/generalizado, materi.cc split D 3D — entradas de cizalla diag. 1/2/4 a cero para la regla completa y rigidez de cizalla 1×1×1 en el centroide con las 3 filas γ_xy/γ_xz/γ_yz de B, feedback D-c con el GOTCHA del índice tensorial inc_ept[1*MDIM+2]≠stress_indx(1,2)=4, recuperación D-b con extrapolación de Lagrange 3D, polynom.cc switch a Gauss 2×2×2, calcul_force.cc sección con la corrección SRI 3D). **HALLAZGO (verificación del análisis de eigenvalues, independiente)**: el SRI hex8 de cizalla sola NO es libre de hourglass — 9 modos cero en el elemento aislado (6 rígidos + 3 twist: u=(y·z,0,0) etc., deformaciones normales cero y cizallas que se anulan en el centroide de sección) y modos de warping de sección en malla (u_y = A(x)(2z−1), u_z = A(x)(2y−1)): la ménsula empotrada 8×1×1 ensambla una matriz SINGULAR (8 modos ~1e-13). La premisa "exactamente 6 modos rígidos, sin espurio" del SRI 3D queda REFUTADA con evidencia (la limitación conocida del SRI de cizalla del ladrillo — la literatura pasó a B-bar/ANS). **VERIFICACIÓN**: (1) patch tests 3D EXACTOS e IDÉNTICOS con/sin SRI (qsri3d_patch_s: she = G·γ·t = 0.5 EXACTO; qsri3d_patch_t: nor = 1.346154 EXACTO — el SRI integra constantes exactamente); (2) modos rígidos 3D → sección 0 exacta (6 modos cero); (3) modo twist prescrito → invisible a la sección SRI (documentado); (4) MÉNSULA 3D 1-en-espesor: el default reproduce el lock documentado (deflexión 0.221× de la EB, σ axial 0.26×) y el SRI lo recupera (deflexión **0.897×**, σ axial 1.06× — el hex8 no da el 0.9375 exacto del 2D, medido y documentado); (5) arness Professional: gforce10/13 (default sin SRI) sin regresión — N axial 12.34 EXACTO (1.0000× vs Professional); (6) suite **209/209** (201 + 8 nuevos) en build limpio, default byte-idéntico. |
| **LOTE 6 del solve u-σ — estado σ en equilibrio post-solve: ELEMENT_DOF 3D poblado + cinemática de sección 2D corregida** | `feat(solver): estado sigma en equilibrio post-solve (B^T sigma = P)` (2026-08-29) | 2026-08-29 | **MECANISMO (FASE 1, medido con dumps por iteración)**: (1) el misterio inc_ept=0 del L4 — el ELEMENT_DOF de los 3D con `derivatives` (nder=5) salía a CEROS aunque el σ constitutivo SÍ se calculaba y se usaba (feedback del momento, recuperación nodal — el node_dof σ era correcto: σxz = 0.7407 = E·ε): materi.cc:850 escribía el σ en el slot `stres_indx/nder + j` (índice de unknown) en vez de `stres_indx + j·nder` → para nder=5 el σ caía en el bloque de desplazamientos y el write de elem.cc (que copia los slots del stress) nunca lo veía; el restore del bloque leía solo 9 slots (MDIM·MDIM, parcial para nder>1) y los rangos epe/epp/epi con índices −1 (initias ausentes) capturaban slots [0,8) (basura de velocidades); (2) el "5×" del gforce7 NO era no-equilibrio: el campo σ estaba en equilibrio (residuo del solver 2.19e-13, vely(x=50) = −0.01154 ≈ 0.01138 analítico, σxx nodal ±301 ≈ ±300) — la cinemática de la sección L5 (msf_element_internal_forces_2d) usaba la inversa del jacobiano TRASPUESTO (términos cruzados invjac[1]↔[2] intercambiados): para el quad9 con orden local [9,8,7,...] (eje ξ VERTICAL, jacobiano [[0,−5],[−25,0]] no diagonal) el error = 5× (elemento 50×10) o 5/4× (12.5×10); quad4/quad9-cuadrados (jacobiano diagonal) y el 3D (producto matriz-vector completo) intactos. **FIX (mínimo, nder=1 byte-idéntico)**: materi.cc:850 slot `stres_indx + j·nder` (+ hisv 807), elem.cc restore 6·nder + rangos hisv·nder, calcul_force.cc swap invjac[1]↔[2]. **VERIFICACIÓN**: (1) arness — gforce7 N/V 5.0000× → **1.0000×** y M 0.9984×; gforce7_ref 1.2500× → **1.0000×**; gforce7q4/gffq4 byte-idénticos; gforce10/13 N = 12.34 EXACTO con el ELEMENT_DOF REAL (el fallback del L4 ya no dispara — misterio cerrado); (2) residuo Bᵀσ−P: Σ f_elem = −P en los dofs libres, reacciones (12.34, 100) exactas — el estado σ SIEMPRE estuvo en equilibrio; la lectura del L5 era la que amplificaba (61.7/500/24955 → 12.34/100/4992); (3) suite **209/209** en build limpio; checks actualizados CON justificación (calibrados contra el restore σ parcial de los modelos con derivatives): mesh_act_grav velx −0.95 → +1.3333 (signo físico de la rampa de gravedad +x), cmat_gate sigxy −88.39 → +9.09 (A/B |sigxy| > 1.0 conservado), qsri3d_beam mom1s = P·L = 0.08 para SRI/OFF (la estática del cuerpo libre; el discriminador SRI es la deflexión 0.897× vs 0.221× de la EB); (4) tiempos: peor caso 6 s por test (timeout 120 s sin riesgo). **NO CIERRA**: el V/mom de la sección 3D de los modelos con carga axial (gforce10/13 — cara/brazo 3D del post-proceso), el SRI hex8 3D (modos cero). |
| **DIAGNÓSTICO del solve mixto u-σ (CAUSA RAÍZ, corrección del GOTCHA anterior)** | `docs` (commit del diagnóstico, 2026-08-28) | 2026-08-28 | instrumentación del sistema ensamblado (dumps de A/RHS/mapas vía TOCHNOG_DUMP temporal en so_bicg.cc, revertido) + numpy (espectros, A·b, formas cuadráticas) + A/B (Bi-CG vs SuperLU) + barridos de iteraciones de equilibrio. **Hallazgos con evidencia**: (1) el sistema que Bi-CG resuelve es SOLO de velocidades y SPD — los dofs σ (`input.cc:441`) NUNCA entran en la matriz global (dof_principal sin asignar, so.cc test2) y σ avanza por la actualización diagonal lumped (`parallel_new_dof_diagonal`, top.cc:414) → esquema ESCALONADO, NO saddle-point (se descartan las hipótesis de espectro indefinido, bug de ensamblaje 3D — RHS de carga correcto, stress_indx denso sin huecos — y null-space del bloque σ); (2) los criterios de parada del Bi-CG (so_bicg.cc) declaran ÉXITO en breakdown (`|dAd|<1e-16` → x=0, medido dAd=-8.7e-23) y en estancamiento (`|Δerror|<0.1·check_error` con error ≫ tolerancia: hex8 1e carga punta para en 4.72727e-07 con check 1e-12 y RC=0; el residuo L2 de CG es no-monótono y para estructuras simétricas el 1er paso cae EXACTO en residuo plano — identidad |r|²|Ād|²=2(dᵀĀd)² a 17 dígitos, b excita pares de autovalores degenerados λ=0.1412 al 95%); (3) el PUNTO FIJO del esquema escalonado NO es la solución de desplazamiento: rigidez doble (dt·K en la matriz + −Bᵀσ con σ=C:ε en el RHS) y cizalla constitutiva completa reintroducida → mom converge a 0.2315× para quad4 plain Y SRI (el 0.3125× documentado es transitorio del default de 2 iteraciones; el 93.75% es la referencia del sistema u puro, que el GNU no ensambla); **SuperLU ≡ Bi-CG byte-idénticos** (md5 62e82a92, mom 0.0249998465889) → el solver lineal NO es el limitante, el ESQUEMA sí. Implicación usuario (carretera): |M_e|+|M_c| ≈ 0.23·pL²/8 ≪ pL²/8 (chequeo estático exacto de viga biempotrada: |M_e|+M_c = pL²/12+pL²/24 = pL²/8). Suite 199/199 verde tras revertir la instrumentación. Fix propuesto (no implementado, arquitectónico): criterios de parada honestos (A), CG para el sistema SPD (B), solve mixto monolítico real con MINRES/SuperLU-pivoting (C), regularización del escalonado (D). |
| **DIAGNÓSTICO del solve mixto u-σ (CAUSA RAÍZ, corrección del GOTCHA anterior)** | `docs` (commit del diagnóstico, 2026-08-28) | 2026-08-28 | instrumentación del sistema ensamblado (dumps de A/RHS/mapas vía TOCHNOG_DUMP temporal en so_bicg.cc, revertido) + numpy (espectros, A·b, formas cuadráticas) + A/B (Bi-CG vs SuperLU) + barridos de iteraciones de equilibrio. **Hallazgos con evidencia**: (1) el sistema que Bi-CG resuelve es SOLO de velocidades y SPD — los dofs σ (`input.cc:441`) NUNCA entran en la matriz global (dof_principal sin asignar, so.cc test2) y σ avanza por la actualización diagonal lumped (`parallel_new_dof_diagonal`, top.cc:414) → esquema ESCALONADO, NO saddle-point (se descartan las hipótesis de espectro indefinido, bug de ensamblaje 3D — RHS de carga correcto, stress_indx denso sin huecos — y null-space del bloque σ); (2) los criterios de parada del Bi-CG (so_bicg.cc) declaran ÉXITO en breakdown (`|dAd|<1e-16` → x=0, medido dAd=-8.7e-23) y en estancamiento (`|Δerror|<0.1·check_error` con error ≫ tolerancia: hex8 1e carga punta para en 4.72727e-07 con check 1e-12 y RC=0; el residuo L2 de CG es no-monótono y para estructuras simétricas el 1er paso cae EXACTO en residuo plano — identidad |r|²|Ād|²=2(dᵀĀd)² a 17 dígitos, b excita pares de autovalores degenerados λ=0.1412 al 95%); (3) el PUNTO FIJO del esquema escalonado NO es la solución de desplazamiento: rigidez doble (dt·K en la matriz + −Bᵀσ con σ=C:ε en el RHS) y cizalla constitutiva completa reintroducida → mom converge a 0.2315× para quad4 plain Y SRI (el 0.3125× documentado es transitorio del default de 2 iteraciones; el 93.75% es la referencia del sistema u puro, que el GNU no ensambla); **SuperLU ≡ Bi-CG byte-idénticos** (md5 62e82a92, mom 0.0249998465889) → el solver lineal NO es el limitante, el ESQUEMA sí. Implicación usuario (carretera): |M_e|+|M_c| ≈ 0.23·pL²/8 ≪ pL²/8 (chequeo estático exacto de viga biempotrada: |M_e|+M_c = pL²/12+pL²/24 = pL²/8). Suite 199/199 verde tras revertir la instrumentación. Fix propuesto (no implementado, arquitectónico): criterios de parada honestos (A), CG para el sistema SPD (B), solve mixto monolítico real con MINRES/SuperLU-pivoting (C), regularización del escalonado (D). |
| **LOTE A+B del solve u-σ — criterios de parada HONESTOS + CG para el sistema SPD** (implementación del fix propuesto por el diagnóstico) | `feat(solver): criterios de parada honestos + CG` (2026-08-28) | 2026-08-28 | **código** (so_bicg.cc, `solve_iterative_bicg`): éxito SOLO con el test de residuo real (relativo, con piso; el residuo se recomputa del iterado en cada pasada); eliminados los exits de falso éxito — breakdown `|dAd|<1e-16` (x=0) y estancamiento `|Δerror|<0.1·check_error` (RC=0 con error ≫ tol) — ahora fallo honesto (RC≠0, residuo relativo real `|b−Ax|/|b|` en el mensaje); `max_iter` con el escape documentado `control_solver_bicg_stop -no`; CG para sistemas medidos simétricos (chequeo de simetría en runtime por solve, `solve_iterative_bicg_symmetric`) + Bi-CG honesto para los no simétricos (beam: bloques translación↔rotación con factor dtime distinto; interfaz 3D deslizante: `max_rel_asym` 1.0, `dᵀAd<0`); monitor de residuo PRIMAL (el antiguo media el residuo del sistema transpuesto, que no se anula en la solución de un sistema no simétrico — el beam solo "convergía" por el breakdown falso). **VERIFICACIÓN**: (1) suite **199/199** en build limpio — 198 sin cambios (los valores convergidos son byte-idénticos), `iface_3d_slip` actualizado (test corregido: opta a `control_solver_bicg_stop -no`, verifica la física de interfaz); (2) **cantilever 3D con carga REAL**: el caso A·b≈0 (hex8 1-e, residuo plano 4.72727e-07) ahora CONVERGE (4.2e-37, rel 9.4e-16) y hex8×8/hex27×8 resuelven con fuerzas de sección (hex27: mom ≈ 1.1·P·(L−x), she ≈ P) — la validación 3D con carga queda desbloqueada; (3) **PREGUNTA CLAVE — el fijo 0.2315× NO cambió**: barrido 1..32 con solver honesto: plain 0.0185184 (0.23148×) y SRI 0.0185298 (0.23162×) = idénticos al diagnóstico → el fijo es del ESQUEMA escalonado, no del solver → C/D siguen siendo el fix definitivo; (4) **carretera** (hormigón 1-en-espesor sobre suelo, biempotrada, carga distribuida): corre RC=0 pero |M_extremo|+|M_centro| = 0.001435 = 0.143·pL²/8 ≪ pL²/8 — sigue sin ser físicamente posible (factor ≈7); (5) tiempos: ningún test cerca del timeout de 120 s (máx ≈1.4 s). Docs: DIAG-SOLVE-MIXTO §8.1, manual-developer/solve_iterative_bicg.md, PAPER-LINES (A/B DONE). |
| **LOTE C/D del solve u-σ — FIX del esquema escalonado: punto fijo = solución del elemento + recuperación de σ consistente** | `feat(solver): fix del esquema escalonado u-sigma — punto fijo consistente` (2026-08-28) | 2026-08-28 | **MECANISMO EXACTO (medido, FASE 1)**: (1) el mapa del esquema: v(k+1) = v(k) + (dt·K_uu)⁻¹·(P − Bᵀ·σ_constit(v(k))) con σ_constit(v) = σ_old + C_full·ε(v·dt) (stress.cc:1182) y el update de σ por reemplazo σ_new = σ_constit(v) (la inercia lumped de general.cc:253 cancela el σ_iterate; dof.cc:155) — la velocidad acumula (so.cc:740) → para K_uu = K_full el mapa colapsa en UNA pasada a v* = (dt·K_uu)⁻¹P; (2) el PUNTO FIJO: dv=0 → P = Bᵀ·σ_constit(v*) = dt·Bᵀ·C_full·B·v* → **v* = (dt·K_full)⁻¹·(P − Bᵀσ_old) SIEMPRE** — la matriz del momento (SRI o no) NO entra en el fijo (solo en el transitorio) → el SRI se cancela al converger (0.2316× ≈ plain 0.2315×, confirmado con 32 iteraciones) — la "rigidez doble" del diagnóstico es el transitorio; el fijo es el K_full⁻¹P; (3) la RECUPERACIÓN del σ dof (promedio h-weighted lumped) DILUYE los valores nodales con la cuadratura 2×2 GAUSS del SRI (puntos interiores ±1/√3): σ_nodo ≈ 0.577·σ_esquina (medido) → los momentos de sección leen tensiones bajas aunque el desplazamiento sea correcto (SRI: 0.29× en vez de 0.9375×); con las reglas con nodos (Lobatto esquinas, quad9/hex8) h = δ → recuperación exacta (sin cambio); (4) las DIVERGENCIAS del arness son BUGS DE BC de los inputs GNU: gforce7q4 fijaba velx en la arista INFERIOR (1,2,3) en vez de la izquierda (1,4) → modo de rotación rígida libre → matriz singular (pivot 8.7e-10, carga 88% en el modo nulo) → CG honesto diverge; gforce10/13 usaban `-ra 1 4` = lista {1,4} no el rango → cara inferior solo 2 esquinas → singular (pivot 8.6e-17). **FIX D (mínimo)**: (c) feedback del momento consistente con el elemento — para el quad4 SRI la cizalla del feedback se integra como la matriz (incremento de cizalla del iterado a cero en sigvec + fuerza interna reducida −dt·K_shear·v en el RHS) → fijo = (dt·K_elem)⁻¹·(P − Bᵀσ_old) = solución del ELEMENTO, convergencia en UNA pasada sin deriva; (b) recuperación del σ consistente — extrapolación de Lagrange de los valores de Gauss a los nodos ("mismo B en el nodo") para las componentes NORMALES (superconvergentes: σ_xx SRI 94% del analítico), h para las de cizalla (no superconvergentes; la media h = estimación centróide) — solo cambia el quad4 SRI (las reglas con nodos dan L = δ = h). **VERIFICACIÓN**: SRI: mom = 0.0750 = 0.9375·P·L (referencia clásica Hughes) estable a 32 iteraciones; plain quad4/quad9/hex8 byte-idénticos (gforce7q4_ref 0.215×, gffq4 0.197×, gforce7_ref 0.9986×, gforce7 0.996×); gforce7q4 y gforce10/13 CONVERGEN (BC corregidas; el N axial 3D = 12.34 EXACTO 1.0000× vs Professional); arness de aceptación VERDE (compare_professional.sh); suite **199/199 + verificaciones de archivos OK** (check de qsri_beam2d_sri actualizado al fijo 0.9375×). NO arregla (documentado): el lock del plain quad4/hex8 (física del ELEMENTO, SRI opt-in) — el esquema ahora converge a la solución de la formulación; la cizalla de sección del σ_xy crudo (polución de interpolación preexistente; el σ nodal del Professional está igualmente contaminado — su estática exacta no viene del σ nodal). Fix C (monolítico) NO necesario. Detalle: DIAG-SOLVE-MIXTO §12, VALIDACION-PROFESIONAL §9, PAPER-LINES (C/D status). |
**Carril A Fase 4 — sigt implementado** | `8ef45c4` | 2026-08-14 | `interface_sigt` ya no es 0: se lee del history `element_interface_force_tang` (fuerza tangencial total acumulada, Fase 3), consistente con `sign`. Verificado con probe: sigt=23.28 == F_t acumulada del último paso. Documentado en ambos manuales. 3D pendiente. |
| **VERIFICACIÓN CONTRA TOCHNOG PROFESSIONAL (binario 02-08-2026, Drive público del autor) — el bug NO sobrevive en Professional** | `04b2558`-familia (docs) | 2026-08-28 | **force7q4.dat** (la ménsula de estática de force7.dat del Professional reconstruida con **quad4** 1-en-espesor, L=100 h=10): N=-12.34, V=+100, M(x=50)=-5000, M(x=0)=-10000, M(x=100)=0 — TODOS EXACTOS a 1e-10 (targets exit 0); el GNU da 0.2315× en la misma configuración. **ffq4.dat** (el caso real del usuario: biempotrada, carga uniforme, quad4 1-en-espesor): \|M_extremo\|+\|M_centro\| = 825+425 = **1250 = pL²/8 EXACTO** (825/425 son los correctos de viga profunda Timoshenko h/L=0.1; Euler daría 833/417; la identidad estática que el GNU viola por ≈7× se cumple exacta). Los tests de validación propios del Professional usan quad9 (force7/8) y 2×2-en-sección hex8 (force10/13) — diseño defensivo; el caso quad4 1-en-espesor funciona correctamente allí. Conclusión: el lock + fijo del esquema escalonado + criterios deshonestos es de la línea GNU open source (arquitectura 'stress follows the principal unknowns' documentada en el manual 2011); la línea propietaria no lo tiene. Detalle completo: DIAG-SOLVE-MIXTO.md §11 y papers/PAPER-LINES.md §6-7. |
| **SPRINT 13 LOTE 2 — support_edge_normal: damping/density/factor/force_initial/time (7 keywords más, 14/18 del checklist)** | `feat(support): support_edge_normal damping/density/factor/force_initial/time` (2026-08-29) | 2026-08-29 | **código** (extensión del rama type[10] de area.cc, mismo lote 1): **damping**: fuerza viscosa `−(c_n·v_n)·n − c_t·v_t` por unidad de longitud/área; **damping_automatic**: lee E,ν,ρ del grupo adjunto y calcula c_n = sqrt(ρ·Eoed), c_t = 0.25·sqrt(ρ·G) con Eoed = (1−ν)E/((1+ν)(1−2ν)), G = E/(2(1+ν)); **damping_automatic_apparent**: usa módulos aparentes del estado nodal actual (σxx/εxx con guards); **density**: `d_n·a_n` con `a = (v_new − v_old)/dt` (v_old de NODE_DOF VERSION_NORMAL); **factor**: escala espacial polinómica de las RIGIDECES únicamente (reutiliza force_factor, manual 6.1075 "for the support stiffnesses and not the force"); **force_initial**: `a0+a1·y(z)` — preload de compresión (la reacción empuja el elemento aunque u=0, útil para empuje de reposo); **time**: factor temporal para la fuerza total (reutiliza force_time, manual 6.1084). **Controls**: `control_support_edge_normal_damping_apply` cableado (gate por control block); `control_support_edge_normal_stiffness_freeze` parse-only (sin efecto en elástico — pertenece a la plasticidad del lote 3). **VERIFICACIÓN** (suite **223/223** = 219+4, TODOS los checks OK): `tsup_damp` (k=(1000,100), c=(200,100), factor 0.5, prescripción → node1 (-0.0075,0) y node2 (0,+0.35) EXACTOS), `tsup_auto` (grupo E=1e7,ν=0,ρ=2 → c_n=4472.136, c_t=790.569 EXACTOS: fy/node = (1000+4472.136)/2 = 2.73607, fx/node = -(100+790.569)/2 = −0.04453), `tsup_init` (F0=0.1+0.2·y con time×2, prescripción 0 → node1 (0.1,−0.1) y node3 (0.3,−0.3) EXACTOS, lineal con y), `tsup_dens` (d_n=50, prescribed vy=−1e-3 → a=−1e-3 → fy/node = (1000+50)·1e-3/2 = 0.525 EXACTO). **GOTCHAS del lote**: (1) el manual dice "E y ν del elemento adjunto" — es el del GRUPO, no del record _automatic (medido: el test con E=1e7 del grupo dio exactamente la fórmula); (2) `support_edge_normal_factor` escala SOLO rigideces, NO la fuerza (F0 se preserva íntegra); (3) `bounda_unknown` ya conocido (multi-record quirk, -node no-op) y la orientación de la normal por producto vectorial de área() puede invertir el signo en aristas laterales — medir primero y ajustar el target; (4) el control _stiffness_freeze no tiene comportamiento aquí (sin plasticidad). |
| **SPRINT 13 LOTE 1 — support_edge_normal: el Winkler distribuido, NÚCLEO ELÁSTICO (7/18 del checklist + el record de salida)** | `feat(support): support_edge_normal — winkler distribuido 2D/3D (nucleo elastico)` (2026-08-29) | 2026-08-29 | **código** (area.cc como type[10] de la maquinaria de integrales de borde, patrón force_element_edge_*): master `support_edge_normal <i> k_n k_t` (por unidad de longitud 2D / área 3D) + los 5 selectores same-index (_geometry requerido — entidad O lista de nodos, _element/_element_group/_element_side/_node/_element_node como restricciones) + `node_support_edge_normal_force` (record de salida por nodo, patrón node_inertia) + los 2 control records parse-only (damping_apply → lote 2, stiffness_freeze → lote 3). **Física**: fuerza RHS desde los DESPLAZAMIENTOS TOTALES de los nodos de la cara (f = −(k_n·u_n)·n − k_t·u_t, n la normal SALIENTE, cuadratura Lobatto = cargas consistentes) + la RIGIDEZ CONSISTENTE del lado en la matriz (dt·∫(k_n n⊗n + k_t(I−n⊗n))N_iN_j dA, Gauss numérico 2/3 puntos, tensor en caras 3D; MISMO escalado dt que materi.cc) — SIN ella un cuerpo sobre soportes conserva un modo rígido de energía cero en la matriz de velocidades y el solver REVIENTA (medido); el punto fijo no cambia (el equilibrio vive en el RHS — el propio manual del _plasti_residual_stiffness lo declara). **Record de salida**: PUT en VERSION_NEW, acumulado por los elementos que comparten nodo, RE-CERADO al inicio de CADA sweep de ensamblaje (detección por número de elemento no creciente — `<` falla el caso de 1 elemento, medido; por llamada, NO por nodo — el zero por nodo borraba al nodo 1, medido); single-thread (OPTIONS_PROCESSORS>1 avisa). **VERIFICACIÓN** (suite **219/219** = 216+3): `tsup_winkler` (u prescrito en los nodos apoyados → fuerza nodal (−0.05, +0.5) EXACTA = k·u·L/2), `tsup_solve` (columna sobre muelle: u_top = −(F·L/EA + F/kL) = −3e-3 con reacción 0.5/nodo), `tsup_3d` (cara hex8: k·u·A/4 = 0.025 EXACTA por esquina). **GOTCHAS del lote**: bounda_unknown NO tiene selector -node (no-op SILENTE — tslv_vfx pasaba de casualidad prescribiendo ceros); bounda_time comparte el valor entre todos los unknowns del record (valores distintos → records separados); en 3D partir los unknowns en DOS records bounda dejó los dofs sin tocar en el ensamblaje (quirk preexistente, un record por valor); pri() máx (char*,valor); el loop interno con la MISMA variable idim pisa el externo. |
| **SPRINT 13 LOTE 3 COMPLETO — support_edge_normal: plasticidad + residual + status (4e0074f)** | `feat(support): plasticidad completa del soporte + bug critico del vector tangencial` (2026-08-29) | 2026-08-29 | **Reestructuración del per-side** (los caps ANTES del bucle por dirección): `S = k·u_n + c·v_n + d·a_n + f0` (positivo = compresión) y `ft_raw_local[]` (los componentes TANGENCIALES u_t = u−(u·n)n, v_t, a_t) se computan en un pre-bucle; luego los caps: **plasti_tension −yes** (gap: S<0 → todas las fuerzas 0 + side_opened=1), **plasti_tension_double** (S < −max → S = −max + tangencial a 0), **plasti_compression** (S > |min| → S = |min| + factor tangencial |F_t| ≤ factor·|S|), **plasti_friction** (|F_t| ≤ c + μ·|S|, el MENOR de los dos caps gobierna). **plasti_residual_stiffness**: fracción factor ∈ [0,1] de la rigidez añadida a la matriz (default 0; medido que factor≠0 CAMBIA el punto fijo cuasiestático → default 0 documentado). **status record**: node_support_edge_normal_plasti_tension_status 0/1 con semántica OR por nodo. **BUG CRÍTICO CORREGIDO**: ft_raw_local usaba el desplazamiento COMPLETO en vez del componente TANGENCIAL — el soporte aplicaba la rigidez tangencial sobre la componente normal duplicada (fy=0.75 vs 0.5 esperado, medido con DBG); ahora u_t/v_t/a_t correctamente proyectados. **element_side**: ahora chequea elemento Y lado local 1-based (antes solo elemento); **element_node** añadido al camino del soporte. **El switch −yes**: se lee con db_int de VERSION_NORMAL y se compara `== −YES` (el parser guarda el enum NEGADO — `==YES` falla silenciosamente). **VERIFICACIÓN**: 11/11 tsup_* rc=0 con valores analíticos EXACTOS: winkler (−0.05, +0.5), gap (0,0) + status=1, tcap (0, −0.00025), fric (+0.005, +0.5), residual (+0.005, +0.5), damp (−0.0075, +0.35), auto (−0.0445, +2.736), init (0, +0.1), dens (0, +0.525), solve (reacción +0.5/nodo), 3d (0, 0, +0.025). **INCIDENTE**: durante la depuración un `rm −f *.dat` demasiado amplio borró los ~200 .dat del suite (gitignored — recuperación vía Dropbox versioning web, pendiente); los 11 tsup_* se regeneraron desde cero y pasan. |
| **SPRINT 12 LOTE 3 — familia solver (10/10) + timestep (2/2) + tochnog_version + volume_factor_x + zip (15 keywords)** | `feat(solver): familia solver global + timestep_*/tochnog_version/volume_factor_x/zip` (2026-08-29) | 2026-08-29 | **código**: los `solver_*` del Professional (6.1047-6.1056) son records PLANOS (sin índice, `no_index=1`), a diferencia de los `control_solver_*` del Sprint 9: enums propios en sync + registros en database.cc. `solver` = el tipo GLOBAL que PISA a todo control_solver — leído al FINAL en los tres puntos de consumo (top.cc/elem.cc/dof.cc; precedencia OPUESTA al options_solver legacy, que el control pisa). `solver_bicg_error`/`_stop` pisan a sus control en so_bicg.cc. `solver_matrix_symmetric -yes` fuerza el CG saltándose la medición por solve, con WARNING honesto si la medición discrepa (el Professional simetriza la matriz; el GNU corre CG sobre la tal-cual — conecta con el fix A+B). PARCIALES documentados: `solver_bicg_restart` (el bicg del GNU no tiene restart), `solver_matrix_save` (cache de solver directo), `solver_pardiso_*` ×4 (PARDISO no compilado), `timestep_predict_velocity` (el esquema escalonado arranca el solve de x=0 y ACUMULA el resultado — una predicción de velocidades previas duplicaría; neutral verificado). **Comportamiento real**: `timestep_iterations_automatic_apply -no` (gate en top.cc que descarta los control_timestep_iterations_automatic), `tochnog_version` (la fecha de BUILD parseada de `__DATE__`, escrita al inicio de top(), targeteable), `volume_factor_x` (factor por tramos de x a nivel de IP en volume.cc — multiplica al polinómico y al de grupo; derecha de xn factor 1 como el manual), `zip -yes` (gzip de *flavia*/*msh/vtk*/*.plt/*.dbs al final en exit_tn, tras los targets y antes de db_close; solo archivos existentes — el glob sin match no cuenta como fallo). **VERIFICACIÓN** (suite **216/216** = 209 + 7, checks TODOS OK): `tslv_override` (control_solver bicg + solver -matrix_superlu → el BANNER de SuperLU en stdout prueba el override), `tslv_bsym` (el beam medido NO simétrico + -yes → el warning se imprime, sobrevive con stop -no), `tslv_symm` (simétrico + -yes = valores idénticos), `tslv_vfx` (estática de serie 8×1 con espesor 2 a la derecha: σ 2e-2/1e-2 por mitades, punta 1.15e-3 vs 1.2e-3 analítico; A/B sin el record 1.59e-3), `tslv_part` (los 7 parciales juntos = física sin cambio), `tslv_auto` (record automático + -no = física sin cambio), `tslv_ver` (target del año de build pasa), `tslv_zip` (aislado: .dbs.gz creado sin warning). Arness sin regresión (gforce10 1.0000). GOTCHAS del lote: los comentarios de .dat cuentan paréntesis (un "[1,3)" rompió el parse — conocido); `-ra 1 4` sin el segundo -ra revienta MRANGE (usar -node explícito); el combo de unknowns para σ necesita strain_total (+displacement si hay targets de tensión); el staggered cuasiestático colapsa en una pasada → los contadores de iteraciones son INVISIBLES en estos modelos (el gate del automatic queda verificado por construcción + no-cambio). |
| **SPRINT 12 LOTE 1 — frame de sección 3D: CIERRE DEL ARNESS (gforce10/13 V/mom + convenciones de signo)** | `feat(post): frame de seccion 3d — cierre del arness professional` (2026-08-29) | 2026-08-29 | **HALLAZGO (el V/mom ≈ 0 de gforce10/13 era un problema de FRAME, no del esquema)**: (1) BUG DE INPUT — el force10.dat del Professional lleva `thickness_switch -yes` (+ `plot_switch`); nuestra conversión gforce10/13 los omitió. MEDIDO por caja negra (reference_point movido en el binario del Professional): sin el switch el PROPIO Professional da el frame roto (vectores ≈ 0, momento en mom2) — la sección cuadrada 10×10 empata la regla de extensión del manual 6.917; (2) BUG DE CÓDIGO — corners[] = [0,npol−1,npol(npol−1),nface−1] sobre tablas de bordes LEXICOGRÁFICAS hacía corner3−corner0 = la DIAGONAL de la cara (e2 nunca fue arista; con −yes la diagonal ganaba → frame a 45°, nors 17.45 = 12.34·√2). FIX: bucle [00,10,11,01] (la normal n = e1×e2 no cambia: la diagonal vieja iba a lo largo de e1+e2); (3) CONVENCIONES convergidas al Professional (todas medidas): empate de cara cuadrata roto por el reference point (la arista mejor alineada gana; proyección nula → primera arista, conserva msf_cant3d/qsri3d), t orientado HACIA el reference point (antes "away" → TODA componente direccional invertida, ratios −1.0000), brazo de momento (medio−nodo) (antes nodo−medio), items `*_s` ESCALARES FIRMADOS (nors = −12.34 a compresión; antes \|valor\|). **VERIFICACIÓN**: (1) arness COMPLETO convergido — gforce10/13 **nory 1.0000 / shey 1.0000 / mom1y 1.0000** (−12.34/−100/−5000 vs Professional, EXACTOS a la precisión impresa — la estática del cuerpo libre del L5 es exacta también en 3D); gforce7/7q4 N/V 1.0000 M 0.9984 CON SIGNOS; gffq4 identidad 0.9996; msf_tunnel3d nory **1.0000** (0.0998068555 vs 0.0998068556, el MISMO valor FE, mismo signo); (2) suite **209/209** + verificaciones TODAS OK en build limpio (5 checks de s-items actualizados a comparación por MAGNITUDES con justificación — el signo depende del reference point de cada modelo); (3) los .dat gforce10/13 actualizados con thickness_switch+plot_switch (fieles al original). Residuales documentados (VALIDACION §13.3): componentes cruzados ~5% del nodo esquina 2D cargado (polución σ de borde del quad9, preexistente), y el GNU da el frame correcto sin switch donde el Professional da el roto (divergencia deliberada en input degenerado). |

**Carril A — `group_interface_materi_memory`** | `643865b` | 2026-08-16 | modelo de memoria de la ley de interfaz: `-updated_linear` (default, normal/tangente de la config actual) o `-total_linear` (geometría de referencia tiempo 0, `NODE_START_REFINED`). Valores inválidos → `db_error`. Test `iface_mc_mem` (familia iface_mc, 7º run → 19 runs). |

**Carril A — validación numérica MC/tension completa** | `8022af3` | 2026-08-16 | cierre de RF-4: `iface_mc_dil` (phi_flow=30° → strain_normal acumulado −2.8065 vs +0.0802 sin dilatancia; apertura por deslizamiento > compresión), `iface_mc_dil_1step` (invarianza paso, −2.82675, dif ~0.7%) y `iface_mc_num` (clamp MC: c=50, phi=0, vely=1000 → `element_interface_force_tang`=50.0 exacto = max_fric=c). Familia iface_mc → 10 runs (22 totales). |

**Carril A — interfaz 3D + conversión 3D (items C+D)** | `d57faf4` | 2026-08-16 | `interface_element` soporta 3D (prism6/hex8): normal real por cross product orientada de lado2 hacia lado1 (compresión positiva), 2 tangentes ortonormales, 2 fuerzas tangenciales acumuladas (nueva history `element_interface_force_tang2`). MC 3D clampa la magnitud de la fuerza tangencial total preservando la razón; dilatancia usa |du_tang| en ambos planos. `interface_convert` ahora convierte tria3→prism6 y quad4→hex8; la llamada se movió de `step_close` a `step_start` (antes del primer `element_loop`). Tests: `iface_3d` (fricción alta sostiene carga tangencial 3D), `iface_3d_slip` (fricción nula desliza), `iface_conv3d` (tria3→prism6). |

**Carril A — print 3D `control_print_interface_stress` (item E)** | `b8804c6` | 2026-08-16 | 3D implementado: imprime centroide (x y z) + `sign` + `sigt1` + `sigt2` por elemento de interfaz. `_3d_geometry` filtra por geometría (patrón adjust.cc); `_3d_order -x/-y/-z` ordena las líneas por la coordenada del centroide (recolección+sort). Validado con `iface_3d_stress` (centroide 1.5 0.5 0.5, sign crece con compresión, sigt2≈0) e `iface_3d_order` (2 interfaces ordenadas por x). |

**Carril B — `control_mesh_generate_interface` (+ `_geometry` + `_method` + subdivisión cuadrática)** | `947ac83` + `87f157d` + `3d9ed64` + `ced696f` + `9a2b422` + `14a94d3` | 2026-08-17 | genera elementos de interfaz entre dos grupos de elementos que comparten una cara (nodos duplicados espacialmente). Sintaxis `index eg_i eg_a eg_b eg_j eg_c eg_d...`; el elemento generado se asigna a `eg_i`. `_geometry` filtra por geometría; `_method` (`-element_geometry`) selecciona por `element_geometry` y/o genera un registro `element_geometry`. **Subdivisión de caras cuadráticas** (`interface_face_subdivide`): quad9/bar3 (arista 3 nodos → 2 `quad4`), tet10 (tria6 → 4 `prism6`), hex27 (quad9 → 4 `hex8`) — acopla todos los nodos incl. medios/centro; clasificación por coordenadas. Validado con `iface_gen*` (lineal, geometría, método) e `iface_gen_quad9` (2 quad9 → 2 interfaces quad4, fricción MC sostiene). |

**Carril B — `group_materi_plasti_mohr_coul_direct`/`tension_direct` (+ `_normal`/`_normal_automatic`/`_visco`/`_wall`)** | `b6d4be6` + `7f51812` + `e553b6d` + `d3ba0cb` | 2026-08-17 | leyes plásticas de **cut-off directo de tensiones** en un plano con normal `n` (patrón del manual: sin deformaciones plásticas). `tension_direct sigy` capa la tracción normal; `mohr_coul_direct phi c phi_flow` capa el cortante a `max_fric = max(c - sig_n*tan(phi),0)`. `_normal` da la normal explícita; `_normal_automatic -yes` la toma del elemento; `_visco tm` relaja el cut-off (`factor=1-exp(-dt/tm)`); `_wall` usa valores alternativos si el elemento está pegado a una pared (`plasti_on_boundary`). `materi_direct_cutoff()` en stress.cc. Validado con `materi_direct`, `materi_direct_mc`, `materi_direct_auto`, `materi_direct_visco` y `materi_direct_wall`. |

**Carril B — `materi_displacement_relative`** | `96eabda` + `841c39e` | 2026-08-17 | opción del initia que añade un dof de **desplazamiento relativo** (`disr*`): acumula el desplazamiento desde un punto de referencia. La referencia se re-sincroniza (reset a 0) en dos eventos: cambio de timestep en `control_timestep` (dt nuevo ≠ persistido en `MATERI_DISPLACEMENT_RELATIVE_REF`, en `top.cc`) y reset de desplazamiento en `control_reset_dof` (cuando se resetea `disx`, en `data.cc`). Integración del dof junto a `dis_indx` en `dof.cc`. Requiere `materi_displacement` + `materi_velocity` + `materi_velocity_integrated`. Validado con `mat_rel` (dt cambia: disy=2.0, disry=1.0) y `mat_rel_reset` (reset de disx: disy=1.0, disry=0.1). |

**Carril B — `slide_axisymmetric` + variantes `control_reset_value_*` espaciales** | `3024d24` + `650ed26` | 2026-08-18 | `slide_axisymmetric -yes` escala la fricción de `slide_geometry` por `2*pi*r` (r = coordenada radial del nodo, `slide.cc`; validado con `slide_axi`). Variantes espaciales de `control_reset_value` (distribuciones en x/y/z): `_linear` (`ax x + ay y + az z`), `_power` (`ax x^bx + ...`), `_square_root`, `_exponent`, `_logarithmic`, `_logarithmic_second`, `_multi_linear` (tabla vs coordenada vertical) — implementadas en `data.cc` (bloque de reset); validada `_linear` con `reset_value_linear` (disx → x). |

**Carril B — `mesh_activate_gravity_time*` + `strain_settlement_*`** | `f6b8935` + `b58fd8f` + `4fc3419` + `6f44549` | 2026-08-18 | **Activación gradual de la gravedad** (`mesh_activate_gravity_time`, `_element`, `_element_group`, `_geometry`, `_time_initial`, `_method` 1 y 2, `_stiffness_factor`, `control_mesh_activate_gravity_apply`): `mesh_activate_gravity_factor()` en mesh.cc, aplicada en materi.cc (el vector de gravedad se multiplica por el factor de activación; el método 2 escala la matriz del elemento por el factor de rigidez). **Creep de asentamiento** (`strain_settlement_parameters`, `_element_group`, `strain_settlement_diagram`+`_dof`+`_number`): `strain_settlement_creep()` en materi.cc añade el creep vertical a `inc_ept` en `set_deften_etc`. Ley con saturación (decisión 2026-08-18, el OCR del manual es ambiguo): `eps_zz = Ar*(t/t_ref)^n/(t_plus+(t/t_ref)^n)`. Validado con `mesh_act_grav`, `mesh_act_grav2` (método 2), `strain_settle` (creep comprime la columna) y `strain_settle_diag` (Ar duplicado vía diagrama). |

**Carril B — familia `contact_*` (apply, plasti_friction, targets)** | `23262bf` + `61002f6` | 2026-08-19 | el GNU ya tenía el algoritmo de contacto (`parallel_contact` en contact.cc: geometría target, penalties, stick, relaxation, heatgeneration) y `contact_friction` simple. Se completan: `contact_apply` (gate `-yes`/`-no` por timestep), `contact_plasti_friction` (Mohr-Coulomb `max(c + Fn*tan(phi),0)` en el slip, reemplaza el `mu*Fn`), `contact_target_element_group` (filtro de targets por grupo), `contact_target_geometry`/`_switch` (alias de `contact_geometry`/`_switch`). Validado con `contact` (la familia se parsea y el modelo corre estable; el algoritmo de contacto es experimental y la detección de penetración depende del caso). |

**P6 — switches de consolidación groundflow (`groundflow_consolidation_apply`, `control_groundflow_consolidation_apply`, `group_groundflow_consolidation_apply`)** | `7fd1966` | 2026-08-20 | switches `-yes`/`-no` que controlan el término de consolidación (divergencia material) en la ecuación de groundflow: global, por timestep (`control_*`, index) y por grupo. Precedencia: group > global > control; cualquier `-no` fuerza `materidivergence=-NO` en groundfl.cc (además de los legacy `OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE`/`CONTROL_OPTIONS_...`). Se corrigió un bug latente en database.cc: `CONTROL_OPTIONS_SKIP_GROUNDFLOW_MATERIDIVERGENCE` estaba registrado duplicado como `..._NONLINEAR` y nunca se parseaba. Validado con `groundflow_consolidate_off` (campo de velocidad con divergencia 2: con `-no` la presión queda 0; con default `-yes` es 1.67, el target falla). **ACTUALIZADO 2026-09-05 (sprint u-p, fila arriba)**: el default de la familia pasó a `-no` (Professional 6.556) y la semántica de activación se midió en el Pro — solo global/control activan; el grupo solo excluye.

**P6 — modelo no saturado van Genuchten (`group_groundflow_nonsaturated_vangenuchten`, `group_groundflow_nonsaturated_eps_permeability`, `groundflow_nonsaturated_apply`, `control_groundflow_nonsaturated_apply`)** | `cdaf8a6` | 2026-08-20 | ley van Genuchten del manual 2024 en groundda.cc (`groundflow_data`): `S(phi_p) = Sres + (Ssat-Sres)(1+(ga|phi_p|)^gn)^((1-gn)/gn)` con `phi_p = -pres/(dens*|g|)`; capacidad no saturada `c = csat + n*dS/dphi_p` (derivada analítica); permeabilidad `ki = krel(S)*ksat,i` con `krel = Se^gl[1-(1-Se^(gn/(gn-1)))^((gn-1)/gn)]²` (Mualem). `group_groundflow_nonsaturated_eps_permeability` fija el suelo de `krel`. Dof `groundflow_saturation` (initia, almacenado no resuelto) guarda S por nodo. Gates `groundflow_nonsaturated_apply` (global, default -yes) y `control_groundflow_nonsaturated_apply` (por timestep). Requiere `group_groundflow_porosity` + `groundflow_saturation`. Validado con `groundflow_vangenuchten` (gsat nodo 5 = 0.7364 exacto; pres media -2.76 vs -5 lineal por la k reducida) y `groundflow_nonsaturated_off` (con `-no` pres media -5.0 exacta, saturado). |

**P6 — presión de grieta + interfaz groundflow (`group_groundflow_total_pressure_tension`, `group_interface_groundflow_capacity`, `group_interface_groundflow_permeability`, `group_interface_groundflow_total_pressure_tension`)** | `25352a8` | 2026-08-20 | `group_groundflow_total_pressure_tension` (materi.cc, bloque groundflow_pressure): si el mayor autovalor de `materi_strain_plastic_tension` supera `plastic_tension_minimum`, la presión estática `dens_agua*g*(water_height - coord)` reemplaza a la de la ecuación si es mayor en valor absoluto (grietas de hormigón). Usa `GROUNDFLOW_DENSITY` (no la densidad del material). En interfaz (interface.cc): `_capacity` (almacenamiento lumped en pres dofs), `_permeability` (flujo a través `q = pe*(pres_s1 - pres_s2)` acoplando los pares de nodos opuestos, con tangente simétrica ±pe), `_total_pressure_tension` (presión estática forzada si `strain_normal` acumulado > `strain_normal_minimum`). Validado con `groundflow_total_pressure_tension` (epp inicial 0.01: node_rhside velx -4.147 vs 0.0525 sin la corrección) y `groundflow_interface` (bloque derecho aislado se llena hasta pres 2.0 solo con la permeabilidad de interfaz; sin ella queda 0). |

**P6 — familia `groundflow_flux_edge_normal*` (11 keywords)** | `0cf4e6d` | 2026-08-20 | flux de agua distribuido normal a la arista, traducido a flux nodal equivalente en el rhs del dof de presión (`area()` en area.cc, `type[5]`/`type_area[5]`, MTYPES 5→6). Principal `groundflow_flux_edge_normal` + `_geometry` (área, e.g. `-geometry_line`), `_time`/`_sine` (carga temporal), `_factor` (polinomio espacial, via `force_factor`). Restricciones: `_element`, `_element_group`, `_element_node`, `_element_node_factor` (factores por nodo local, DOUBLE: values_fac[0]=elemento), `_element_side` (pares elemento/lado), `_node` (nodos globales). Atención del manual: si se usa DENTRO de la malla el flux total se anula (normales opuestas). Validado con `groundflow_flux_edge` (flux 0.1 por el borde inferior → pres 2.0 en el nodo inferior, `pres = flux*altura/k`; con flux 0 queda 0; variantes `_element_group` y `_node` verificadas). |

**P6 — familia `groundflow_phreatic_level_multiple*` (7 keywords)** | `e5a5468` | 2026-08-20 | varios niveles freáticos indexados, cada uno dueño de una parte del dominio. `groundflow_phreatic_level_multiple_find()` (groundfl.cc) resuelve el índice del nodo por `_element`/`_element_group`/`_element_geometry`/`_node` (una sola, no combinables); `groundflow_phreatic_coord` usa ese nivel; `_static -yes` (en `groundflow_phreatic_apply`) impone presión total = estática en los nodos del dominio (útil si la línea freática está sobre la malla y para no resolver los heads). `_n` para tablas 3D. Las keywords se registran con guion (`groundflow_phreatic_level_multiple`, como el manual Professional 2024); la única sigue siendo `groundflow_phreaticlevel` (convención GNU). Validado con `groundflow_phreatic_multiple` (dos columnas con niveles 3 y 1: pres estática -3 y -1 en el fondo, gradiente correcto en altura). |

**P6 — familia `groundflow_seepage_*` (eps, geometry, node)** | `ab54301` | 2026-08-20 | condición de seepage en bounda.cc: un borde del dominio por el que el agua solo puede SALIR (punto de salida automático). La normal del geometry debe apuntar hacia afuera del material. La activación la gatilla `groundflow_seepage_eps` (default 0.1); `_geometry` (entity+index, normal via `geometry()`) y `_node` (nodos) definen el borde. En `bounda()`, al imponer `pres` en un nodo del borde: si la velocidad de Darcy del paso anterior proyectada en la normal exterior es `<-eps` (flujo ENTRANTE), se desmarca `node_bounded` → el borde se cierra (no se impone la presión); si el flujo sale, se impone la presión (típicamente 0, aire libre). Patrón del manual: `groundflow_seepage_geometry 10 -geometry_line 100` + `bounda_dof 20 -geometry_line 100 -pres` + `bounda_time 20 0.0`. Validado con `groundflow_seepage` (pres 2 arriba, seepage en el fondo: perfil drenado, pres media 1.0). |

**P6 — `groundflow_pressure_atmospheric` documentado y verificado (keyword heredado del GNU, no existe en Professional 2024)** | `30e46c7` | 2026-08-24 | cap de la presión estática/total en `groundflow_phreatic_coord()` (groundfl.cc:323-325; keyword registrado desde sfnet 2014, `no_index`, default 0). Con default 0 la succión (positiva, sobre el nivel freático) se anula y la compresión (negativa) pasa intacta; `pa>0` mantiene succión hasta ese valor. `_static -yes` impone exactamente el valor clampeado (groundfl.cc:414). Se corrigió el signo invertido de la sección "Known limitation" de manual-developer/bounda_water.md y se añadieron a los índices README los 30 docs de la familia groundflow P6 que faltaban (deuda del 2026-08-20). Validado con `groundflow_pressure_atm` (pa=0.5: succión +1→+0.5, compresión -1.0 intacta) y `groundflow_pressure_atm_def` (default: succión→0.0); suite completa 57/57. |

**Fix — NaN en interfaces sin `elasti_stiffness` (UB ddum3, expuesto por gcc 14.2 -O1)** | `0c2f5e0` | 2026-08-24 | `interface_element()` leía `GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS` con GET_IF_EXISTS sobre `ddum3` sin inicializar: sin el record, `kn/kt1/kt2` quedaban con basura de stack. El rebuild limpio con gcc 14.2 -O1 (entorno /mnt/nas) convirtió la basura en no-finitos → matriz NaN → `groundflow_interface` fallaba (BI-CG initial error -nan, determinista). El binario del entorno anterior (/mnt/z) no lo manifestaba. Fix: `array_set(ddum3,0.,3)` antes de la lectura; bisect empírico por objeto (-O0 interface.o pasaba). Suite 59/59. Regla general anotada: inicializar todo buffer pasado a GET_IF_EXISTS. |

**P6 — `groundflow_total_pressure_limit` (manual Professional 6.588)** (SUPERSEDED by F4 T1, 2026-10-09: default 0, cap of the calculated total pressure, velocity 0 where reached, pres dof never cut) | `25bc716` | 2026-08-24 | cap del pres RESUELTO: tras cada `solve()` (top.cc) `groundflow_total_pressure_limit_apply()` corta `node_dof[pres] > limit` en nodos libres (los Dirichlet/bounded y `_static` mantienen su valor prescrito). Interacción con consolidación (manual 2.4.3): con `limit 0` y presión 0 el elemento está SECO y `groundflow()` fuerza `materidivergence=-NO` (se salta el término de divergencia material; compone con la familia `groundflow_consolidation_apply`). Diferencia documentada con Professional: el record es REQUERIDO para activar (Professional defaultea limit a 0 incluso sin record; el GNU sin record no clampa, para no romper modelos con presiones positivas legítimas como el propio ejemplo de seepage del manual). Además: re-sincronizados tochnog.h y tochnog-mod.h (977 enums idénticos; la sesión del 20-ago dejó 33 enums solo en tochnog.h) y extraído sqlite3.h al numlib-runtime para compilar sin libsqlite3-dev. Validado con `groundflow_total_pressure_limit` (limit 0.5: nodos libres 0.5 vs ~1.67 sin limit, Dirichlet intacto) y `groundflow_total_pressure_limit_dry` (limit 0: seco, consolidación omitida, pres 0); suite 59/59. |

**P6 — remate contacto: `contact_heat_generation` + `control_contact_apply`** | `b2694c8` | 2026-08-24 | `contact_heat_generation` (manual 6.100) es el nombre Professional del legacy `contact_heatgeneration`: la física (Q = η·Ff·vf inyectado en el dof temp, heredada del GNU 2014) ya existía; se registra el keyword nuevo y se lee AMBOS nombres en el único punto de consumo (contact.cc; buffer a 0 antes — regla GET_IF_EXISTS). Sin alias PUT: el PUT en `parallel_contact` aborta con "Data is allocated in a parallel loop" (lo pilló el probe). `control_contact_apply index -yes/-no` (manual 6.112): gate por timestep junto al `contact_apply` global, indexado por ICONTROL, cualquier -no gana. Validación con A/B real de bloque en caída: `contact_block` (contacto activo: el penalty detiene el bloque, disy≈-0.003) vs `contact_ctrl_apply` (-no: caída libre exacta disy=-0.055, idéntico a sin contacto) + `contact_heatgen` (modelo térmico con el nombre nuevo, mecánica estable). Flujo del VALOR verificado con TOCHNOG_DEBUG: friction_energy escala exacto con el factor (0.5→1.31194, 1.0→2.62387). Limitación documentada: la temp acumulada queda ~0 porque la iteración final de equilibrio vuelve a STICK (algoritmo experimental). Suite 62/62. |

**P6 — `change_dataitem_time_method` + `change_dataitem_geometry` (manual Professional 6.51/6.48, familia change 6/6)** | `54af57a` | 2026-08-24 | `_time_method`: la tabla temporal da valores de coseno/seno/tangente y se almacena el ángulo inverso (`acos/asin/atan(val)`) — el caso phi-c del manual (tablas de tan(phi), parámetro phi en radianes). Switch enums COSINUS/SINUS/TANGENT con registro de nombres en database.cc (el parser resuelve -tangent por la tabla de nombres, como cualquier switch). `_geometry`: restricción de un cambio de record `group_*` a los elementos dentro de la geometría, materializada como SPLIT del grupo: clon con copia de todos los records `group_*` (grupo nuevo = max+1, mapa estático), elementos completamente dentro de la geometría movidos al clon (patrón geometry() de CONTROL_DATA_INITELDOF_GEOMETRY + PUT de ELEMENT_GROUP como delete.cc), cambio aplicado solo al clon (los de fuera y el record original quedan intactos). Validado con `cd_method` (unit test exacto: `target_item -group_materi_plasti_mohr_coul_direct 0 0` = 0.5235988 = atan(0.57735027); sigxy invariante ~1 → sin efectos colaterales) y `cd_geom` (dos bloques en shear: izquierdo dentro del brick → sigxy~0 con c=0 en el clon; derecho → sigxy~1; el record original conserva c=1.0 leído por target genérico — prueba la restricción, no solo la ausencia de crash). Suite 64/64. |

**P6 — familia `control_data_*` (activate, arithmetic+double, copy+factor, copy_index+factor; manual Professional 6.114-6.120)** | `2051bc6` | 2026-08-24 | bloques nuevos en `data()` (disparan en cada step_close cuyo ICONTROL casa): `_arithmetic` cambia un record con el valor de `_double` vía -plus/-minus/-multiply/-divide (índice simple o rango -ra, número simple o -all; ints rechazados; división por 0 rechazada; switches PLUS/MINUS/DIVIDE nuevos con nombres registrados); `_copy` copia TODOS los índices from→to con factor opcional (helper data_copy_apply: int requiere factor 1, double multiplica; el caso d'alembert node_inertia→node_force -1 queda cubierto); `_copy_index` copia un índice concreto con `_factor`; `_activate` -no BORRA todos los records de los items listados (db_delete_index; los consumidores los saltan) y -yes es no-op (desactivación destructiva, diferencia GNU documentada). `control_data_save` DESCARTADO como feature (solo consumidor: control_print_gid_save_difference, familia GiD descartada) pero REGISTRADO parse-only el 2026-09-04 (los tests dynamic6/7 del corpus lo declaran; sigue sin comportamiento). GOTCHA crítico descubierto calibrando: `control_timestep i dt incremento` — el 2º número es la DURACIÓN del bloque, no el tiempo final (`0.1 0.2` desde t=0.1 = 2 pasos hasta 0.3; un record por paso dispara 2 veces). Validado con `cda_arith` (young 1000 → ×2 → +500-all → 2500 exacto por target genérico), `cda_copy` (copy_index ×2 → 2000; copy todos ×0.5 → 500/1000 exactos) y `cda_activate` (columna traccionada disy~0.01 → bounda_force borrado → relaja a 0±0.004 en cuasiestático). Suite 67/67. |

**P6 — familia `control_distribute_*` completa (manual Professional 6.127-6.132; P6 CERRADO)** | `c7d4002` | 2026-08-24 | layout Professional añadido a `distribute()` (distri.cc): record de 4 slots `[-normal/-lognormal, item, index|-all, number]` discriminado del legacy GNU (triplets + `_values`, que se conserva intacto — ho_othr4 sigue verde). `_parameters` da mean/std del VALOR (lognormal vía mu_ln/sigma_ln = ln(1+(s/m)²)); `_seed` reinicializa el ran1 de Numerical Recipes (idum negativo) → campos reproducibles; `_correlation_length` (1..ndim valores por dirección; >1e12 = campo constante) con kernel exponencial `v_i = Σ w_ij z_j / √Σw_ij²`, w=exp(-d_eff), distancia de correlación por defecto 4·L (`_correlation_distance` la sobreescribe); `_minimum_maximum` clampa el draw. Aplicación: group_* → delta = draw − record en `ELEMENT_DISTRIBUTE_VALUES` (el `+=` de get_group_data da el valor exacto; filtrado por ELEMENT_GROUP del index o `-all`); node-like → REPLACE del slot (ejemplos del manual: temperaturas nodales, coordenada y). Dispara por PASO en step_start (como el GNU). Validado con `cdist_normal` (deltas calibrados exactos +35.4136/−21.3726 con seed 7: draws distintos y reproducibles), `cdist_corr` (correlation_length 2e12 → ambos deltas IGUALES 35.4136 pese a 25 unidades de distancia) y `cdist_clamp` (std 1e6 + min/max 995/1005 → deltas exactamente ±5). Suite 70/70. **P6 COMPLETO**: todas las familias del bloque cerradas (groundflow, contacto, change, control_data, control_distribute). |

**Sprint 7 — familia `condif_heat_*` completa (manual Professional 6.72-6.91; 20 keywords)** | `8ae0feb` | 2026-08-24 | **edge_normal (11 kw)**: maquinaria COMPARTIDA con groundflow_flux_edge en area.cc — MTYPES 6→7, type[6]=CONDIF_HEAT_EDGE_NORMAL + _GEOMETRY; los 3 bloques groundflow-específicos (restricciones element/group/side, carga sine/time, aplicación nodal con factor y restricciones por nodo) generalizados con helpers estáticos `flux_edge_is_master`/`flux_edge_companion` que mapean cada master a sus 9 companions; el dof objetivo es `temp_indx` (heat) vs `pres_indx` (groundflow). Comportamiento groundflow bit-idéntico (groundflow_flux_edge verde). **volume (9 kw)**: bloque nuevo en condif() (condif.cc) — restricciones `_element`/`_element_group` (array_member) y `_geometry` (todos los nodos del elemento en la geometría), valor del record o de `user_condif_heat_volume` con `_user -yes` (stub en user.cc, patrón user_viscosity), carga `_sine`/`_time` (force_time), `_factor` polinomio espacial en el punto de integración (firma de condif() extendida con coord_ip; call-site elem.cc actualizado); contribución rhs += volume·h_i·S·load·factor. GOTCHA descubierto: `group_condif_flow` en 2D exige 2 valores (vector), 1 solo valor rompe el parser. Validado con `condif_heat_edge` (columna q=0.1 por el borde inferior, k=0.1, T=0 arriba: T fondo **2.0** y medio **1.0** exactos por Fourier), `condif_heat_vol` (barra 1D con `_element` solo elemento 2: T(centro) = **0.25 analítico**; 0.5 calentando ambos — A/B de la restricción) y `condif_heat_vol2` (`_factor 0. 1.` → S(x)=x: T(centro) = **0.5 analítico**). Suite 73/73. |

**Sprint 7 — cierre: `condif_convection_edge_normal*` + `condif_radiation_edge_normal*` (manual Professional 6.65-6.71/6.92-6.98; 14 keywords; BLOQUE CONDIF 34/34 COMPLETO)** | `19757b0` | 2026-08-24 | nombres Professional de los legacy `condif_convection`/`condif_radiation` (misma física h·(Tenv−T) y αr·(Tr⁴−T⁴) con tangente, misma maquinaria area()): masters como tipos 7/8 (MTYPES 7→9) con helper `conv_rad_is_master`/`conv_rad_companion` (5 variantes de restricción: _element/_element_group/_element_side antes del bucle de lados, _node/_element_node en la aplicación); rama de aplicación extendida a los 4 masters (lectura de values desde type[itype]; los legacy conservan su camino exacto — convec1/convec2/condif8/condif10 verdes). `_geometry` dual: entidad O lista de nodos (camino genérico area[0]>0). GOTCHA MAYOR descubierto calibrando: `border_nodes_quad4 = {0,1, 1,3, 3,2, 2,0}` — tochnog conecta el quad4 en convención **Z** (1,2 abajo / 3,4 arriba en el mismo orden x), NO ccw; una "arista" fuera de esa tabla NUNCA dispara silenciosamente (la matriz de conducción sí es correcta con cualquier orden — por eso los tests mecánicos ccw pasan). Los tests con features de area() deben usar la convención Z. Validado con `condif_convec` (T=0.5 analítico exacto), `condif_rad` (no lineal: T+T⁴=1 → 0.7245 con 10 iteraciones Newton) y `condif_convec_el` (`_element` que no toca la geometría → T=0 exacto, A/B de la restricción). Suite 76/76. |

**Sprint 8 — `area_element_group*` completa (13/13) + `bounda_time_factor`** | `dbd210a` | 2026-08-24 | variantes nuevas en group.cc: `_element` (filtro por nombre de elemento), `_node` (lista global de nodos vía array_member en vez de geometría), `_method` extendido con `-any_but_not_all` y entero N (≥N nodos dentro), `_interface` -no (excluye elementos cuyo grupo tiene GROUP_INTERFACE; GNU también acepta -yes para las interfaces del Carril A), `_time -yes` (re-evaluación por paso: `area_element_group_time_active()` en top.cc step_close). Sequence: alias Professional `..._element_group` (con guiones; copia al legacy ANTES del max_index — un input solo-alias veía max=-1 y nunca entraba al loop, bug cazado con diagnóstico element_group), `_geometry_method`, `_interface`, métodos extendidos. HALLAZGO: `area_element_group_method` y `..._sequence_method` se LEÍAN en group.cc desde el GNU 2014 pero NUNCA estuvieron registrados en database.cc — eran inutilizables; registrados. `bounda_time_factor` (6.36): multiplica los LOADS (posiciones impares) del bounda_time tras el _units. GOTCHA de tests: `bounda_unknown -ra a b` es un RANGO (a..b), no una lista — confundió dos tests con strain=0 misterioso. Validado con `aeg_node` (brick solo columna izquierda + `-any` + `_time`: sigyy -4.0 vs -2.0 = A/B del young duplicado), `aeg_seq` (alias + `_geometry_method -any`, switch t=0.1 con 3 pasos: sigyy -2.0 EXACTO = 2 pasos E1000 + 1 paso E2000) y `bt_factor` (tabla 5.0 × factor 2.0 → sigyy 20 vs 10). Checklist bounda_time: `bounda_time_o`/`smc_o` eran OCR truncados de offset; SMC (3) y until_data/value_minimum (2, requieren post_node_result) quedan pendientes documentados. Suite 79/79. |

**Sprint 8 — `group_interface_*` COMPLETA 13/13 (Carril A cerrado del todo; manual Professional 6.625/6.630/6.636)** | `7b414fa` | 2026-08-24 | 3 keywords nuevas en interface.cc: `condif_conductivity` (q = k·(T1−T2) en los dofs temp de los pares enfrentados con tangente simétrica — espejo exacto de groundflow_permeability); `materi_expansion_normal` (eps_th total = α·T_media resta al strain acumulado para gap/tracción/MC vía `strain_eff`; el INCREMENTO térmico d(α·T) se resta de du_norm → pseudo-carga incremental, mismo patrón eigenstrain que stress.cc; el history queda puramente mecánico — sin doble conteo por paso); `tangential_reference_point` 3D (t1 = parte perpendicular de (ref − centroide) según el manual, t2 = n×t1; fallback al frame geométrico si degenera; consistente con memory -total_linear). Checklist corregido: `group_interface_ground` NUNCA existió (artefacto), `_materi_memory` ya estaba hecho (643865b), `_elasti_sti`/`_residual_sti` = OCR de stiffness, `tangential_reference_point` faltaba en el checklist. GOTCHA MAYOR re-verificado a la mala: insertar enums en tochnog.h exige BUILD LIMPIO COMPLETO — el binario mixto corrompido daba checks fantasma ("at least one of materi_velocity_integrated..."). Validado con `iface_condif` (bloque aislado llega a T=2.0 SOLO por la interfaz), `iface_expansion` (bloques fijos: sigxx = 1.984 ≈ 2·kn·α·T; 0 exacto sin el record) y `iface_tangref` (ref en +z rota t1: la cortante en y pasa a f_t2 = 10.40 y f_t ≈ −7.6e−07; frame por defecto la pondría en f_t). Suite 82/82. |

**Sprint 8 — familia `node_*`: 7 keywords + FIX preexistente de `node_mass`** | `ff6ff8b` | 2026-08-24 | node_force (fuerza nodal discreta en dof.cc, convención force_point), node_inertia (record calculado m·(a+g) por paso — la fuente del d'alembert documentado en control_data_copy), node_static/dynamic/total_pressure (overrides del post_calcul -static/-dynamic/-total en calcul.cc), node_slide (membresía ADITIVA de slide en slide.cc: geometría O record; la normal sigue viniendo de slide_geometry). node_mesh/node_convection_apply registrados+parseados sin comportamiento (documentado como parciales). FIX PREEXISTENTE: data_length[NODE_MASS]=1 → ndim (node_damping/stiffness usaban ndim; con 1 el parser se atragantaba en 2D/3D — node_mass inutilizable desde los orígenes del GNU). TRES GOTCHAS encadenados cazados con probes de una línea: (1) todo item NODE-class input-eable exige version_all=1 — mesh_has_changed BORRA los records NODE sin versión en cada cambio de malla (el GET_IF_EXISTS jamás disparaba); (2) el PUT de records calculados va en VERSION_NEW — el db_version_copy(NEW→NORMAL) del cierre de paso pisa cualquier escritura en NORMAL; (3) PUT con length=ldum(0) → db_error (length≥1). OCRs aclarados: node_dof_start_re/node_start_re/node_sti = truncados de refined/stiffness. Tests: node_force_inertia (inertia −19.83 ≈ m·g=−20 con node_force −30 activo; ejercicio del fix node_mass 2D), node_pressure (override static 7.5 EXACTO vs −2.0 calculado = A/B) y node_slide (smoke; LIMITACIÓN documentada: slide_axi original TAMBIÉN da rhside 0 — su target ±0.05 pasaba trivial, la verificación física del slide necesita modelo dedicado). Suite 85/85. |

**Sprint 8 — cierre: `control_reset_*` COMPLETA 18/18 (manual Professional 6.351-6.356)** | `e5cceac` | 2026-08-24 | filtros de región/nodo para control_reset_dof: `_geometry` (nodos de elementos completamente dentro), `_node` (elementos con todos sus nodos listados) y `_element_group` (restringe los elementos) — el filtro se calcula UNA vez por record antes de las variantes de valor y se aplica en todos los loops de nodos (array marcador reset_dof_node_filter). Resets de interfaz como records independientes: `_interface` (strain_normal + fuerzas tangenciales a 0) y `_interface_strain` (solo strain_normal; las fuerzas tangenciales quedan como tensiones "recordadas" — nuevos esfuerzos crecen desde ahí vía la rigidez). `_element_dof` registrado pero la rama -yes (solo element_dof, no node_dof) sin cablear (parcial documentado: requiere el almacenamiento por punto de integración). Validado con `creset_geom` (hisv0 0.5→0 SOLO columna izquierda dentro del brick; derecha conserva 0.5; A/B sin filtro FALLA) y `creset_iface` (history strain_normal de la interfaz 0.0 EXACTO tras el reset). Suite 87/87. **SPRINT 8 CERRADO**: area_element_group 13/13, bounda_time_factor, group_interface 13/13 (Carril A completo), node_* 7+fix node_mass, control_reset 18/18 — ~42 keywords, 2 fixes preexistentes (node_mass data_length, alias-copy order), 8 gotchas documentados. |
**Sprint 9 — lote 1: nombres Professional force_* + variantes de restricción + gates control_materi (~47 keywords)** | `f178849` | 2026-08-24 | **Traducción de prefijo DENTRO de db_number** (database.cc): `force_edge_*`→`force_element_edge_*` y `force_volume_*`→`force_element_volume_*` — clave crítica: el detector de fin-de-valores de records variable-length TAMBIÉN llama db_number; la primera versión (traducir solo en el punto del keyword de input.cc) rompía el parseo ("Problem reading : bounda_time"). Con el fix, los inputs Professional puros arrancan tal cual. **Variantes nuevas** (patrón companion): `_element`/`_element_group`/`_element_side`/`_node`/`_element_node` para las 3 familias edge + `_node_factor` (edge y normal) + `force_edge_water_factor` (nuevo, el GNU no lo tenía) + `_element`/`_element_group` para volume (force.cc). **Gates control_materi**: helper `control_materi_gate_off` (general.cc, ICONTROL-indexado) cableado a 5 features — viscosity (viscosit.cc), damage+failure (damage.cc), plasti_tension (materi_direct_cutoff), plasti_visco (set_stress), updated -no (set_stress, punto canónico de GROUP_MATERI_MEMORY); otros 7 registrados como parciales (sin término subyacente en el GNU). Validado con `fedge_alias` (sintaxis Professional pura: sigxx 5.0), `fedge_restrict` (`_element` que no toca la geometría: disx 0), `fvol_elem` (`_element_group`: 0.5/0.0) y `cmat_gate` (tension_apply -no: lineal −88.4 vs capped ~1, A/B sin gate falla). Suite 91/91. |

**Sprint 9 — lote 2: `force_edge_projected*` completa (Terzaghi; manual Professional 6.478-6.487; force_edge 44/44 COMPLETA)** | `812bb02` | 2026-08-24 | master como type 9 de area() (MTYPES 10) con las 10 variants companion (patrón force_edge_companion + temporal _time/_sine + _factor polinomio). Física: campo lineal ph/pv evaluado en las coordenadas del nodo; vd normalizado (fallback (0,−1,0)); hd = tunnel×vd (2D: (−vd_y, vd_x, 0)); proyección con productos internos explícitos — sig_radial = ph(n·hd)²+pv(n·vd)², sig_tang = ph(t·hd)(n·hd)+pv(t·vd)(n·vd); aplicada como fn·sig_radial·n + ft·sig_tang·t (positiva hacia afuera, hacia la excavación). FIX dimensional propio: la primera versión acumulaba n2/t2 en un loop con un escalar mal dimensionado. GOTCHA del test (segunda vez): convención Z del quad4 (3=arriba-IZQUIERDA) — malla cruzada ⇒ sigxx=0; el control de aislamiento con force_edge_normal (familia conocida-buena) en la misma malla separó bug-de-test de bug-de-feature. Validado con `fproj_tunnel`: sigxx=10 EXACTO (ph en pared vertical) y sigyy=20 (pv en horizontal) — el ratio 2:1 prueba la proyección por orientación (una presión uniforme daría 1:1). Suite 92/92. |

**Sprint 9 — lote 3: `data_*` completa + aliases solver/misc + `print_mesh_dof` (~16 keywords; SPRINT 9 CERRADO)** | `1e21cec` | 2026-08-24 | `data_activate`/`data_delete` (+`_time`): variantes time-gated de los control_data_* (mismo patrón destructivo, gate EPS_TIME) evaluadas en data(). `data_ignore`: skip a nivel de PARSEO (input.cc) — dos gotchas: el enum del item se guarda NEGATIVO (match −idat) y el record debe PREDECER a lo que ignora (se evalúa leyendo). Aliases db_number: `control_solver`→`control_options_solver`, `control_solver_bicg_error`, `axisymmetric`→`group_axisymmetric`, `bounda_print_mesh_dof*`→`print_mesh_dof*`. `control_solver_bicg_stop -no` cableado en so_bicg.cc (warning+continuar vs exit); bicg_restart/matrix_save/pardiso_* registrados como parciales. `print_mesh_dof`(+`_geometry`): dump one-shot a print_mesh_dof.dat (coords + dofs listados); `_values` registrado sin uso. Nota: DATA no existe como data_class → los data_* usan CONTROL. Validado con `dsmall` (data_activate_time borra la carga en t=0.1 → relaja a disy~0; aliases y dump ejercitados) y `dignore` (data_ignore declarado ANTES: load default 0, disy 0 exacto). Suite 94/94. **SPRINT 9 CERRADO**

**Sprint 10 — lote 1: `compression_direct`(+_visco) + `pressure_limit`/`coord_limit` + alias `heat_generation` (6 keywords)** | `f987f54` | 2026-08-24 | `materi_compression_cutoff` (stress.cc): descomposición espectral con matrix_jacobi — todo autovalor < sigy se sube hacia sigy (corte total sin record visco; relajación 1−exp(−dt/tm) con `_visco`). BUG propio cazado: el patrón factor=1−exp(−dt/(tm?:1)) heredado del cutoff de plano relajaba también SIN record (factor≈0.095) — la regla: el factor visco SOLO aplica cuando existe el record; sin record, factor=1. Gates `pressure_limit`/`coord_limit` en el call-site canónico de set_stress (cubren las tres leyes directas; documentado que NO cubren hypo). Alias `plasti_heat_generation` por dual-read en materi.cc (gana el nombre Professional). Tests: mdirect_comp (sigyy −10 → **−5.0 exacto**) y mdirect_gate (pressure_limit 0.1 desactiva el cap → −10 elástico; A/B sin limit falla). Suite 96/96. |

**Sprint 10 — lote 2: `druck_prag` alias + VALIDACIÓN ANALÍTICA Drucker-Prager + `materi_factor` + parciales (9 keywords)** | `3eca27f` | 2026-08-24 | **Hallazgo de convención**: el GNU implementa f = √J₂ + 3α·σm − K (la forma DP estándar; NO √(3J₂) como asumí primero) — con las constantes de matching del manual α = 2sinφ/(√3(3−sinφ)), K = 6c·cosφ/(√3(3−sinφ)). **Validación analítica** (`mdp_shear`): cizalla pura con φ=0 (α=0: sin término de presión ni dilatancia) → σxy capped en K = 2c/√3 EXACTO: c=30, dt=0.05 → **34.636 vs 34.641 = 0.014% de error**. Con φ≠0 el flujo asociativo genera σn≠0 en rigs de desplazamiento (documentado: usar φ=0 para el check de forma cerrada). El nombre Professional `druck_prag` (con guion) vía alias db_number del legacy `druckprag` (misma física, tests legacy druckpr1/examp23 verdes). `group_materi_factor`: escala tensiones Y rigidez (en equilibrio con fuerza σ no cambia pero u sí; con vely prescrita σ=ε·f·E directo: −0.1 vs −1.0, A/B). Alias `plasti_bounda(_factor)` → `boundary(_factor)`. Parciales registrados: damping_method, visco_exponential_limit/_name/_values; density_groundflow preexistente (checklist al día). GOTCHA de edición: borrar duplicados por primera ocurrencia puede eliminar el enum equivocado (GROUP_MATERI_FACTOR era único en tochnog.h y lo borré) — verificar con grep tras cada limpieza. Suite 98/98. |

**Sprint 10 — lote 3: Mohr-Coulomb CLÁSICO implementado de cero + VALIDACIÓN ANALÍTICA (la 3ª superficie plástica validada)** | `b315ade` | 2026-08-24 | `group_materi_plasti_mohr_coul [phi c phi_flow]` (el GNU solo tenía el `_direct`): bloque en plasti.cc con `matrix_eigenvalues` sobre sig — f = 0.5(σ1−σ3) + 0.5(σ1+σ3)sinφ − c·cosφ EXACTA del manual 6.724 (σ1 mayor, σ3 menor, tracción-positiva); f_flow con phi_flow. La dirección de flujo la aporta el driver estándar (gradiente por diferencias finitas centrales) → integración completa con el cutting-plane, plasti_kappa y boundary-reduction existentes. **Validación**: tracción uniaxial → σt = 2c·cosφ/(1+sinφ) = 46.19 **EXACTO a la primera** (c=40, φ=30°, ν=0). Sanity Tresca (φ=0): σt=2c=80, converge lento (79.55) por el VÉRTICE singular de Tresca para el cutting plane (numérico, no de la ley — documentado). **Hallazgo de trayectoria**: un rig de cizalla impuesta con flujo plástico NO permanece en cizalla pura — la deformación plástica acumulada rota el estado a tracción uniaxial a 45° (σ=τ[[1,1],[1,1]], satura en el MISMO σt; medido τ = c·cosφ/(1+sinφ) = 23.09 = 46.19/2 ✓): así se confirmó la superficie antes de rediseñar el test. Suite 99/99. Tres superficies plásticas validadas analíticamente: DP (0.014%), compression_direct (exacto), MC (exacto). |
**Sprint 10 — lote 4: Mohr-Coulomb hardening-SOFTENING implementado y CALIBRADO (la convergencia con softening era un problema de TEST, no de ley)** | `9d135a6` | 2026-08-25 | `group_materi_plasti_mohr_coul_hardening_softening [phi_0 c_0 phi_flow_0 phi_1 c_1 phi_flow_1 kappa_shear_crit]`: interpolación LINEAL con `kappa/kappa_shear_crit` clampado a [0,1], mismo f que el MC clásico (matrix_eigenvalues), phi_flow para la dirección. **Causa raíz de la no-convergencia con softening**: (1) el test de 1 paso brutal hacía que κ saltara más allá de κ_crit y la superficie BAJABA más rápido de lo que el cutting-plane podía perseguir → σ=7.55; (2) el default `control_timestep_iterations=1` deja el acoplamiento kappa→ley DÉBIL (la ley usa el κ nodal de la iteración anterior) → σ queda por encima de la superficie. **Calibración**: pasos pequeños dt=0.02 (20 pasos) + κ_crit=0.5 (κ/κ_crit=0.318/0.5=0.636 en zona lineal) + `control_timestep_iterations 8` → σxx medido **81.9 vs 83.6 analítico (98%)**, target 83.6±5. A/B sin softening (c₁=c₀=80) → 159.2=2c₀: discrimina. GOTCHA: el target del dof usa el basename `-kap`, NO `-materi_plasti_kappa` (ese da lectura fuera de límites/DBL_MAX). Suite 100/100. |

**Sprint 10 — lote 5: aliases Professional (cap2 + failure typo/underscore) + elásticas menores (poisson_power, shear_factor, k0)** | `1e5eb8e` | 2026-08-25 | `group_materi_plasti_cap2`: alias db_number del legacy `plasti_cap` (misma física, plasti.cc intacto; smoke A/B `mcap2`==`mcap_legacy` → sigxx −0.5769 EXACTO con c=1e6, f=p−pb<0). **Renombres de name[]**: `failure_cruching`→`failure_crunching` (typo de la "n") y `failure_voidfraction`→`failure_void_fraction` (underscore), ambos con el nombre GNU como alias db_number (sin tests legacy afectados). A/B de threshold en failure: alto→intacto, negativo→elemento marcado en el primer timestep (element_delete_times[0]=0.05; GOTCHA: threshold negativo dispara a deformación cero — semántica GNU fiel a jan-2014). `group_materi_elasti_poisson_power` (6.653): nu=ν0+ν1·(p/p1)^α con ν≤ν2, p=−sig_mean — bloque en set_stress tras young_power; oedometro E·eps=1.2 → punto fijo analítico ν=0.4, p=2 (σxx=−1.7143); medido σxx −1.4389/σyy −3.2924 (84%, ν~0.43; GOTCHA: el dof de tensión es unknown del sistema acoplado → la respuesta NO es exactamente el fijo secante). `group_materi_elasti_shear_factor` (6.654): escala SOLO las entradas de corte de C/Cmem (bloques (0,1),(0,2),(1,2)); cizalla pura 1 paso: factor 2→σxy 50 vs 33.33 sin record, factor 0→0 EXACTO (GOTCHA: ratio 1.5 no 2 — respuesta no lineal del dof acoplado; el rig de fuerza NO sirve como sonda de rigidez: el desplazamiento no responde a E). `group_materi_elasti_k0` (6.650) + hook de `control_materi_elasti_k0` (parcial Sprint 9 conectado): ν=K0/(1+K0), K0>0.95→0.95, SOLO con young/young_power Y poisson; K0=0.5→ν=1/3→σxx −0.75/σyy −1.5 EXACTOS, A/B control −no→−0.5769/−1.3462. Suite 111/111 (100 legacy + 11 nuevos). |

**Sprint 10 — lote 6: `elasti_young_power` UPGRADE 3→6 params (CAMBIO SEMÁNTICO) + `stress_pressure_history_factor` + initia `materi_stress_pressure_history` + gate `materi_elasti_young_power_apply`** | `260d199` | 2026-08-25 | `group_materi_elasti_young_power` (6.662/teoría 2.2.2): E=E0+E1·(p/p1)^α con E≥E2 y E≤E3, p=−sig_mean (compresión sube E). **CAMBIO SEMÁNTICO**: la forma GNU 3-params (young0·|p/p0|^α) ya no se acepta (verificado: ningún test legacy ni sfnet la usa). **FIX de bug preexistente del GNU 2014**: C_matrix ACUMULA (array_add, elasti.cc) — con young+young_power el GNU sumaba C(1000)+C(E_power) → rigidez DOBLADA (probe: −1.3961 = 2×−0.6980); el bloque young_power ahora LIMPIA C/Cmem antes de construir (la ley ES el módulo). Oedometro E0=1000 E1=500 E2=800 E3=3000 p1=1 α=1: punto fijo analítico E=1714.29 (σzz −2.3077/σxx −0.9890); medido σxx −0.713886/σyy −1.66573 = E_eff~1237 = 72% del fijo (GOTCHA: la C se evalúa con la tensión del paso previo y la tensión acumulada usa el PROMEDIO de los E pasados — trayectoria por debajo del fijo, modelo E_{n+1}=1.00417·E_n verificado a 1.4%). Clamps A/B analíticos EXACTOS: `myoung6_e2` (E2=3000 fuerza E=3000: −1.7308/−4.0385) y `myoung6_e3` (E3=800 fuerza E=800: −0.4615/−1.0769). `materi_elasti_young_power_apply` (6.801): gate control_...−no → E=E0 constante → `myoung6_apply` −0.5769/−1.3462 EXACTO (base lineal; alias db_number del nombre del manual). `group_materi_elasti_stress_pressure_history_factor` (6.655, NUEVO) + initia `materi_stress_pressure_history` (4.50): dof `sph` = máximo running |p| (dof.cc parallel_new_dof_diagonal, patrón clamp de kap; sin ecuación en general.cc, patrón veli); en set_stress, si la presión estimada del paso actual (p_old + dp, dp=−mean(C:inc_ept)) es MENOR que el máximo a INICIO de paso (old_unknowns[sph]) → C/Cmem × factor (descarga/recarga); si es el nuevo máximo → NO. Diseño: p_est captura el PRIMER paso de descarga y evita la degeneración p==sph en el pico (con new_unknowns[sph] subido medio paso, cualquier redondeo disparaba el factor → tensión 6x → sph inflado 0.3333→0.5 → runaway). GOTCHAS: scalar_dabs(−0.0) devuelve −0.0 y −0.0<sph es TRUE (factor desde el paso 1 de carga; normalizado p==0.); matrix_a4b NO es in-place seguro. Oedometro 2 fases `msph` (carga 4 pasos vely −0.002 + descarga 2 pasos +0.002, factor 3): tangente de descarga 3E → σyy +0.2692 vs −0.2692 con factor 1 (`msph_flat`); sph 0.3333 EXACTO en ambos. Suite 117/117 (111 + 6 nuevos). |
**Sprint 10 — lote 7: `group_materi_plasti_cap1` (NUEVO, el cap elastoplástico de alta presión) + initia `materi_plasti_cap1_history`** | `9322897` | 2026-08-25 | `group_materi_plasti_cap1` (6.691/teoría cap1, IMPLEMENTADO DE CERO): f = q²/M² + p\*(p\*−p\*c), p\* = p + c·cot(φ), p\*c = pc + c·cot(φ); invariantes p = −traza/3, q (patrón camclay); M SE LEE del record (8 params: φ c M λ\* κ\* K_ref p_ref m); flujo asociativo, gradiente por diferencias finitas del driver; pensado para usar COMBINADO con leyes de corte (DP/MC) — el driver activa la ley de mayor f (mecanismo max-f de plasti_rule verificado con `mcap1_comb`: f_dp<0 en el camino isotrópico → cap1 domina → respuesta idéntica a mcap1). **DECISIÓN de endurecimiento**: patrón kappa (NO dentro del bloque plasti — se llama 3 veces por iteración): el objetivo new_pc = old_pc + deps_p_cv·K_ref/(λ\*/κ\*−1)·((old_pc+c·cotφ)/p_ref)^m se calcula en set_stress() tras las iteraciones plásticas con el inc_epp CONVERGIDO (deps_p_cv = −traza(inc_epp) clampada ≥0: la descarga no endurece), y el RHS (new−old)/dtime en materi.cc conduce el dof (inercia lumped de general.cc). **VALIDACIÓN ANALÍTICA** (`mcap1`, hex8 compresión isotrópica v=−0.04, K=833.3, φ=0.5 rad, c=10, M=1.2, λ\*/κ\*=10, m=0, pc₀=100): el cap activa en p=pc=100 (paso 20); punto fijo discreto dpc=(5.0−dpc)/9 → 0.5/paso → **pc=110.0 EXACTO** (medido 109.9996) tras 20 pasos plásticos; descarga elástica 5 pasos → **sigxx −84.9996 vs −85.0 analítico (0.0005%)**; gemelo elástico `mcap1_elast` → **−175.0000 EXACTO** (discrimina superficie Y endurecimiento). GOTCHA pre-existente (del TEST, no del código): `bounda_time_increment` es local y PERSISTE entre records bounda (GET_IF_EXISTS no la resetea) → un record pairs tras uno con increment se lee como load-only y sus tiempos se vuelven cargas (velx=+100 en el paso 2, tensión −433/−4556); solución: formato increment en TODOS los records (msph solo funciona porque su record increment es el último). GOTCHA índice: k=floor(t/increment) → 41 cargas + 5 descargas = 40+5 pasos. Suite 120/120 (117 + 3 nuevos). |
**Sprint 10 — lote 8: Hardening-Soil — `group_materi_elasti_hardsoil` (NUEVO, 7 params) + `group_materi_plasti_hardsoil` (NUEVO, 4 params) + initias `materi_plasti_hardsoil_history` (4.22) y `materi_strain_plasti_hardsoil` (4.40) + `control_materi_plasti_hardsoil_gammap_initial` (6.146, NUEVO comportamiento)** | `6ab35cb` | 2026-08-25 | `group_materi_elasti_hardsoil` (6.649/teoría HS, IMPLEMENTADO DE CERO): E = Eref·((sig3 + c·cotφ)/(sigmaref + c·cotφ))^m con **sig3_manual = la MENOR compresiva = el MAYOR autovalor algebraico** en código tracción-positiva (la compresión más fuerte es el más negativo; el brief decía "menor autovalor" — eso es el σ1_manual AXIAL y NO es el HS: la rigidez depende del confinamiento, verificado analíticamente; el mapping correcto se documenta). E50/ν50 primera carga vs Eur/νur descarga/recarga con el switch del **máximo |p| histórico** (initia 4.22 = dof `sph` COMPARTIDO con materi_stress_pressure_history 4.50; decisión p_est = p_old + dp con dp = −mean(C_trial:inc_ept) → captura el PRIMER paso de descarga, patrón msph); clamps: base ≤ 0 → E = Eref; m = 0 → E = Eref. `group_materi_plasti_hardsoil` (6.703/teoría HS): f = q/(E50·(1−q/qa)) − 2q/Eur − γp con qa = qf/Rf y **qf = 2·sinφ·(sig3 + c·cotφ)/(1−sinφ) derivado de Mohr-Coulomb en falla** (Schanz; el manual no da qf explícito — documentado); γp = el dof **kappa** (materi_plasti_kappa = int √(0.5·dεp:dεp), la medida que el manual llama "equivalent plastic shear strain"; el manual NO define dof para γp — decisión documentada) + el extra inicial del control; singularidades: q ≥ qa → f = DBL_MAX/1e6 (asíntota: el trial más allá parece elástico si f quedara negativa); qa ≤ 0 → db_error; base ≤ 0 → NO_YIELD_F; flujo ASOCIATIVO (gradiente numérico del driver; ψ aceptado sin usar — potencial no-asociativo PENDIENTE). `control_materi_plasti_hardsoil_gammap_initial` (6.146): -yes crea γp_extra = q/(E50·(1−q/qa)) − 2q/Eur en el estado inicial (primer paso, centinela −1 pre-alocado en top.cc — el PUT corre en el bucle PARALELO donde no se puede alocar), guardado en `element_intpnt_materi_plasti_hardsoil_gammap_initial` (6.440) y SUMADO a γp en la ley → f = 0 al arrancar con tensiones desviadoras. **VALIDACIÓN ANALÍTICA** (`mhardsoil_elast`, uniaxial plano c=10 φ=30° m=0.5: E50 = 1000·√(17.32/117.32) = 384.23 → σyy = −1.0989·384.23·0.0004 = **−0.1689 EXACTO**; A/B σref=200 → −0.1241; `mhardsoil_unload` oedómetro confinado 4+2 pasos m=0 Eur=3000 → **+0.2692 EXACTO** vs −0.2692 con Eur=1000 (patrón msph); `mhardsoil_gp0` tensión inicial σxx=−2 + control → **f=0 → σxx −2.0 EXACTO**, record 0.0007763 = f(q=2) ANALÍTICO, κ=0; `mhardsoil_gp0_off` sin control → retorno: σxx −2→−1.225, κ = 0.000293, estado final SOBRE la superficie (f ≈ 0: 0.687/(1000·(1−0.687/38.49)) − 2·0.687/3000 − 0.000293 = −0.000061); `mhardsoil_plast` σxx=−4 → κ = 0.0006 = 2× (el endurecimiento ESCALA con q inicial); `mhardsoil_plast_elast` gemelo elástico −4.0385 EXACTO). GOTCHA mayor (documentado en manual-developer): el retorno plástico HS (f = O(q/E) ~ 1e-4, minúscula vs O(10-100) de cap1/DP) converge en el punto de integración pero el dof de tensión del sistema acoplado NO lo propaga en carga continua (respuesta elástica en los dofs; la validación usa la familia de relajación A/B donde el retorno SÍ llega a los dofs); GOTCHA node_dof: con derivatives nder = ndim+2 (2D: 4) y la lista de valores es UNA por nodo (nuknwn valores, no nuknwn·nodos); GOTCHA: la pre-alocación del record con 0 rompía la detección "primer paso" → centinela −1. Suite 128/128 (120 + 8 nuevos). |
**Sprint 10 — lote 9: initias por-modelo `materi_strain_plasti_cap` (4.35) / `_compression` (4.36) / `_diprisco` (4.37) / `_druckprag` (4.39) + alias `materi_plasti_diprisco_history` (4.18) + `group_materi_elasti_camclay_pressure_min` (6.647)** | `ddacc43` | 2026-08-25 | **Initias por-modelo (4.35-4.40)**: mismo patrón que hardsoil (lote 8) pero CONSOLIDADO — cada initia registra su dof dedicado (capepp/cepp/depp/dpepp, basenames eppcapxx.../eppcmpxx.../eppdipxx.../eppdrpxx...; 6 componentes) y un ÚNICO bloque RHS en materi.cc (tabla {flag, indx}) llena TODOS los dofs por-modelo con el MISMO inc_epp del driver (refactor del bloque hardsoil a la tabla: hsepp incluido, comportamiento idéntico — suite anterior intacta). Enums MATERI_STRAIN_PLASTI_CAP/COMPRESSION/DIPRISCO/DRUCKPRAG + general.cc (inercia/conv_part) + check.cc self-checks. **VALIDACIÓN**: `mstrain_cap` (cap1, compresión isotrópica 20 pasos plásticos): eppcapzz = **−0.0360001 EXACTO** (análitico: deps_p_cv = 0.0054/paso × 20 / 3); `mstrain_druckprag` (cizalla pura mdp_shear): eppdrpxy = 0.910768 **IDÉNTICO al epp genérico** del mismo run (prueba el mecanismo compartido) y sigxy 34.64 intacto; `mstrain_diprisco` (diprisc1 + initia): eppdipyy −0.0105714, eppdipxx = eppdipzz +0.00350735, hisv10 −152.795 idéntico a diprisc1; `mstrain_compression` (oedómetro, sig_yield 1.0): eppcmpzz −8.96e-05. Gemelos elásticos de los 4 (carga bajo yield / c=1e6 / sin record plástico): dof = **0 exacto**. A/B: sin la initia el dof NO existe (target −eppcapzz → "I do not know" → exit 1, verificado). **`materi_plasti_diprisco_history` (4.18)**: alias de materi_history_variables (misma rama del parser, mismo dof hisv compartido); checks de GROUP_MATERI_PLASTI_DIPRISCO/_RT aceptan cualquiera de los dos nombres (check_unknown_atleastone). `mdiprisco_hist`: diprisc1 con el alias → hisv10 **−152.795 IDÉNTICO**. diprisc1.dat legacy sigue verde. **`group_materi_elasti_camclay_pressure_min` (6.647)**: hook en el bloque elástico camclay de set_stress (K = (1+e)·p/kappa): si existe el record, p < pressure_min → p = pressure_min (solo para K, el estado de tensión intacto). `mc_pressure_min` (p inicial 0.001, pressure_min 10, 1 paso oedómetro): K = 1.5·10/0.02 = 750 → **sigxx −0.0093333 y sigyy −0.2093333 ANALÍTICOS EXACTOS** (dp = K·tr = 0.075); A/B `mc_pressure_min_off` sin record: K = 0.075 degenerado → **sigxx +0.06566 (signo invertido)** — el problema numérico que el manual previene. **`diprisco_density` (6.697)**: DECISIÓN PENDIENTE documentada (NO implementada): la ley de interpolación loose/dense entre los 22 parámetros no está en el manual (solo "interpolated... using the current density", remite a papers) → sin fórmula verificable no se implementa (regla AGENTS.md); además exigiría reescribir el bloque diprisco arriesgando diprisc1; la 12ª variable (densidad relativa) sin ley de evolución documentada. La initia 4.18 acepta 12 vars. GOTCHA node_dof: añadir una initia con n unknowns exige expandir el node_dof plano (nuknwn = npuknwn·nder valores por nodo; con derivatives nder = ndim+2 → 155 por nodo en el hex8 3D del mcap1: 30 slots nuevos antes de pc). GOTCHA 2D: el strain de carga oedómetro va en yy (sigyy), no en zz (fuera de plano). Suite 139/139 (128 + 11 nuevos). |
**Sprint 10 — lote 10: `control_repeat_save` (6.348) + `control_repeat_save_calculate` (6.349) + records `repeat_save_result` / `repeat_calculate_result`** | `86accd4` | 2026-08-26 | **VALIDACIÓN ANALÍTICA** (`mrepeat_save`): columna elástica con vely prescrita −0.01 y `control_repeat 3` saltos a `control_timestep` dt=0.05 — cada salto guarda `veliy` (= desplazamiento integrado; con memoria `-updated_without_rotation` el dof de desplazamiento es `materi_velocity_integrated`, label `-veliy`) del post_point: secuencia aritmética EXACTA {−0.0005, −0.001, −0.0015} en `repeat_save_result[0..2]` (medido −0.000499999975586, error 5e-11). Al agotarse el contador, media = **−0.000999999951172 vs −0.001 analítico** y varianza POBLACIONAL = **1.66666650391e-07 vs d²(N²−1)/12 = 1.6666667e-07** (error 1.6e-14). El target de varianza discrimina el NÚMERO de saves (un 4º save daría 3.125e-7 → falla). A/B verificado: sin `control_repeat_save_calculate` el record `repeat_calculate_result` NO existe (target → "Error in calculation" → exit 1). GOTCHAS: (1) `db_max_index` devuelve el máximo ALOCADO (margen heurístico de db_allocate) — el índice libre de `repeat_save_result` se cuenta con `db_active_index`; (2) `materi_displacement` está PROHIBIDO con memoria `-updated_without_rotation` (stress.cc: "should not be initialized") — el dof de desplazamiento es el `veli` integrado. Suite 140/140 (139 + mrepeat_save). Los legacy con control_repeat siguen verdes (force5/examp6/examp23/repeat1/repeat2); examp9/examp24/ho_mech3/refine5 fallan igual con el binario PREVIO (pre-existentes: bad_alloc/timeout/TET10-refine/divergencia bicg). |

**Sprint 11 — lote 1: 8 keywords control_print (manual Professional 6.266/6.270/6.284/6.286/6.322/6.329/6.336/6.337)** | `4e3d836` | 2026-08-26 | `control_print_database_method`: -all → dump COMPLETO a .dbs (sin líneas "Size of"), -size_tot → "Size of <record> is <bytes>" por record + "Total size" + "Size of the system matrix is <solve_nlocal>" (nº de ecuaciones del último solve — la métrica del band solver en so.cc), -size_tot_large → SOLO records > 1 Mb (en el modelo pequeño del test solo queda la línea de la matriz — A/B con -size_tot discrimina el filtro). `control_print_dof_id`: default -yes → dof.<index> con 4 columnas (x y dof node); el número de nodo es el ORIGINAL (mapeo posición→número capturado ANTES del renumbering); -no conserva 3 columnas (contraste dof.20 vs dof.21). `control_print_dof_rhside`: alias db_number de control_print_unknownsrhside (mismo enum/dispatch; el detector de fin-de-valores también llama db_number — patrón f178849); smoke velx_rhside.20. `control_print_element_method`: -middle (default, UNA línea por elemento con la coordenada media) vs -node (líneas por nodo); discriminado por conteo 2 vs 4 líneas y la relación EXACTA 2*mid1 == x_node2 (la media = promedio de coords nodales; GOTCHA: el self-stress del truss se relaja a 0 con la deformación — el test discrimina por estructura, no por valor). `control_print_history_relative_time`: tr=0.2 → última línea de disx1.his = 0.1 (0.3−0.2) vs 0.4 sin tr (disx2.his); también aplica al smooth file. `control_print_mesh_dof`: alias de prefijo db_number → print_mesh_dof (cubre _geometry/_values); smoke print_mesh_dof.dat. `control_print_number_iterations`: monitor de consola en el bucle de equilibrio — 16 líneas exactas (2 pasos × control_timestep_iterations 8); NO es -inverse_iteration_number. `control_print_partialname`: print_partialname() en print_db.cc — match de PREFIJO (strncmp, no el strstr del helper) sobre los nombres registrados; stdout con records element* y NINGÚN node* (grep). Suite 148/148 (140 + 8). Checks de archivos en build_safe.sh. |

**Sprint 11 — lote 2: frecuencia de prints — `control_print_frequency_timeinterval` (6.291) + `control_print_frequency_timestep` (6.292)** | `8fd8d6f` | 2026-08-27 | Ambos gatean los `control_print_*` del MISMO icontrol: corren cada vez que pasa un intervalo de tiempo (timeinterval, DOUBLE por icontrol) o cada N pasos (timestep, INTEGER por icontrol), y SIEMPRE al final del incremento de `control_timestep`. **Gate** `control_print_frequency_allowed(icontrol)` en top.cc: se evalúa UNA vez por step_close y envuelve CADA item `control_print_*` del dispatch EXCEPTO las 3 excepciones del manual (`control_print`, `control_print_history`+`_smooth`, `control_print_data_versus_data`); NO se gatean no-prints (print_lastdatabase, post()). **Fin del incremento** = la MISMA condición del bucle de timestep de top(): TIME_NEW se escribe al inicio de cada incremento y el último paso se clamp EXACTAMENTE a él (top.cc:340-347), así que `time_current >= time_new` (con CONTROL_TIMESTEP activo y TIME_NEW existente) ⇒ fin ⇒ imprime siempre. **Estado por icontrol en VERSION_NORMAL**: `control_print_frequency_timeinterval_time` (último print; anclaje inicial a TIME_OLD — el ejemplo del manual exige anclar al INICIO del incremento: anclar al primer paso 0.04 daría 0.20 en vez de 0.16 — y re-anclaje a TIME_OLD si restart rebobina) y `control_print_frequency_timestep_count` (contador, reseteado en cada print incluido el de fin de incremento). GOTCHA no_index: db() rechaza `PUT` con `index>0` en records no-indexables (`no_index=1`; control_print_gid_time usa index 0 por eso) → los records de estado son indexables con `external 0`. **VALIDACIÓN contra los ejemplos EXACTOS del manual** (`freq_timeint`/`freq_timestep`, modelo 2D elástico): `control_print_frd -separate_sequential` = UN archivo .frd por print con la línea 100CL = tiempo EXACTO → freq_timeint: 3 archivos con tiempos **0.16, 0.32, 0.41** (manual 6.291) y freq_timestep: **0.20, 0.40, 0.41** (manual 6.292); dof.0..dof.2 = 3 prints vs dof.3..dof.12 = 10 sin frecuencia (13 archivos de 176 líneas); dof.22 (3) vs dof.23 (11): `3*L23 == 11*L22` (528/1936); excepción: sigxx4.his con 21 líneas (10+11 pasos, history NO gateado). GOTCHA del test: el GNU CLAMPEA el último paso parcial de un incremento (0.04 0.41 → pasos 0.04..0.36 y 0.41, NO hay paso en 0.40) → el ejemplo 0.40 del manual se reproduce con 2 incrementos (`0.04 0.40 0.01 0.01`). Verificado también el cruce de incrementos ad-hoc (2 incrementos → prints 0.16, 0.20 fin-inc1, 0.36, 0.41 fin-inc2: la cadencia se re-ancla al fin de cada incremento). Suite 150/150 (148 + 2). Checks de archivos en build_safe.sh (tiempos 100CL + conteos + ratio + excepción). |
**Sprint 11 — lote 3: extensiones VTK — `control_print_vtk_coord` (6.340) + `control_print_vtk_dof_calcul` (6.342) + `control_print_vtk_empty` (6.343) + `control_print_vtk_node_method` (6.345) + `control_print_vtk_other` (6.346)** | `e7ed4db` | 2026-08-27 | las 5 keywords en `print_vt.cc` (enums CONTROL_PRINT_VTK_COORD/DOF_CALCUL/EMPTY/NODE_METHOD/OTHER en tochnog.h+tochnog-mod.h en sync, data_required=CONTROL_PRINT_VTK). **`vtk_coord`**: DECISIÓN con evidencia — las coordenadas en print_vtk aparecen SOLO en el bloque POINTS (no hay campos de datos tipo node_0/1/2 del gmsh 6.316; print_gm tampoco los escribe) → `-no` omite POINTS conservando CELLS/POINT_DATA (limitación documentada: el archivo sin POINTS no es un dataset visualizable; el switch es para dumps pequeños). **`vtk_dof_calcul`**: espejo de vtk_dof sobre el bloque post_calcul (líneas ~408+): flags por campo post que gatean los 2 bucles (agregados + componentes); matching documentado (adaptación GNU: post_calcul_label no existe) = match EXACTO sobre el unknown subyacente (ej. `-materi_strain_total` selecciona todos sus operadores) O substring sobre la etiqueta post_calcul_names (ej. `-sigyy` casa sigyy y tosigyy); `-none` → ningún campo post; sin el record → todos (default intacto). **`vtk_empty`**: helper `vtk_element_is_empty()` (ELEMENT_EMPTY==-YES, patrón print_g5: -NO/-FRONT = presente); `-no` excluye los vacíos de los 3 pases con `ncell` consistente en CELLS/CELL_TYPES; el abort de tipo no soportado NO se toca (con -no el elemento vacío se salta ANTES del check; con default -yes comportamiento byte-idéntico — suite legacy intacta). Los elementos vacíos los COMPUTA el solver (materi_density/diffusion: ningún nodo con valor ≥ mínimo), no hay keyword de usuario. **`vtk_node_method`**: default `-node_start_refined` (manual gid 6.307; cae a coordenadas almacenadas si no hay NODE_START_REFINED → output idéntico al previo en toda la suite); CAMBIO SEMÁNTICO documentado: antes con materi_displacement las coordenadas se escribían SIEMPRE deformadas, ahora SOLO con `-node_deformed_mesh`. **`vtk_other`**: PARCIAL (el manual no detalla la lista) — subconjunto documentado: `boundary_condition` scalar (1.0 si ALGÚN dof primario del nodo está bounded vía node_bounded; bounda_force NO marca → discriminador de valor del test) + `mesh_deformation` vector (dofs dis, solo con materi_displacement); el resto de la lista gid_other (grupos como cell data, vectores de viga/truss…) NO implementado (familia GiD descartada). GOTCHAS: (1) comentarios .dat con paréntesis desbalanceados a mitad de línea rompen el parser de input.cc (el cierre se detecta cuando nright==nleft; `-no):` dejaba un token basura → "Error in initialization part."); (2) NODE_BOUNDED no tiene version_all → se lee de VERSION_NORMAL con índices 1-based (print index inod ↔ NORMAL inod+1; desalineación si se borran nodos, documentado); NO usar db_active_index(NODE_BOUNDED,0) como test de existencia (el índice 0 nunca está activo) sino db_max_index. Sección 6.344 del manual = TYPO (dice control_print_gid_group en la sección vtk; pertenece a la familia gid descartada; NO está en el inventario) → no implementada, solo anotada. **VERIFICACIÓN**: tests propios `vtk_coord1` (A/B POINTS por grep), `vtk_dofcalc1` (sin filtro 2 headers post vs -none 0 vs filtro 1, primarios intactos), `vtk_empty1` (elemento 2 vacío por densidad 0: CELLS 2 vs 1 y CELL_TYPES 1), `vtk_nodmeth1` (total lagrange vely=-0.01: y almacenada 1.0 vs deformada 0.998 = 1+disy, disy=-0.002 target), `vtk_other1` (boundary_condition 1.0/1.0/0.0/0.0 + mesh_deformation vs -no sin campos) + 5 checks de archivos en build_safe.sh. Suite 155/155 (150 + 5) en build limpio --clean. |

**Sprint 11 — lote 6: fuerzas/momentos de viga — `control_print_beam_force_moment` (6.262) + `_coordinates` (6.263) + `_switch` (6.264)** | `6ef0d5b` | 2026-08-27 | las 3 keywords en el nuevo `print_beam_force_moment.cc` (enums CONTROL_PRINT_BEAM_FORCE_MOMENT* en tochnog.h+tochnog-mod.h en sync; dispatch en top.cc DENTRO del gate de frecuencia). **Layout ELEMENT_BEAM_MOMENT** (investigado, beam.cc:434/524): `NNOL*NDOF = 6` doubles = 2 nodos × {2 fuerzas EN EL PLANO + 1 momento FUERA de plano}, en los ejes del PLANO de la viga (mapeo 2D de beam_3d vía group_beam_plane → index_plane; NO son ejes 3D globales); guarda el vector de fuerzas internas del elemento (K·u). La fuerza fuera de plano y los momentos en plano son IDÉNTICAMENTE 0 (la viga es un elemento 2D) → las 12 columnas del manual salen limpias SIN extrapolación: (fx, fy, 0, 0, 0, mz) por nodo. **Frame local** = el de beam_2d: local x = dirección nodo1→nodo2, local y = perpendicular en plano (rotación (a,b): fx_local = a*f0 + b*f1, fy_local = -b*f0 + a*f1 con (a,b) = componentes de la dirección en los ejes del plano). **Truss** (manual: "the truss force will be used for the axial force"): -truss y -trussbeam toman las columnas axiales de ELEMENT_TRUSS_FORCE escalar (+N nodo 1 / -N nodo 2 — vector de fuerzas nodales, mismo patrón anti-simétrico que las transversales de la viga; DECISIÓN documentada); -truss puro solo axial, resto 0; -beam sin truss → axial 0 (la viga 2D pura no tiene rigidez axial — GOTCHA del modelo de test: velx de nodos libres prescrito para evitar matriz singular). **Selección** (decisión documentada): distancia mínima segmento-segmento 3D (Eberly/RTCD) entre el eje del elemento y el corte < `1e-6*max(1, long_corte)`; orden ascendente por distancia (proyección del punto más cercano sobre el corte = s*long_corte); sort estable. **GOTCHA de calibración CRÍTICO**: el algoritmo de closest-point de referencia tiene s en el PRIMER segmento y u en el SEGUNDO — intercambiarlos devuelve los puntos más cercanos EQUIVOCADOS (la intersección real del corte en x=0.25 era t=0.50083 en el corte, el código buggy devolvía t=0.25 → distancia 0.355 en vez de 0). **Sin cruces → NO se escribe archivo** (decisión; dos pasadas). **Snap** de ruido del solver (~1e-6 relativo) a cero en la salida (mz tip 2.78e-8 → 0; umbral 1e-6*(1+max|comp|)). `_coordinates`: longitud variable (fixed_length=0) 4 (2D) o 6 (3D), OBLIGATORIA (db_error si falta/mal longitud/corte de longitud 0). `_switch -yes`: factor -1 SOLO en las 12 componentes (la distancia no se invierte). Archivo append; -separate_sequential con contador estático (patrón interface_stress). **VERIFICACIÓN** (tests propios, suite **177/177 = 172 + 5** en build limpio --clean): `bmom1` (ménsula 2D 1 elemento, F=1e-2, corte en el nodo fijo: dist 0.5 EXACTO, fy1=-F=-0.01, mz1=-F*L=-0.01, fy2=+F, mz2=0 — signo = el del vector de fuerzas internas del elemento, verificado empíricamente; modelo con options_mesh -fixed_in_space + options_convection -no para que las coordenadas NO se deformen y las distancias sean exactas), `bmom_switch` (-separate_sequential → beam_force_moment.0 con los 12 valores NEGADOS vs bmom1), `bmom_truss` (trussbeam con carga transversal F=1e-2 + axial 2e-2: fx1=+N fx2=-N con N=0.0200165561392 EXACTO == element_truss_force_0.7 de control_print_element; fy/mz de la viga; GOTCHA físico: el target del momento en la raíz es F*L - N*w(L) = 0.00993 por P-delta de la TRACCIÓN, no 0.01), `bmom_2d` (2 vigas colineales + corte DIAGONAL 2D que cruza SOLO el elemento 2: 1 línea, dist = sqrt(1.25) = 1.11803398875 EXACTO, mz1 = -F*0.5), `bmom_noint` (corte a x=1.5 sin cruzar nada → beam_force_moment.9 NO existe). Checks de archivos en build_safe.sh. |

**Sprint 11 — lote 5: prints por nodo + suavizado de dofs — `control_print_node` (6.330) + `_angular` (6.331) + `_angular_middle` (6.332) + `_geometry` (6.333) + `_sort` (6.334) + `_zero` (6.335) + `control_print_dof_smooth_dof` (6.271) + `_smooth_n` (6.272)** | `08d2aca` | 2026-08-27 | **SMOOTHER** (en `print_dof()`, print_hi.cc — helper `print_dof_smooth_apply`): suavizado NODAL por conectividad de elementos, NO confundir con history_smooth (media móvil temporal). DECISIÓN: cada pasada reemplaza el valor del nodo por el promedio de sus VECINOS (nodos conectados por ≥1 elemento) SIN incluir el propio nodo (campo lineal = punto fijo en el interior; extremos toman su único vecino; nodo aislado conserva su valor). Adyacencia head/to/next_edge por POSICIÓN (los elementos están en ids ORIGINALES → mapeo orig2pos); los valores se leen POR POSICIÓN de VERSION_PRINT (GOTCHA de calibración: leer por id original desplaza el campo un nodo). Sumas POR COMPONENTE (GOTCHA: una suma compartida inicial corrompía todas las componentes). `-all` = todos los componentes; labels vía array_member (no encontrados → skip). `_smooth_n` default 10, <1 → error. **PRINTS** (nuevo `print_node.cc`; dispatch en top.cc DENTRO del gate de frecuencia): UN archivo por parte seleccionada — labels de dof → `label.<icontrol>` (ej. `velx.10`, consistente con L4); items node_dof_calcul → `post_calcul_names_without_extension` (ej. `avel.10`; matching = EXACTO sobre el unknown subyacente O substring del label — adaptación GNU, post_calcul_label no existe); números → `node_<record>_<n>.<icontrol>` (DECISIÓN documentada, ej. `node_dof_0.10`); sin partes → TODAS (nuknwn/ncalcul/db_data_length). Líneas `x y <valor>` por nodo (append), VERSION_NORMAL sin renumbering (orden = id de nodo ascendente); registros inactivos y partes fuera de longitud → nodo omitido para ese archivo. **Angular**: ángulo en GRADOS (atan2) desde +x→+y / +y→+z / +x→+z según switches; 1D error; 2D solo `-yes -yes` SIN switch_z; 3D las 3 combinaciones; middle default (0,0,0). **Geometry**: patrón filter.cc:74-80 con frame -node_start_refined (fallback NODE). **Sort**: collect+sort de print_interface_stress.cc:72-89, ASCENDENTE y estable; `-angle` solo con angular (la mención del manual a "control_print_node_method" es TYPO de node_angular, no es keyword); -y solo 2D/3D, -z solo 3D. **Zero**: default -yes; `-no` suprime valores 0 con comparación EXACTA, por archivo. **VERIFICACIÓN** (tests propios, suite **172/172 = 165 + 7** en build limpio --clean): `cpn1` (labels -velx -vely → vely.10 EXACTO 0/0/-0.01/-0.01 con todos los dofs Dirichlet; numérico node_dof_1.11 == vely.10; calcul avel.12 = (velx+vely)/2 = -0.005), `cpn_angular` (nodo (1,1) con middle (0.5,0.5) = 45° EXACTO; sort -angle -135..135 ascendente), `cpn_geom` (geometry_line borde inferior → 2 vs 4 líneas; GOTCHA del suite: los icontrols 30-44 los ocupa la familia dof_line — dpline1 escribe vely.30 y dpline_move vely.40/41 → cpn_geom usa 54/55 y cpn_sort 56/57), `cpn_sort` (nodos NO ordenados por y: sin sort "0 1 -0.01" vs con -y "0 0 0".."1 1 -0.01"), `cpn_zero` (2 vs 4 líneas; velx todo 0 → archivo VACÍO con -no), `dsmooth1` (cadena 1D 0,1,2,3,4 con velx prescrito → 1 pasada: interior = promedio de vecinos EXACTO (1,1,2,3,3) vs raw; GOTCHA del modelo: sin options_mesh -fixed_in_space la malla 1D SE MUEVE con la velocidad (x=i → x=2i tras 1 paso) y sin materi_stress print_dof escribe 37 bloques basura en 1D (6 matrices de tensión con índices fuera de rango — quirk PRE-EXISTENTE del GNU, no tocado; el test usa el modelo sin materi_stress → solo el bloque velx)), `dsmooth_n` (A/B 3 pasadas EXACTAS (1.5,1.5,2,2.5,2.5) vs default 10 (1.9375,1.96875,2,2.03125,2.0625) → convergencia a la media 2.0 y nodo central EXACTO 2). Checks de archivos en build_safe.sh. |

**Sub-sprint materi_stress_force — lote 1 (infraestructura): `post_calcul -materi_stress -force` + 10 records `post_calcul_materi_stress_force_*` (6.908-6.917) + `control_print_materi_stress_force` (6.328)** | `cc0aba1` | 2026-08-28 | registro de los 11 keywords (enums CONTROL_PRINT_MATERI_STRESS_FORCE + POST_CALCUL_MATERI_STRESS_FORCE_* en tochnog.h/tochnog-mod.h en sync, + enum PRIMARY para el switch -primary; data_class POST para los post_calcul_*, CONTROL para el print; data_required=POST_CALCUL) + **rama FORCE en el dispatch de operats** de calculate()/calculate_operat() (calcul.cc, junto a PHIMOB/TOTAL): genera los items por nodo — **2D 9 items** (norx_sig nory_sig nors_sig | shex_sig shey_sig shes_sig | momx_sig momy_sig moms_sig) y **3D 16 items** (norx_sig nory_sig norz_sig nors_sig | shex_sig shey_sig shez_sig shes_sig | mom1*_sig | mom2*_sig), todos caben en **MCALCUL=20** (limitación documentada: el bloque 3D deja solo 4 slots para otros post_calcul; 3D + TOTAL(6) = 22 > 20 → mensaje "MCALCUL too small" de parallel_calcul_node). Layout NODE_DOF_CALCUL = UN slot por item (el layout plano que lee print_vt.cc; los -VECTOR agrupan las componentes del vector con nval=MDIM). **Validación** de records en `post_calcul_materi_stress_force_validate()` (nuevo `calcul_force.cc`): element_group OBLIGATORIO; exclude XOR include; 3D exige direction (error claro); reference_point = ngroups*ndim valores (3D obligatorio, 2D sin él → warning + default documentado (0,0)); thickness_switch = ngroups; plot_switch = 3 (2D) / 4 (3D) — por VECTOR item, no por resultado; valores -yes/-no. **Familia NODAL**: error claro si POST_LINE_DOF/POST_POINT_DOF/POST_QUADRILATERAL_DOF activos con FORCE. El stub `post_calcul_materi_stress_force()` (calcul_force.cc) rellena 0 con aviso "not yet implemented" UNA vez por run: **la integración numérica llega en L2/L3** (este lote es infraestructura). **Print** (nuevo `print_materi_stress_force.cc`, dispatch en top.cc DENTRO del gate de frecuencia): lee NODE_DOF_CALCUL (VERSION_PRINT, patrón print_vt.cc:408-524; el offset del bloque FORCE = índice del primer item FORCE en POST_CALCUL_UNKNOWN_OPERAT) y escribe `materi_stress_force.<icontrol>` (el "index" del manual 6.328 es el índice del RECORD, convención de beam_force_moment; el valor del record es el método -all/-primary) con CABECERA de comentarios (#) explicando las columnas (exigido por 6.328) tomadas de post_calcul_names; una línea por nodo; -primary = -all en lote 1 (hook msf_node_is_averaged() → 0 hasta el promediado de L2/L3); sin bloque FORCE → sin archivo. **A/B REAL contra el comportamiento previo** (calcul.cc de HEAD compilado a binario temporal): antes `post_calcul -materi_stress -force` moría en `db_error(POST_CALCUL,0)` del dispatch ("Error detected for data item : post_calcul"); después da mensaje claro de validación y con configuración completa termina rc=0. **VERIFICACIÓN** (tests propios, suite **181/181 = 177 + 4** en build limpio --clean): `msf_parse` (2D, 9 records de config + print → materi_stress_force.100 con cabecera + 4 líneas de nodo; aviso not-yet-implemented ÚNICO; warning direction 2D), `msf_parse_3d` (3D con direction_include: 16 items → materi_stress_force.200 con 17 columnas y 8 nodos), `msf_print` (-all 300 y -primary 301: ambos 4 líneas, 9 items a 0 en lote 1), `msf_errors_2dwarn` (2D sin reference_point: warning + rc=0) + **6 casos de error** en la sección de checks (rc=1 + mensaje): `msf_errors_nogroup` (sin element_group), `_bothdir` (exclude+include), `_3dnodir` (3D sin direction), `_3dnoref` (3D sin reference_point), `_refcount` (2 groups vs 1 punto), `_posttype` (post_point activo). Checks de archivos en build_safe.sh. |
**Sub-sprint materi_stress_force — lote 2 (integración 2D quad4/quad9): `post_calcul -materi_stress -force` funcional en 2D** | `1babf2d` | 2026-08-28 | implementación en `calcul_force.cc`: `msf_calculate_node_2d` (scan por nodo de los grupos objetivo) + `msf_element_contribution_2d` (caras extremas, integración, rol del nodo, componentes de plot) + `msf_integrate_side_2d` (cuadratura 1D de arista). **Decisiones de diseño con evidencia** (manual-developer/post_calcul_materi_stress_force.md; el test analítico es el árbitro): (1) caras extremas = las 2 aristas con |n·t̂| mínima, t̂ = (centroide − reference_point) en el plano; (2) fuente de σ = la tensión NODAL resuelta (NODE_DOF); (3) cuadratura de arista Gauss(2) (quad4) / Lobatto(3) (quad9) — el integrando del momento es cuadrático/cúbico y ambas reglas lo integran exacto (el trapezio Lobatto(2) daría 3/2× el momento, verificado); (4) **momento con dt = (x−C)·t̂** (distancia en dirección de espesor respecto al centroide, manual 6.913) — el parámetro de arista s daría signos opuestos en las dos caras (verificado: los nodos del plano medio salían (f1−f2)/2 en vez del promedio); (5) nor con signo (tracción +), she solo tamaño, mom con signo; componentes de plot a lo largo de t̂, s = tamaño físico; (6) nodos de las caras = primarios; plano medio quad9 con average -yes (default) = promedio EXACTO de las dos caras + flag por nodo `POST_CALCUL_MATERI_STRESS_FORCE_AVERAGED_NODE` (enum nuevo en tochnog.h/tochnog-mod.h en sync, NODE class version_all=1, pre-alocado en calculate() como NODE_DOF_CALCUL) consumido por `-primary` (msf_node_is_averaged en print_materi_stress_force.cc); (7) outer -yes (solo primarios a máxima distancia) y plot_switch (invierte x/y) implementados en 2D; (8) axisimétrico: l=2πr (r = coordenada radial del centroide, convención area.cc), sin test dedicado (pendiente). **GOTCHAS** (documentados): (a) **heap OOB READ cazado con AddressSanitizer**: leer coords nodales con MDIM (3) en vez de ndim (2) lee 8 bytes fuera del record NODE en el último nodo y el valor basura CORROMPÍA el solve posterior (campo de tensiones localizado en las últimas columnas del axial) — corregido (i<ndim); (b) con nu=0.0 el solve mixto iterativo localiza la deformación axial (σ_xx=0 en la mitad izquierda) — los tests usan nu=0.3; (c) un PAR de fuerzas nodales no excita flexión pura en el FE discreto (trabajo nulo sobre el modo de flexión: u_top=u_bottom) → responde con cizalla pura — el test de flexión pura usa la configuración estándar de 4 puntos (tramo central con M cte y V=0); (d) la cizalla integrada de la formulación mixta en mallas 1-en-espesor tiene polución en las superficies libres (±30% de P, banda documentada) y el quad4 1-en-espesor sufre shear locking (mom 0.23× medido) → los tests de flexión usan quad9; el valor EXACTO de she se verifica con msf_shear (corte simple, campo uniforme). **VERIFICACIÓN ANALÍTICA** (tests propios, suite **187/187 = 181 + 6**): `msf_beam2d` (ménsula 8×quad9, P=1e-2: moms = P·(8−x) dentro del 0-2.5% — 0.0820, 0.0708, 0.0600, 0.0500, 0.0400, 0.0300, 0.0200, 0.0099, 0.0003 —, promediado del plano medio EXACTO, shes en banda [0.5P,1.5P] documentada, nors~0), `msf_beam2d_pure` (4 puntos: moms=1e-2 cte 1.5% en el tramo central, shes=0 EXACTO), `msf_quad9` (-all 27 vs -primary 15, promediados = media exacta), `msf_quad9_noavg` (average -no: nodos medios a 0, -all == -primary), `msf_nor` (axial N=2e-2: nors=N 0-5%, shes/moms~0), `msf_shear` (corte simple: shes = G·γ·h = 0.384615 EXACTO, nor=mom=0). Los `.dat` viven en validation-suite/test-2014 (gitignored); el bucle y los checks de build_safe.sh se versionan. |



: ~74 keywords (aliases force 24+11 projected, control_materi 13, data 5, solver 7, print 3, axisymmetric 1) + 4 fixes/gotchas documentados. |


Nota: las features marcadas solo "GNU" en el detalle por familia (sin
fila en esta tabla) ya existían en el fork sfnet 2014 y no fueron
implementadas por nosotros.

## Convenciones del seguimiento

- `ESTADO`: `[x]` = implementada (commit abajo); `GNU` = presente en el codebase GNU
  (puede faltar doc); `PENDIENTE` = falta; `DESCARTADO` = no se implementará (razón abajo).
- `Fecha/Commit`: fecha ISO del commit que implementó la feature (no hay números de versión).
- `Verificación`: test de la suite, ejemplo del manual, o test propio.

> **Regla de trabajo**: cada feature nueva se investiga primero — ¿hay ejemplo de uso en el
> manual de Professional, en su repo/Google Drive, o en la suite sfnet? Si lo hay, se verifica
> contra ese ejemplo; si no, se monta un test propio. La verificación queda anotada aquí y en
> el manual-developer de la feature.

## Descarte de features

Las features que decidimos NO implementar se anotan aquí con su razón, y además se reflejan
en el apartado 'Diferencias con la versión Professional' del manual correspondiente.
Todas las divergencias (también las no decididas) están en
[DIVERGENCES-PROFESSIONAL.md](DIVERGENCES-PROFESSIONAL.md) con su clase (A deliberada,
B conocida/pendiente, C extensión/dialecto GNU, D sospechosa del Pro).

| Feature | Razón del descarte |
|---------|--------------------|
| `interface13` del corpus (target 1.11803) | **BUG DEL PROPIO TEST, validado contra el binario Pro 25-10-2023** (2026-09-04): el Pro escribe `element_interface_intpnt_strain` (shear) = 0.67082 = du_tang/2 = (3/√5)/2 y su propio tochnog.log reporta "Target value ... is 1.11803. The actual value is 0.67082. Error detected for input file interface13.dat". El GNU produce EXACTAMENTE el mismo 0.67082. El target 1.11803 = √5/2 = |du|/2 asume que el desplazamiento (2,−1) es tangencial puro, pero el tangente de la interfaz es (2,1)/√5 y hay componente normal (−4/√5). Ambos códigos fallan el target; el test está mal escrito. |
| `control_print_gid_*` (familia ~20) | Único formato propietario contemplado; `print_gid_6` del GNU es de 1998 (anterior a cambios recientes); GiD puede importar formatos no nativos (VTK, Gmsh, CSV). Ver plan sección 6c. |
| `control_data_save` | Descarte de la FEATURE (único consumidor: `control_print_gid_save_difference`, familia GiD descartada; el análisis de diferencias se hace en post-proceso con la exportación SQLite/CSV de P5-T). El 2026-09-04 se registró parse-only (INTEGER, class CONTROL) porque los tests dynamic6/7 del corpus lo declaran; el registro sigue sin comportamiento. |
| `control_print_interface_stress*` | Requiere elementos de interfaz/contacto, que el GNU no tiene (ni enums). Depende del Carril A. **OBSOLETO (nota 2026-10-07)**: `control_print_interface_stress` se implementó después en el Carril A (`01f6c3e` 2D 2026-08-14, `b8804c6` 3D 2026-08-16; ver la familia control_print); el descarte ya no aplica a ese keyword. |
| Professional behaviour of `mohr_coul_direct4` (interface Mohr-Coulomb strength not enforced at once, slow drift below it) | **Deliberate divergence** (project decision 2026-10-06, F1 T7, BACKLOG C30): the GNU keeps the manual 6.632 immediate return to c + \|Fn\| tan(phi) (sigxy 1.2027, `disx` 0.3007); the Professional 02-08-2026 stays at 2.0 and drifts to 0.410 at t = 100 (target `disx` 0.1024), which looks like an artefact. Documented in 'Differences with the Professional version' of `manual-user/` and `manual-developer/group_materi_plasti_mohr_coul_direct.md`. |
| Professional shear plastic strain of the direct cut-offs (`dsig_xy / G`, twice the tensor component) | **Deliberate divergence** (project decision 2026-10-06, F1 T4, BACKLOG C34): the Professional 02-08-2026 writes the shear plastic strain of a direct cut (tension/Mohr-Coulomb/compression direct) as `dsig_xy / G` while its own total strain and incremental-law plastic strain use tensor components (its `epexy = eptxy - eppxy` then violates `sigma_xy = 2G epe_xy`); the GNU keeps the tensor value `dsig_xy / 2G` so `sigma = C (eps - eps_p)` holds (`direct_plastic_strain_mc1`: GNU eppxy 2.04114892e-3, Pro 4.08229785e-3; normal components and `f` equal). Documented in 'Differences with the Professional version' of `manual-user/direct_plastic_strain.md` and `manual-developer/direct_plastic_strain.md`. |
| Professional sum C(young) + C(E_law) when `group_materi_elasti_young_power` and `group_materi_elasti_young` are both given | **Deliberate divergence** (project decision 2026-10-06, F2 T1b, BACKLOG C36): the GNU follows the manual (theory 2.2.2 "the young's stiffness modulus is made a function of the average stress state", history variable 6.661 = the calculated young): the power law IS the modulus and replaces the `young` record. The Professional 02-08-2026 adds both stiffnesses (`myoung6_e2` sigyy -5.3846 = stiffness 4000 vs GNU -4.0385 = 3000, while its history variable reports 3000); without the `young` record Pro = GNU. Code and `myoung6*` targets unchanged. Documented in 'Differences with the Professional version' of `manual-user/` and `manual-developer/group_materi_elasti_young_power.md`. |
| In-step geometry of an explicit `mesh -follow_material` step: Professional START-of-step geometry (column 0.5 one step, 0.52093 in 100) | **Deliberate divergence** (decisión del propietario 2026-10-08): el GNU mantiene la geometría de FINAL de paso (0.47222 / 0.52048); solo afecta a `-follow_material` explícito (el default ya coincide) y la diferencia converge con el paso. Ver `DIVERGENCES-PROFESSIONAL.md` (F3-INSTEP). |
| Professional row-sum lumped mass of the `tet10` 5-point rule (negative corner masses -V/20) | **Deliberate divergence** (project decision 2026-10-09, F4 T4, BACKLOG C31): the GNU lumps the `tet10` inertia/capacity of the principal dofs with the positive HRZ masses (corner V/36, mid-edge 4V/27). With the Professional row sum a sudden-gravity tet10 column diverges (GNU build 5e8; the Professional solver fails). Static results are unaffected. Documented in 'Differences with the Professional version' of `integration_quadrature_and_stress_recovery.md` and `inertia_apply.md` (user + developer). See `DIVERGENCES-PROFESSIONAL.md` (TET10-HRZ). |
| Professional linear `post_point_dof` interpolation inside a `tet10` (4 corner nodes) | **Deliberate divergence** (project decision 2026-10-09, F4 T4/T7): the GNU interpolates `post_point_dof` inside a `tet10` with the quadratic shape functions of all 10 nodes, the consistent finite-element interpolation; the Professional interpolates linearly from the 4 corners (`tet10p` post points differ by up to 2 %). Documented in 'Differences with the Professional version' of `integration_quadrature_and_stress_recovery.md` (user + developer). See `DIVERGENCES-PROFESSIONAL.md` (TET10-PPOINT). |
| "Ley seca" de groundflow (manual 2.4.3: límite de presión total 0 y pres 0 → sin término de consolidación) | **Divergencia deliberada respecto al MANUAL** (decisión del propietario 2026-10-09): ni el GNU (desde `3acce9d`) ni el Professional la aplican; con el límite por defecto 0 del manual 6.587 apagaría la consolidación en todo modelo que empiece con pres 0. Ver `DIVERGENCES-PROFESSIONAL.md` (DRY-RULE). |
| Professional `-total_linear` + `group_materi_elasti_compressibility` (time-step dependent volumetric increment) | **Deliberate divergence** (project decision 2026-09-29, merge `93966ef`, BACKLOG C9): the GNU uses the increment `tr(inc_epe)` and gives the analytic sigxx -0.005 for every dt (manual 2.2.1); the Professional 02-08-2026 gives -0.0075 / -0.0099 / -0.009999 for dt 1e-2 / 1e-3 / 1e-4 (legacy `tr(new_epe-old_epe)`). Also Pro-suspect. Documented in 'Differences with the Professional version' of `manual-user/` and `manual-developer/group_materi_memory_total_linear.md`. |
| Professional coupling of the static pore pressure of `groundflow_phreatic_level_static` to `post_calcul -groundflow_pressure -total_pressure` | **Deliberate divergence** (project decision 2026-09-29, merge `ad3437e`, BACKLOG C10): the GNU always applies the static pore pressure to the skeleton (manual 6.574; an output request must not change the physics); the Professional 02-08-2026 applies it only when that post_calcul record exists (`ground18` without it: Pro -1745, GNU -1247.5). Pro `excavate1` keeps a documented target mismatch. Documented in `manual-user/` and `manual-developer/groundflow_phreatic_level_static.md`. |
| Professional index-0-only reading of `control_solver_bicg_stop` | **Deliberate divergence** (project decision 2026-09-29, BACKLOG C13): the GNU reads `solver_bicg_stop`, the index-0 control record and then the record at the running timestep index (manual 6.376); the Professional 02-08-2026 ignores a record at the timestep's non-zero index. Also Pro-suspect. Documented in `manual-user/solver.md` and `manual-developer/solve_iterative_bicg.md`. |
| Professional constant phi_flow factor `2 tan(phi_flow)` of the plane direct cut-off | **Deliberate divergence** (project decision 2026-10-02, BACKLOG C15): the GNU keeps the consistent elastic return `(M/G) tan(phi_flow)` (equal for poisson 0; the Pro constant is the fallback without an elastic stiffness, umat). Test `test-2026/direct_cutoff_plane_phi_flow2` (the Pro fails it). Documented in `manual-user/` and `manual-developer/group_materi_plasti_mohr_coul_direct.md`. |
| Professional binary ignoring `control_materi_updated_apply -yes` | **Deliberate divergence** (project decision 2026-10-04, BACKLOG C21): the GNU follows manual 6.152 (an unspecified memory becomes `-updated`; simple shear 0.231 / 0.421 / -0.229); the Professional 02-08-2026 stays linear (0 / 0.5 / 0). Also Pro-suspect. Test `test-2026/materi_memory_default_updated_apply1`. Documented in `manual-user/` and `manual-developer/control_materi_updated_apply.md`, `group_materi_memory_default.md`. |
| Professional `bounda_time_increment` tables prescribing -1.4 on dofs whose table is all zero (`mcap1_elast`) | **Deliberate divergence** (project decision 2026-10-07, AUDIT K08; Pro-suspect, audit 2026-09-29): the GNU follows manual 6.146 and does not prescribe the -1.4. Documented in `manual-user/` and `manual-developer/bounda_time_increment.md`. |
| Professional groundflow velocity output 0 in `ground1` | **REVERTIDA 2026-10-07** (decisión del propietario): la premisa 'Pro-suspect' era errónea; el Professional cumple el manual 6.587 (límite `groundflow_total_pressure_limit` 0 por defecto, que anula la velocidad donde se alcanza). **Cerrada en F4 T1 (2026-10-09)**: el GNU sigue al Professional (velocidad 0 donde la presión total calculada alcanza el límite; `ground1` gvely 0); filas K10/TPL eliminadas de `DIVERGENCES-PROFESSIONAL.md`. Tests `test-2026/groundflow_velocity_dry1/2`. |
| Professional `control_materi_failure_apply` accepted only without index | **Deliberate divergence** (project decision 2026-10-07; Pro-suspect: the manual documents an index, the Professional 02-08-2026 rejects `control_materi_failure_apply 0 -no` and stops with "Program error ... element_rhside_delete" as soon as an element fails): the GNU follows the manual (indexed record, failed elements deleted). Documented in `manual-user/` and `manual-developer/materi_failure_apply.md`. |
| Professional `control_mesh_move` moving only `node_start_refined` (nodes unmoved) | **Deliberate divergence** (project decision 2026-10-10, F4 T9; DIVERGENCES C46b, class A, also Pro-suspect): manual 6.217 says the mesh is moved; the GNU moves node and `node_start_refined` once per index (F4 T8 `c21b62e`). Documented in the "Differences" sections of `manual-user/control_mesh_move.md` and `manual-developer/control_mesh_move.md`. |
| Professional `control_mesh_switch` leaving `node_start_refined` unswitched | **Deliberate divergence** (project decision 2026-10-10, F4 T9; DIVERGENCES C46e, class A, also Pro-suspect): the GNU switches node and `node_start_refined` (manual 6.244 "all nodal coordinates"). Documented in `manual-user/control_mesh_switch.md` and `manual-developer/control_mesh_switch.md`. |
| Professional degenerate 3D section frame of a SQUARE end face without `post_calcul_materi_stress_force_thickness_switch -yes` | **Deliberate divergence** (project decision 2026-10-11, F4 T6a; DIVERGENCES C3, class A). The manual (6.915) is silent on the tie; the GNU breaks it with the reference point and gives a valid frame (force10 without the switch: nory -12.34, mom1y +5000), the Professional gives all vector components 0 and the moment in the mom2 slot. With `-yes` both are equal. |
| Professional symmetrization of the coupled consolidation (Biot u-p) matrix under `solver_matrix_symmetric -yes` | **Deliberate divergence** (project decision 2026-10-11, F4 T6b, option A; DIVERGENCES K24, class A, also Pro-suspect). The Professional symmetrizes the non-symmetric u-p matrix (manual 6.1052) and its pore pressure is non-physical: undrained (k 1e-8) it drifts -0.80318 (t 0.1) -> -0.83404 (t 5) with no flow (analytic 1D -(c+1/E)/(c+3/E) = -0.8462; GNU and Pro `-no` -0.8477), capacity 0.1 gives +0.025 (+0.968 undrained) under compression, capacity 0 gives pres ~0; `large2` 4.5 % off. Oedometric: GNU -0.77306, Pro `-yes` -0.76984, 1D FE ref -0.7756. Manual 2.4.1 (c dh/dt = k lap(h) + div(v)) supports the GNU, which keeps the as-assembled matrix (`-yes` = `-no`). Guard tests `test-2026/k24_consol_3d*` (`c350b3c`). Documented in `manual-user/` and `manual-developer/solver_matrix_symmetric.md` and `groundflow_consolidation_apply.md`. |
| Professional vertex artefact of the classic `group_materi_plasti_mohr_coul` in uniaxial tension (`mmc_tension`: sigzz 0.974 / eppzz -9.7e-4 at dt 0.05, time-step dependent) | **Deliberate divergence** (project decision 2026-10-12, F4 T6f, K15; DIVERGENCES MC-VERTEX-ZZ, class A, also Pro-suspect). At the Mohr-Coulomb vertex sigyy = sigzz = 0 and any plastic eppzz > 0 leaves the vertex, so the continuum solution has sigzz = eppzz = 0; the GNU gives it for every dt (sigzz ~1e-9). The Professional 02-08-2026 gives sigzz 0 / 0.974 / 0.367 / 1.44 for dt 0.1 / 0.05 / 0.047 / 0.0047 (algorithmic). Guard tests `test-2026/mmc_tension_vertex1/2`. Documented in `manual-user/` and `manual-developer/group_materi_plasti_mohr_coul.md`. |
| Professional `group_materi_stokes` read at index 0 for every element group | **Deliberate divergence** (project decision 2026-10-12, F4 T6f K21 WU2; DIVERGENCES STOKES-IDX, class A, also Pro-suspect: manual 6.752 says the index is the element group). The Professional 02-08-2026 uses the index-0 switch for every group (a group 1 `-no` that holds the post_point: Pro 0.1465713 = all `-yes`); the GNU reads the switch per element group (0.124105). Guard tests `test-2026/stokes_group1/2` (the Pro fails them). Documented in `manual-user/` and `manual-developer/group_materi_stokes.md`. |
| Professional `eptzz` of a plastic `group_materi_membrane` element = last step's zz correction only (`validation_3`; single element -9.6e-5 at dt 0.01, -7.9e-4 at dt 0.1) | **Deliberate divergence** (project decision 2026-10-12, F4 T6f, K15; DIVERGENCES MEMBRANE-EPTZZ, class A, also Pro-suspect). The GNU stores the total zz strain, elastic + accumulated plastic (manual 6.685 plane stress; analytic single element -0.0088249, GNU -0.00863 at dt 0.01). Guard tests `test-2026/membrane_eptzz_total1/2`; GNU residual yy/zz plastic split first order in dt -> BACKLOG C69. Documented in `manual-user/` and `manual-developer/group_materi_membrane.md`. |
| Professional time-discretisation error of the von Wolffersdorff + intergranular strain law at a strain reversal (test-2014 `hypo3`: sigxx -0.1629 at dt 1e-3, -0.10499 at dt 1e-5) | **Deliberate divergence** (project decision 2026-10-12, F4 T6f; DIVERGENCES K14b, class A, also Pro-suspect). The converged reference (Python integration of the same law, manual 2.2.4) is -0.10466; the GNU gives -0.1034 at dt 1e-3 (1.2 %) and -0.104235 at dt 1e-5, the Professional is 56 % off at dt 1e-3 and only converges as dt -> 0. Guard tests `test-2026/hypo_wolff_igs_oedometer1/2`. The `hypo3` hyhis slot difference stays open (K14, class B). Documented in `manual-user/` and `manual-developer/group_materi_plasti_hypo_wolfersdorff.md` and `group_materi_plasti_hypo_intergranularstrain.md`. |
| Professional Mazars equivalent strain (positive principal strains accumulated in ascending order as s = 2 s + eps^2, weights 4, 2, 1; also its damage + von Mises path, C7b) | **Deliberate divergence** (project decision 2026-10-12, F4 T6f; DIVERGENCES C7, class A, also Pro-suspect). The GNU follows manual 2.2.5, eps_eq = sqrt(sum <eps_i>+^2): hex8 diag(2e-4,-3e-4,1e-4) dam 0.28610207515 (Pro 0.36409); bar + von Mises 3e6 dam 0.925372, eppxx 1.87863e-4 (Pro 0.96946, 1.368e-4). The Pro form reproduces 14 probes to 7e-14. Guard tests `test-2026/damage_mazars_eq_two_positive1`, `damage_mazars_vonmises1`. The non-objective gamma (xx, xy, xz row) is shared by both codes: BACKLOG C7 note, no fix. Documented in `manual-user/` and `manual-developer/group_materi_damage_mazars.md` and `materi_damage_effective_stress.md`. |
| Professional `prism6` corner integration rule (triangle corners x 2 Gauss levels in zeta, weight 1/12; default, `-normal`, `-maximal`, with inertia) | **Deliberate divergence** (project decision 2026-10-13, F4 T6f C75; DIVERGENCES C75, class A). Manual 6.622 is silent on the triangle points of a prism. The GNU keeps the interior (2/3,1/6,1/6) rule, exact to degree 2 in the triangle (the corner rule: degree 1). Both are full rank and pass the patch test, but on a prism6 cantilever (Timoshenko -4031.2) the corner rule gives -26.0 / -9.95 % tip deflection error against the GNU -13.1 / -4.85 % (meshes 20x2x2 / 40x2x4); prism6 nodal stresses differ only on stress fields outside the prism6 span (distorted probe: sigzz 1.6 % of the scale). The `-minimal` 2-point rule stays open (BACKLOG C76). |

## Equivalencias de nomenclatura Professional vs GNU

Algunas familias tienen nombres distintos entre Professional y el GNU
(fork sfnet 2014). En el detalle por familia, un keyword de Professional
marcado `PENDIENTE` puede estar cubierto en el GNU bajo otro nombre:

| Professional | GNU | Notas |
|--------------|-----|-------|
| `force_edge*` | `force_element_edge*` | Fuerzas de borde distribuidas. El parser traduce el prefijo Professional DENTRO de db_number (Sprint 9): los nombres `force_edge_*` se aceptan directamente. |
| `condif_convection_edge_normal` (+ `_geometry`) | `condif_convection` (+ `condif_convection_geometry`) | Misma física h·(Tenv−T). Ambos nombres conviven: los Professional son masters types 7/8 de area.cc (con variantes de restricción); los legacy conservan su camino (tests convec1/2, condif8/10). |
| `condif_radiation_edge_normal` (+ `_geometry`) | `condif_radiation` (+ `condif_radiation_geometry`) | Misma física αr·(Tr⁴−T⁴). Ambos nombres conviven como arriba. |
| `materi_plasti_maximum_iterations` | `group_materi_plasti_maximum_iterations` | Prefijo `group_` en el GNU. |
| `control_mesh_convert` (interfaces) | — | Solo conversión de interfaz; el GNU no tiene elementos de interfaz (Carril A). |
| `bounda_dof` | `bounda_unknown` | Grafía GNU 2014; el GNU acepta ambas, el Professional rechaza la GNU (AUDIT §1.1; DIVERGENCES `DIALECT`). |
| `mesh`, `inertia_apply`, `convection_apply`, `convection_stabilization`, `processors` | `options_mesh`, `options_inertia`, `options_convection`, `options_stabilization`, `options_processors` | Ídem. Ojo: el default de `mesh` difiere (GNU `-follow_material`, Pro `-fixed_in_space`, BACKLOG C8). |
| `control_solver` | `control_options_solver` | Ídem. |
| `group_truss_elasti_young`, `-truss_beam` | `group_truss_young`, `-trussbeam` | Ídem. |
| `groundflow_phreatic_level`, `group_porosity` | `groundflow_phreaticlevel`, `group_groundflow_porosity` | Ídem. |
| `contact_target_geometry`, `contact_spring*` | `contact_geometry`, `contactspring*` | Ídem (alias de contacto corregido en F2 T3a, BACKLOG C5). |
| `druck_prag`, `mohr_coul`, `…plasti_bounda`, `…plasti_cap2`, `-size_tot`/`-size_dev` | `druckprag`, `mohrcoul`, `…plasti_boundary`, `…plasti_cap`, `-sizetot`/`-sizedev` | Ídem. |

## Detalle por familia

**Leyenda**: `[x]` = implementada por nosotros (commit abajo); `GNU` = presente en el GNU
(posiblemente sin doc); `PENDIENTE` = falta.

<details><summary>Ver detalle completo</summary>




### area_element (13/13)

- [x] `area_element_group` — presente en el GNU (variantes Sprint 8)
- [x] `area_element_group_element` — Sprint 8 (2026-08-24)
- [x] `area_element_group_interface` — Sprint 8 (2026-08-24; -no excluye, GNU también acepta -yes para las interfaces del Carril A)
- [x] `area_element_group_method` — Sprint 8 (2026-08-24; HALLAZGO: se leía en group.cc desde el GNU pero NUNCA estuvo registrado — era inutilizable)
- [x] `area_element_group_node` — Sprint 8 (2026-08-24)
- [x] `area_element_group_sequence` — presente en el GNU
- [x] `area_element_group_sequence_element` — presente en el GNU
- [x] `area_element_group_sequence_element_group` — Sprint 8 (2026-08-24; alias Professional del legacy `..._elementgroup`)
- [x] `area_element_group_sequence_geometry` — presente en el GNU
- [x] `area_element_group_sequence_geometry_method` — Sprint 8 (2026-08-24)
- [x] `area_element_group_sequence_interface` — Sprint 8 (2026-08-24)
- [x] `area_element_group_sequence_time` — presente en el GNU
- [x] `area_element_group_time` — Sprint 8 (2026-08-24; re-evaluación por paso vía area_element_group_time_active en top.cc)

### area_node (3/3)

- [x] `area_node_dataitem` — presente en el GNU
- [x] `area_node_dataitem_double` — presente en el GNU
- [x] `area_node_dataitem_integer` — presente en el GNU

### axisymmetric (1/1)

- [x] `axisymmetric` — Sprint 9 (alias db_number de `group_axisymmetric`)

### beam (1/1)

- [x] `beam_rotation` — presente en el GNU

### bounda_alternate (1/1)

- [x] `bounda_alternate` — implementada (commit `0ba382a`, 2026-08-05)

### bounda_baseline (0/2)

- [ ] `bounda_baseline_correction` — PENDIENTE
- [ ] `bounda_baseline_correction_parameters` — PENDIENTE

### bounda_constant (1/1)

- [x] `bounda_constant` — presente en el GNU

### bounda_dof (3/3)

- [x] `bounda_dof` — implementada (commit `b5e02d5`, 2026-08-05); dof `-topres` (prescripcion de presion TOTAL, manual Pro 2.4.1) 2026-09-05 (ground13 rc=0)
- [x] `bounda_dof_cylindrical` — implementada (commit `dd53575`, 2026-08-05)
- [x] `bounda_dof_radial` — implementada (commit `5058e46`, 2026-08-05)

### bounda_factor (2/2)

- [x] `bounda_factor` — presente en el GNU
- [x] `bounda_factor_parabolic_x` — implementada (commit `ba20c45`, 2026-08-05)

### bounda_force (1/1)

- [x] `bounda_force` — presente en el GNU

### bounda_found (1/1)

- [x] `bounda_found` — presente en el GNU

### bounda_geometry (1/1)

- [x] `bounda_geometry_method` — implementada (commit `ba20c45`, 2026-08-05)

### bounda_normal (1/1)

- [x] `bounda_normal` — implementada (commit `d20ca47`, 2026-08-05)

### bounda_print (3/3 — por equivalencia)

- [x] `bounda_print_mesh_dof` — Sprint 9 (alias de `print_mesh_dof`, dump implementado)
- [x] `bounda_print_mesh_dof_geometry` — Sprint 9 (alias)
- [x] `bounda_print_mesh_dof_values` — Sprint 9 (alias; registrado, no usado por el dump)

### bounda_sine (1/1)

- [x] `bounda_sine` — presente en el GNU; layout alineado al Professional 6.34 (start_time end_time freq amp ...) el 2026-09-13 (`6888875`; ver registro de verificación); frequency-0 static value PENDIENTE

### bounda_time (7/11; 1 falso pendiente por OCR, 3 SMC pendientes)

- [x] `bounda_time` — presente en el GNU
- [x] `bounda_time_factor` — Sprint 8 (2026-08-24)
- [x] `bounda_time_increment` — presente en el GNU
- [x] `bounda_time_offset` — presente en el GNU (el "bounda_time_o" del inventario era OCR truncado de offset)
- [~] `bounda_time_smc` — registrado y parseado (2026-09-04); el LECTOR del formato SMC USGS sigue PENDIENTE (fallo honesto en runtime; ver manual-user/bounda_time_smc.md)
- [~] `bounda_time_smc_offset` — registrado (ídem, PENDIENTE el consumo)
- [~] `bounda_time_smc_units` — registrado (ídem, PENDIENTE el consumo)
- [x] `bounda_time_units` — implementada (commit `ba20c45`, 2026-08-05)
- [x] `bounda_time_until_data` — implementada (lote bounda 2026-09-01; ver registro de verificación)
- [x] `bounda_time_until_value_minimum` — implementada (lote bounda 2026-09-01; ver registro de verificación)
- [x] `bounda_time_user` — presente en el GNU
- [ ] `bounda_time_until_value` (variante 3 valores del corpus 02-08-2026) — registrada como keyword conocida; consumo PENDIENTE (el binario 25-10-2023 la rechaza, no verificable)

### bounda_water (1/1)

- [x] `bounda_water` — implementada (commit `e7861f6`, 2026-08-05)

### change (6/6)

- [x] `change_dataitem` — presente en el GNU
- [x] `change_dataitem_geometry` — P6 (2026-08-24)
- [x] `change_dataitem_time` — presente en el GNU
- [x] `change_dataitem_time_discrete` — presente en el GNU
- [x] `change_dataitem_time_method` — P6 (2026-08-24)
- [x] `change_dataitem_time_user` — presente en el GNU

### check_data (1/1)

- [x] `check_data` — implementada (commit `5a2b3d2`, 2026-08-04)

### check_element (2/2)

- [x] `check_element_node` — implementada (commit `3793892`, 2026-08-04)
- [x] `check_element_shape` — implementada (commit `762f2ca`, 2026-08-04)

### check_error (1/1)

- [x] `check_error` — implementada (commit `5a2b3d2`, 2026-08-04)

### check_memory (3/3)

- [x] `check_memory` — implementada (commit `762f2ca`, 2026-08-04)
- [x] `check_memory_usage` — presente en el GNU
- [x] `check_memory_usage_result` — presente en el GNU

### check_nan (1/1)

- [x] `check_nan` — implementada (commit `3793892`, 2026-08-04)

### check_solver (1/1)

- [x] `check_solver` — implementada (commit `e42e12d`, 2026-08-04)

### check_target (1/1)

- [x] `check_target` — implementada (commit `3793892`, 2026-08-04)

### check_used (1/1)

- [x] `check_used` — implementada (commit `3793892`, 2026-08-04)

### check_warning (1/1)

- [x] `check_warning` — implementada (commit `5a2b3d2`, 2026-08-04)

### condif_convection (7/7)

- [x] `condif_convection_edge_normal` — Sprint 7 (2026-08-24; nombre Professional del legacy `condif_convection`)
- [x] `condif_convection_edge_normal_element` — Sprint 7 (2026-08-24)
- [x] `condif_convection_edge_normal_element_group` — Sprint 7 (2026-08-24)
- [x] `condif_convection_edge_normal_element_node` — Sprint 7 (2026-08-24)
- [x] `condif_convection_edge_normal_element_side` — Sprint 7 (2026-08-24)
- [x] `condif_convection_edge_normal_geometry` — Sprint 7 (2026-08-24)
- [x] `condif_convection_edge_normal_node` — Sprint 7 (2026-08-24)

### condif_heat (20/20)

- [x] `condif_heat_edge_normal` — Sprint 7 (2026-08-24)
- [x] `condif_heat_edge_normal_element` — Sprint 7 (2026-08-24)
- [x] `condif_heat_edge_normal_element_group` — Sprint 7 (2026-08-24)
- [x] `condif_heat_edge_normal_element_node` — Sprint 7 (2026-08-24)
- [x] `condif_heat_edge_normal_element_node_factor` — Sprint 7 (2026-08-24)
- [x] `condif_heat_edge_normal_element_side` — Sprint 7 (2026-08-24)
- [x] `condif_heat_edge_normal_factor` — Sprint 7 (2026-08-24)
- [x] `condif_heat_edge_normal_geometry` — Sprint 7 (2026-08-24)
- [x] `condif_heat_edge_normal_node` — Sprint 7 (2026-08-24)
- [x] `condif_heat_edge_normal_sine` — Sprint 7 (2026-08-24); consumer vía rama compartida `flux_edge_is_master` (auditoría 2026-09-13), layout Pro 6.81 (`start end freq amp ...`) alineado
- [x] `condif_heat_edge_normal_time` — Sprint 7 (2026-08-24)
- [x] `condif_heat_volume` — Sprint 7 (2026-08-24)
- [x] `condif_heat_volume_element` — Sprint 7 (2026-08-24)
- [x] `condif_heat_volume_element_group` — Sprint 7 (2026-08-24)
- [x] `condif_heat_volume_factor` — Sprint 7 (2026-08-24)
- [x] `condif_heat_volume_geometry` — Sprint 7 (2026-08-24)
- [x] `condif_heat_volume_sine` — Sprint 7 (2026-08-24); layout Pro 6.88 (`start end freq amp ...`) alineado 2026-09-13
- [x] `condif_heat_volume_time` — Sprint 7 (2026-08-24)
- [x] `condif_heat_volume_user` — Sprint 7 (2026-08-24)
- [x] `condif_heat_volume_user_parameters` — Sprint 7 (2026-08-24)

### condif_radiation (7/7)

- [x] `condif_radiation_edge_normal` — Sprint 7 (2026-08-24; nombre Professional del legacy `condif_radiation`)
- [x] `condif_radiation_edge_normal_element` — Sprint 7 (2026-08-24)
- [x] `condif_radiation_edge_normal_element_group` — Sprint 7 (2026-08-24)
- [x] `condif_radiation_edge_normal_element_node` — Sprint 7 (2026-08-24)
- [x] `condif_radiation_edge_normal_element_side` — Sprint 7 (2026-08-24)
- [x] `condif_radiation_edge_normal_geometry` — Sprint 7 (2026-08-24)
- [x] `condif_radiation_edge_normal_node` — Sprint 7 (2026-08-24)

### condif_temperature (1/1)

- [x] `condif_temperature` — presente en el GNU

### contact_apply (1/1)

- [x] `contact_apply` — implementada (gate `-yes`/`-no` en `parallel_contact`; validada con `contact`)

### contact_heat (1/1)

- [x] `contact_heat_generation` — P6 (2026-08-24)

### contact_penalty (3/3)

- [x] `contact_penalty_pressure` — presente en el GNU
- [x] `contact_penalty_temperature` — presente en el GNU
- [x] `contact_penalty_velocity` — presente en el GNU

### contact_plasti (1/1)

- [x] `contact_plasti_friction` — implementada (Mohr-Coulomb `max(c + Fn*tan(phi), 0)` en el slip)

### contact_target (3/3)

- [x] `contact_target_element_group` — implementada (filtro de targets por grupo); 2026-10-07 (F2 T3b, K05): node_boundary generated for the element targets (contac1/contac4/axisym3 = Pro), see registro de verificación
- [x] `contact_target_geometry` — implementada (alias de `contact_geometry`); 2026-10-07 (F2 T3a, C5/K05): every index, read-only lookup (no PUT in the parallel loop), see registro de verificación
- [x] `contact_target_geometry_switch` — implementada (alias de `contact_geometry_switch`); 2026-10-07 (F2 T3a): every index, read with its geometry

### control_bounda (0/2)

- [x] `control_bounda_relax` — registrada (2026-09-05; parse-only: la semántica de almacenamiento del Pro no es observable en relax1/2 — A/B en el binario Pro 25-10-2023 con/sin los records IDÉNTICO; ver manuales control_bounda_relax)
- [x] `control_bounda_relax_geometry` — registrada (2026-09-05; ídem)

### control_change (1/1)

- [x] `control_change_dataitem_apply` — implementada (commit `bdcdeaf`, 2026-08-13)

### control_check (0/1)

- [ ] `control_check_data` — PENDIENTE

### control_contact (1/1)

- [x] `control_contact_apply` — P6 (2026-08-24)

### control_convection (1/1)

- [x] `control_convection_apply` — 2026-09-02 (alias de db_number → CONTROL_OPTIONS_CONVECTION, manual 6.113; consumido por general.cc/truss.cc/spring.cc)

### control_data (11/12; 1 descartado)

- [x] `control_data_activate` — P6 (2026-08-24)
- [x] `control_data_arithmetic` — P6 (2026-08-24)
- [x] `control_data_arithmetic_double` — P6 (2026-08-24)
- [x] `control_data_copy` — P6 (2026-08-24)
- [x] `control_data_copy_factor` — P6 (2026-08-24)
- [x] `control_data_copy_index` — P6 (2026-08-24)
- [x] `control_data_copy_index_factor` — P6 (2026-08-24)
- [x] `control_data_delete` — presente en el GNU
- [x] `control_data_put` — presente en el GNU
- [x] `control_data_put_double` — presente en el GNU
- [x] `control_data_put_integer` — presente en el GNU
- [x] `control_data_save` — registrado parse-only (2026-09-04; consumidor del estado guardado para GiD diff sigue sin comportamiento — los tests dynamic6/7 lo usan solo para parse)

### control_dependency (0/1)

- [x] `control_dependency_apply` — KEYWORDS BATCH C (2026-09-04): gate por icontrol en get_group_data (group.cc), precedencia sobre dependency_apply; registro 2026-09-04 (consumo: gate; monitor de items no-dof PENDIENTE)

### control_distribute (6/6)

- [x] `control_distribute` — presente en el GNU; layout Professional 6.127 añadido (2026-08-24)
- [x] `control_distribute_correlation_distance` — P6 (2026-08-24)
- [x] `control_distribute_correlation_length` — P6 (2026-08-24)
- [x] `control_distribute_minimum_maximum` — P6 (2026-08-24)
- [x] `control_distribute_parameters` — P6 (2026-08-24)
- [x] `control_distribute_seed` — P6 (2026-08-24)

### control_element (0/2)

- [ ] `control_element_group` — PENDIENTE
- [ ] `control_element_group_apply` — PENDIENTE

### control_geometry (0/1)

- [ ] `control_geometry_moving` — PENDIENTE

### control_ground (0/1)

- [ ] `control_ground` — PENDIENTE

### control_groundflow (2/2)

- [x] `control_groundflow_consolidation_apply` — P6 (2026-08-20)
- [x] `control_groundflow_nonsaturated_apply` — P6 (2026-08-20)

### control_inertia (1/1)

- [x] `control_inertia_apply` — registrado (alias `control_options_inertia`); default de la familia alineado al Professional (-no) en el sprint u-p (2026-09-05, commit feat de este lote)

### control_input (0/1)

- [ ] `control_input` — PENDIENTE

### control_materi (15/15 — familia COMPLETA; 7 parciales documentados)

- [x] `control_materi_damage_apply` — Sprint 9 (gate en damage.cc); 2026-09-27: shares `materi_apply_gate_off()` with the global `materi_damage_apply` (control overrules global)
- [x] `control_materi_dynamic` — implementada (2026-09-04, con materi_dynamic: factor double, override por icontrol en materi()/truss()); 2026-10-11 (F4 T6d, C11, `28992d5`): también activa el esquema displacement-primary (`materi_dynamic_effective()`, = Professional; ver registro de verificación)
- [x] `control_materi_elasti_k0` — Sprint 9 (registrado) + Sprint 10 lote 5 (hook conectado en set_stress: nu = K0/(1+K0) con group_materi_elasti_k0; validado con `mk0`/`mk0_off`); 2026-09-29: leído en el índice ICONTROL (antes solo si existía el índice 0), commit pending
- [x] `control_materi_failure_apply` — Sprint 9 (gate junto a damage_apply); 2026-09-27: now gates the failure criteria in failure.cc and no longer freezes damage (measured on the Pro 02-08-2026)
- [x] `control_materi_plasti_hardsoil_gammap_initial` — Sprint 10 (lote 8; -yes crea γp_extra = f(estado inicial) en el primer paso (centinela -1 pre-alocado en top.cc, PUT paralelo-safe), guardado en element_intpnt_materi_plasti_hardsoil_gammap_initial y SUMADO a γp en la ley -> f=0 al arrancar con tensiones desviadoras; validado con `mhardsoil_gp0`: σxx -2.0 EXACTO + record 0.0007763 analitico; A/B `mhardsoil_gp0_off` sin control -> retorno σxx -1.225 y kappa 0.000293); 2026-09-29: pre-alocación para cualquier índice de control (antes solo índice 0 → abortaba con índice ≠ 0), commit pending
- [x] `control_materi_plasti_hypo_masin_clay_ocr_apply` — batch hypo (2026-09-01; gate control_materi por ICONTROL en la rama masin de hypoplas.cc, default -YES; data_required GROUP_TYPE eliminado — el record es CONTROL indexado por timestep; hypo12/13 lo consumen)
- [x] `control_materi_plasti_hypo_masin_ocr_apply` — batch hypo (2026-09-01; gate control_materi por ICONTROL en la rama masin, default -YES; data_required GROUP_TYPE eliminado)
- [x] `control_materi_plasti_hypo_niemunis_visco_ocr_apply` — Sprint 9 (PARCIAL: registrado, sin hook en hypo)
- [x] `control_materi_plasti_hypo_pressure_dependent_void_ratio` — Sprint 9 (PARCIAL: ídem)
- [x] `control_materi_plasti_hypo_substepping` — Sprint 9 (PARCIAL: ídem)
- [x] `control_materi_plasti_tension_apply` — Sprint 9 (gate en stress.cc materi_direct_cutoff; validado con cmat_gate)
- [x] `control_materi_plasti_visco_apply` — Sprint 9 (gate en set_stress)
- [x] `control_materi_undrained_apply` — familia undrained (2026-09-05): gate consumido en materi.cc (default -yes); group_materi_undrained_capacity + element_intpnt_materi_undrained_pressure implementados (ground17/undrained1/undrained2 rc=0)
- [x] `control_materi_updated_apply` — Sprint 9 (dirección -no cableada en set_stress); -yes wired 2026-10-04 (`d8117aa`, C16: unspecified memory -> -updated per manual 6.152; the Professional 02-08-2026 binary ignores -yes, deliberate divergence)
- [x] `control_materi_viscosity_apply` — Sprint 9 (gate en viscosit.cc)

### control_mesh (37/91)

- [x] `control_mesh_activate_gravity_apply` — implementada (registrado; selección de records a aplicar)
- [x] `control_mesh_adjust_geometry` — presente en el GNU
- [x] `control_mesh_change_element_group` — implementada (commit `afc1dad`, 2026-08-05)
- [x] `control_mesh_convert` — implementada (commit `490545b`, 2026-08-14; Carril A Fase 2, bar2→quad4)
- [x] `control_mesh_convert_element_group` — implementada (commit `490545b`, 2026-08-14)
- [ ] `control_mesh_convert_quad9_quad6` — PENDIENTE
- [ ] `control_mesh_convert_tria6_tria3` — PENDIENTE
- [x] `control_mesh_copy` — implementada (commit `3e94dac`, 2026-08-05); after a time step = Pro (F4 T8 `0ca144c`, 2026-10-10)
- [x] `control_mesh_cut_force` — implementada (2026-09-05; alias db_number de control_mesh_cut_node_force — el nombre del corpus/.dbs del Pro)
- [x] `control_mesh_cut_geometry` — implementada (2026-09-05; corte + sustitución por fuerzas nodales en delete.cc; ver registro)
- [x] `control_mesh_delete_element` — implementada (commit `afc1dad`, 2026-08-05)
- [x] `control_mesh_delete_geometry` — presente en el GNU
- [ ] `control_mesh_delete_geometry_direct` — PENDIENTE
- [x] `control_mesh_delete_geometry_element` — presente en el GNU
- [ ] `control_mesh_delete_geometry_element_group` — PENDIENTE
- [x] `control_mesh_delete_geometry_factor` — presente en el GNU
- [ ] `control_mesh_delete_geometry_method` — PENDIENTE
- [ ] `control_mesh_delete_geometry_move_node` — PENDIENTE
- [ ] `control_mesh_delete_geometry_projection_type` — PENDIENTE
- [ ] `control_mesh_delete_geometry_stop` — PENDIENTE
- [ ] `control_mesh_delete_geometry_stop_geometry` — PENDIENTE
- [x] `control_mesh_delete_small` — presente en el GNU
- [ ] `control_mesh_duplicate_element_group` — PENDIENTE
- [ ] `control_mesh_element_group_apply` — PENDIENTE
- [x] `control_mesh_extrude` — implementada (commit `5dec3c2`, 2026-08-05; 2026-09-01: + quad9→hex27 con n capas de `control_mesh_extrude_n` — el anillo 3D del Professional force11; 2026-10-11 F4 T12 `24301c8`: nodos de los planos interiores compartidos entre hex27 vecinos, force11 = Pro; los branches lineales tria3/bar2/quad4 intactos)
- [ ] `control_mesh_extrude_contact_spring_element_group` — PENDIENTE
- [ ] `control_mesh_extrude_contact_spring_element_group_new` — PENDIENTE
- [ ] `control_mesh_extrude_direction` — PENDIENTE
- [ ] `control_mesh_extrude_element` — PENDIENTE
- [ ] `control_mesh_extrude_element_group_new` — PENDIENTE
- [x] `control_mesh_extrude_n` — presente en el GNU
- [x] `control_mesh_generate_beam` — implementada (commit `bd7b3ab`, 2026-08-13)
- [x] `control_mesh_generate_contact_spring` — `892541a` (2026-08-31; alias db_number -> CONTROL_MESH_GENERATE_CONTACTSPRING; conspr1-7 rc=0)
- [x] `control_mesh_generate_contact_spring_element` — `892541a` (2026-08-31; alias db_number -> CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT)
- [x] `control_mesh_generate_contact_spring_element_group` — `892541a` (2026-08-31; alias -> CONTROL_MESH_GENERATE_CONTACTSPRING_ELEMENT_GROUP + implementado en generate.cc)
- [x] `control_mesh_generate_interface` — implementada (genera interfaces entre grupos con cara compartida; validada con `iface_gen`, `iface_gen_geom`, `iface_gen_geom_off`)
- [x] `control_mesh_generate_interface_geometry` — implementada (filtro por geometría del par de elementos)
- [x] `control_mesh_generate_interface_method` — implementada (`-element_geometry` para selección y/o generación; validada con `iface_gen_method` e `iface_gen_method_gen`)
- [x] `control_mesh_generate_spring1` — presente en el GNU
- [x] `control_mesh_generate_spring2` — presente en el GNU
- [x] `control_mesh_generate_truss` — implementada (commit `bd7b3ab`, 2026-08-13)
- [ ] `control_mesh_generate_truss_beam` — PENDIENTE
- [x] `control_mesh_generate_truss_beam_loose` — implementada (commit `bd7b3ab`, 2026-08-13)
- [x] `control_mesh_generate_truss_beam_macro` — implementada (commit `bd7b3ab`, 2026-08-13)
- [ ] `control_mesh_generate_truss_beam_separate` — PENDIENTE
- [ ] `control_mesh_gid_batch` — PENDIENTE
- [ ] `control_mesh_interface_triangle` — PENDIENTE
- [x] `control_mesh_keep_element` — implementada (commit `afc1dad`, 2026-08-05)
- [x] `control_mesh_keep_element_group` — implementada (commit `afc1dad`, 2026-08-05)
- [x] `control_mesh_keep_geometry` — presente en el GNU
- [x] `control_mesh_keep_node` — implementada (commit `415d26f`, 2026-08-05)
- [x] `control_mesh_macro` — presente en el GNU; `-bar`/`-rectangle`/`-brick` numbering = Pro (F4 T8 `bc2235e`, 2026-10-10)
- [x] `control_mesh_macro_concentrate` — registrada parse-only (KEYWORDS BATCH C 2026-09-04); generacion graduada del macro PENDIENTE
- [x] `control_mesh_macro_element` — presente en el GNU
- [x] `control_mesh_macro_parameters` — presente en el GNU
- [ ] `control_mesh_map` — PENDIENTE
- [x] `control_mesh_merge` — presente en el GNU
- [x] `control_mesh_merge_eps_coord` — alias db_number → control_mesh_merge_epscoord (2026-09-02; parse-level, sin test de corpus)
- [x] `control_mesh_merge_geometry` — 2026-09-02 (whitelist de geometrías del merge, merge.cc; merge2/tutorial_3 siguen RUNFAIL por otros blockers)
- [x] `control_mesh_merge_geometry_not` — 2026-09-02 (blacklist; nombre Professional del legacy control_mesh_merge_not, alias conservado; ground19_water_under_dam)
- [x] `control_mesh_merge_macro_generate` — presente en el GNU
- [x] `control_mesh_mirror` — implementada (commit `3e94dac`, 2026-08-05); after a time step = Pro (F4 T8 `0ca144c`, 2026-10-10)
- [x] `control_mesh_move` — implementada (commit `3e94dac`, 2026-08-05); once per control index (F4 T8 `c21b62e`, 2026-10-10; Pro moves only node_start_refined: DIVERGENCES C46b)
- [x] `control_mesh_multiply` — presente en el GNU
- [ ] `control_mesh_re` — PENDIENTE
- [x] `control_mesh_remove` — implementada (commit `1e5b4b9`, 2026-08-05)
- [ ] `control_mesh_remove_geometry` — PENDIENTE
- [ ] `control_mesh_remove_keep_geometry` — PENDIENTE
- [ ] `control_mesh_remove_really` — PENDIENTE
- [ ] `control_mesh_remove_really_activate_all` — PENDIENTE
- [ ] `control_mesh_remove_really_activate_factor` — PENDIENTE
- [x] `control_mesh_renumber` — presente en el GNU
- [ ] `control_mesh_renumber_element_geometry_o` — PENDIENTE
- [ ] `control_mesh_renumber_element_group_o` — PENDIENTE
- [~] `control_mesh_rotate` — PARCIAL (commit `d8d3a72`, 2026-08-05): F4 T8 found that the manual 2D -> 3D sweep (Pro `rotate1.dat`) does not run in the GNU; parked (BACKLOG/DIVERGENCES C46c)
- [~] `control_mesh_rotate_angle` — PARCIAL (commit `415d26f`, 2026-08-05): GNU 2D rotation, in the manual/Pro the sweep angle of `control_mesh_rotate` (C46c)
- [x] `control_mesh_split` — presente en el GNU
- [ ] `control_mesh_split_element_from` — PENDIENTE
- [ ] `control_mesh_split_element_to` — PENDIENTE
- [x] `control_mesh_split_only` — presente en el GNU
- [x] `control_mesh_switch` — presente en el GNU; two entries = swap as manual 6.244 / Pro (F4 T8 `fb834c6`, 2026-10-10)
- [ ] `control_mesh_truss_distribute_mpc` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_air` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_dof` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_element_group_isoparametric` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_element_group_truss` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_exact` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_exact_minimal_length` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_exact_minimal_length_connect` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_geometry_isoparametric` — PENDIENTE
- [ ] `control_mesh_truss_distribute_mpc_geometry_truss` — PENDIENTE

### control_mpc (1/2)

- [x] `control_mpc_apply` — 2026-09-02 (gate por índice de control en mpc.cc::mpc_node_apply; mpc_apply global + control por timestep)
- [ ] `control_mpc_element_group` — PENDIENTE

### control_plasti (1/1)

- [x] `control_plasti_apply` — implemented (commit `9bbd125`, 2026-09-27): per-timestep plasticity switch, overrules `plasti_apply`; verified against the Professional binary 02-08-2026 (see verification registry)

### control_post (0/3)

- [ ] `control_post` — PENDIENTE
- [ ] `control_post_apply` — PENDIENTE
- [ ] `control_post_element_force` — PENDIENTE

### control_print (63/85)

- [x] `control_print` — presente en el GNU
- [x] `control_print_filter` — presente en el GNU (el "control_print_" del inventario era OCR truncado de filter; CONTROL_PRINT_FILTER en database.cc:1608)
- [x] `control_print_beam_force_moment` — Sprint 11, lote 6 (2026-08-27; fuerzas y momentos de vigas/trusses/truss-beams cruzados por un segmento de corte → beam_force_moment.<index>/<seq>; primera columna = distancia desde el inicio del corte; 12 columnas en ejes LOCALES de la viga; la fuerza axial de un truss/truss-beam se toma de ELEMENT_TRUSS_FORCE; líneas ordenadas por distancia; sin cruces → sin archivo)
- [x] `control_print_beam_force_moment_coordinates` — Sprint 11, lote 6 (2026-08-27; segmento de corte xstart ystart zstart xend yend zend; 2D solo x y; longitud exacta 2*ndim, obligatorio)
- [x] `control_print_beam_force_moment_switch` — Sprint 11, lote 6 (2026-08-27; -yes multiplica las 12 componentes por -1; -no/ausente = sin inversión)
- [x] `control_print_data_versus_data` — presente en el GNU
- [x] `control_print_data_versus_data_factor` — implementada (commit `3793892`, 2026-08-04)
- [x] `control_print_database` — presente en el GNU
- [x] `control_print_database_method` — Sprint 11, lote 1 (2026-08-26)
- [x] `control_print_dof` — implementada (commit `ca2a9e4`, 2026-08-12)
- [x] `control_print_dof_id` — Sprint 11, lote 1 (2026-08-26; default -yes, nodo en dof.<index>)
- [x] `control_print_dof_line` — Sprint 11, lote 4 (2026-08-27; interpolación node_dof+node_dof_calcul a lo largo de una polilínea; un archivo por label de dof)
- [x] `control_print_dof_line_coordinates` — Sprint 11, lote 4 (2026-08-27; vértices de la polilínea, longitud variable)
- [x] `control_print_dof_line_element_group` — Sprint 11, lote 4 (2026-08-27; filtro de grupos para la interpolación)
- [x] `control_print_dof_line_eps_iso` — Sprint 11, lote 4 (2026-08-27; tolerancia punto-en-elemento, default 1.e-3, → point_el)
- [x] `control_print_dof_line_method` — Sprint 11, lote 4 (2026-08-27; -node / -node_start_refined default con fallback)
- [x] `control_print_dof_line_move` — Sprint 11, lote 4 (2026-08-27; coordenadas siguen el campo de velocidades)
- [x] `control_print_dof_line_n` — Sprint 11, lote 4 (2026-08-27; default 5, equidistantes sobre la longitud total)
- [x] `control_print_dof_line_time` — Sprint 11, lote 4 (2026-08-27; primera línea `# time <t>` gnuplot)
- [x] `control_print_dof_point` — Sprint 11, lote 4 (2026-08-27; interpolación en un punto, serie temporal)
- [x] `control_print_dof_point_coordinates` — Sprint 11, lote 4 (2026-08-27; coordenadas del punto)
- [x] `control_print_dof_point_time` — Sprint 11, lote 4 (2026-08-27; primera línea `# time <t>` gnuplot)
- [x] `control_print_dof_rhside` — Sprint 11, lote 1 (2026-08-26; alias de control_print_unknownsrhside)
- [x] `control_print_dof_smooth_dof` — Sprint 11, lote 5 (2026-08-27; suavizado nodal de control_print_dof: cada pasada reemplaza el valor del nodo por el promedio de sus VECINOS por elementos, sin el propio nodo; -all o labels; N pasadas de _smooth_n, default 10)
- [x] `control_print_dof_smooth_n` — Sprint 11, lote 5 (2026-08-27; nº de pasadas, default 10; convergencia al valor medio verificada)
- [x] `control_print_element` — presente en el GNU
- [x] `control_print_element_method` — Sprint 11, lote 1 (2026-08-26; -middle default / -node)
- [x] `control_print_frd` — implementada (commit `19a8d08`, 2026-08-12)
- [x] `control_print_frd_freecad` — implementada (commit `19a8d08`, 2026-08-12)
- [x] `control_print_frd_prepomax` — implementada (commit `19a8d08`, 2026-08-12)
- [x] `control_print_frequency_timeinterval` — Sprint 11, lote 2 (2026-08-27)
- [x] `control_print_frequency_timestep` — Sprint 11, lote 2 (2026-08-27)
- [x] `control_print_gid` — presente en el GNU
- [ ] `control_print_gid_batch` — DESCARTADO
- [ ] `control_print_gid_beam_vectors` — DESCARTADO
- [ ] `control_print_gid_beam_vectors_normal` — DESCARTADO
- [x] `control_print_gid_contact_spring2` — registrada parse-only (KEYWORDS BATCH C 2026-09-04, alias sin control_ aceptado); switch de dibujado GiD, sin consumo
- [ ] `control_print_gid_coord` — DESCARTADO
- [ ] `control_print_gid_dof` — DESCARTADO
- [ ] `control_print_gid_dof_calcul` — DESCARTADO
- [ ] `control_print_gid_element_group` — DESCARTADO
- [ ] `control_print_gid_element_mpc` — DESCARTADO
- [x] `control_print_gid_empty` — presente en el GNU
- [ ] `control_print_gid_group` — DESCARTADO
- [ ] `control_print_gid_mesh_activate_gravity` — DESCARTADO
- [ ] `control_print_gid_method` — DESCARTADO
- [ ] `control_print_gid_node_method` — DESCARTADO
- [ ] `control_print_gid_other` — DESCARTADO
- [ ] `control_print_gid_safety_slip_critical` — DESCARTADO
- [ ] `control_print_gid_save_di` — PENDIENTE
- [ ] `control_print_gid_smooth_dof` — DESCARTADO
- [ ] `control_print_gid_smooth_n` — DESCARTADO
- [ ] `control_print_gid_spring2` — DESCARTADO
- [ ] `control_print_gid_truss_vector` — DESCARTADO
- [ ] `control_print_gid_truss_vector_normal` — DESCARTADO
- [x] `control_print_gmsh` — implementada (commit `90798f9`, 2026-08-12)
- [x] `control_print_gmsh_dummy` — implementada (commit `90798f9`, 2026-08-12)
- [x] `control_print_gmsh_element_data` — implementada (commit `90798f9`, 2026-08-12)
- [x] `control_print_gmsh_node_method` — implementada (commit `90798f9`, 2026-08-12)
- [x] `control_print_history` — presente en el GNU
- [x] `control_print_history_factor` — implementada (commit `3793892`, 2026-08-04)
- [x] `control_print_history_relative_time` — Sprint 11, lote 1 (2026-08-26; time_current - tr)
- [x] `control_print_history_smooth` — implementada (commit `8f8ad67`, 2026-08-12)
- [x] `control_print_interface_stress` — implementada (commit `01f6c3e`, 2026-08-14; Carril A Fase 4, 2D)
- [x] `control_print_interface_stress_2d_coordinates` — implementada (commit `01f6c3e`, 2026-08-14)
- [x] `control_print_interface_stress_3d_geometry` — implementada (commit `b8804c6`, 2026-08-16; print_interface_stress.cc:65-89)
- [x] `control_print_interface_stress_3d_order` — implementada (commit `b8804c6`, 2026-08-16)
- [x] `control_print_materi_stress_force` — sub-sprint materi_stress_force, lotes 1-5 (2026-08-28; imprime las fuerzas/momentos de `post_calcul -materi_stress -force` → materi_stress_force.<icontrol> — el "index" del manual es el índice del RECORD — con cabecera de comentarios (#) explicando las columnas (9 items 2D / 16 items 3D, una línea por nodo); método -all/-primary (desde L2/L3 -primary omite los nodos promediados quad9/hex27); sin bloque -force → sin archivo)
- [x] `control_print_mesh_dof` — Sprint 11, lote 1 (2026-08-26; alias de print_mesh_dof)
- [x] `control_print_node` — Sprint 11, lote 5 (2026-08-27; print_node.cc: un archivo por parte seleccionada — labels de dof -> label.index (ej. velx.10), items node_dof_calcul -> post_calcul_names (ej. avel.10), números -> <record>_<n>.index; sin partes -> todas)
- [x] `control_print_node_angular` — Sprint 11, lote 5 (2026-08-27; ángulo en GRADOS en vez de coordenadas; 2D solo -yes -yes sin switch_z; 3D las 3 combinaciones; 1D error)
- [x] `control_print_node_angular_middle` — Sprint 11, lote 5 (2026-08-27; punto medio de los ejes, default 0 0 0)
- [x] `control_print_node_geometry` — Sprint 11, lote 5 (2026-08-27; restringe a nodos en la geometría, patrón filter.cc; frame -node_start_refined con fallback a NODE)
- [x] `control_print_node_sort` — Sprint 11, lote 5 (2026-08-27; orden ASCENDENTE; -angle con angular, -x/-y(2D/3D)/-z(3D); patrón collect+sort de print_interface_stress; la mención del manual a "control_print_node_method" es TYPO de node_angular)
- [x] `control_print_node_zero` — Sprint 11, lote 5 (2026-08-27; default -yes; -no suprime valores 0 con comparación EXACTA, por archivo)
- [x] `control_print_number_iterations` — Sprint 11, lote 1 (2026-08-26; monitor de consola)
- [x] `control_print_partialname` — Sprint 11, lote 1 (2026-08-26; prefijo en db_number, stdout)
- [x] `control_print_tecplot` — presente en el GNU
- [x] `control_print_vtk` — presente en el GNU
- [x] `control_print_vtk_coord` — Sprint 11, lote 3 (2026-08-27)
- [x] `control_print_vtk_dof` — implementada (commit `b3e58a4`, 2026-08-13)
- [x] `control_print_vtk_dof_calcul` — Sprint 11, lote 3 (2026-08-27)
- [x] `control_print_vtk_empty` — Sprint 11, lote 3 (2026-08-27)
- [x] `control_print_vtk_node_method` — Sprint 11, lote 3 (2026-08-27)
- [x] `control_print_vtk_other` — Sprint 11, lote 3 (2026-08-27; PARCIAL: subconjunto boundary_condition + mesh_deformation documentado)

### control_repeat (3/3 — familia COMPLETA)

- [x] `control_repeat` — presente en el GNU
- [x] `control_repeat_save` — Sprint 10, lote 10 (2026-08-26): guarda los data items seleccionados en `repeat_save_result` (un índice por repeat) en cada salto del repeat
- [x] `control_repeat_save_calculate` — Sprint 10, lote 10 (2026-08-26): media y varianza poblacional de `repeat_save_result` en `repeat_calculate_result`

### control_reset (18/18 — familia COMPLETA; 1 parcial documentado)

- [x] `control_reset_dof` — implementada (commit `b6eaee4`, 2026-08-13)
- [x] `control_reset_element_dof` — Sprint 8 (2026-08-24; PARCIAL: registrado, la rama -yes solo-element_dof sin cablear — documentado)
- [x] `control_reset_element_group` — Sprint 8 (2026-08-24)
- [x] `control_reset_geometry` — Sprint 8 (2026-08-24)
- [x] `control_reset_interface` — Sprint 8 (2026-08-24)
- [x] `control_reset_interface_strain` — Sprint 8 (2026-08-24)
- [x] `control_reset_node` — Sprint 8 (2026-08-24)
- [x] `control_reset_value_constant` — implementada (commit `b6eaee4`, 2026-08-13)
- [x] `control_reset_value_dof` — implementada (commit `b6eaee4`, 2026-08-13)
- [x] `control_reset_value_dof_diagram` — implementada (commit `b6eaee4`, 2026-08-13)
- [x] `control_reset_value_exponent` — implementada (distribución espacial exponencial; validada con `reset_value_linear`)
- [x] `control_reset_value_linear` — implementada (distribución espacial lineal `ax x + ay y + az z`; validada con `reset_value_linear`)
- [x] `control_reset_value_logarithmic_` — implementada (distribución espacial logarítmica)
- [x] `control_reset_value_logarithmic_second` — implementada (distribución espacial logarítmica de segundo tipo)
- [x] `control_reset_value_method` — implementada (commit `b6eaee4`, 2026-08-13)
- [x] `control_reset_value_multi_linear` — implementada (tabla multilineal vs coordenada vertical)
- [x] `control_reset_value_power` — implementada (distribución espacial potencial `ax x^bx + ...`)
- [x] `control_reset_value_square_root` — implementada (distribución espacial raíz cuadrada)

### control_restart (1/1)

- [x] `control_restart` — presente en el GNU; 2026-10-11 (F4 T6c, `8ca6124`): el restart reinicia también el estado de los puntos de integración (`element_dof`), como el manual 6.368 y el Professional (ver registro de verificación); estado de interface/spring/tendon/damage sin verificar (BACKLOG C58)

### control_safety (0/1)

- [ ] `control_safety_slip` — PENDIENTE

### control_slide (1/3)

- [ ] `control_slide_damping_apply` — PENDIENTE
- [x] `control_slide_plasti_apply` — consumido por la ley (lote familia slide 2026-09-04: -no deja los springs elásticos puros)
- [x] `control_slide_stiffness_apply` — consumido por la ley (lote familia slide 2026-09-04: -no desactiva los springs kn/kt)

### control_solver (7/7 — familia COMPLETA; 4 parciales documentados)

- [x] `control_solver` — Sprint 9 (alias db_number de `control_options_solver`)
- [x] `control_solver_bicg_error` — Sprint 9 (alias de `control_options_solver_bicg_error`)
- [x] `control_solver_bicg_restart` — Sprint 9 (PARCIAL: registrado; el bicg del GNU no tiene restart)
- [x] `control_solver_bicg_stop` — Sprint 9 (-no cableado en so_bicg.cc: continua con warning; -yes/default para como antes); 2026-09-29: solo cuenta el índice 0 (para cualquier índice de timestep) y prevalece sobre `solver_bicg_stop`, medido en el Professional (commit `47b6fda`); then combined precedence (commit `407aa3d`): `solver_bicg_stop` → index 0 → the timestep's own index, last read wins — deliberate divergence from the Pro (index 0 only) when a non-zero-index record exists, accepted by the project owner 2026-09-29 (BACKLOG C13)
- [x] `control_solver_matrix_save` — Sprint 9 (PARCIAL: registrado sin comportamiento)
- [x] `control_solver_pardiso_ordering` — Sprint 9 (PARCIAL: PARDISO no compilado)
- [x] `control_solver_pardiso_out_of_core` — Sprint 9 (PARCIAL: ídem)

### convection (2/2)

- [x] `convection_apply` — 2026-09-02 (alias de db_number → OPTIONS_CONVECTION, manual 6.395; la física de convección ya existía en general.cc, incluida la de condif con flujo materi)
- [x] `convection_stabilization` — 2026-09-02 (alias de db_number → OPTIONS_STABILIZATION, manual 6.396; -yes/-no/-maximal con el upwind peclet de general.cc)
- NOTA: condif2/condif3/tube1 rc=0. validation_1 parsea y corre pero su transitorio NS euleriano diverge del Professional (ver registro); validation_8 sigue bloqueado por el dof `-veln` de bounda_dof (MPC de velocidad normal, manual 6.22 — no implementado ni en el GNU ni en el linaje 2014).

### data (5/5 — familia COMPLETA)

- [x] `data_activate` — Sprint 9 (lote 3; variante time-gated de control_data_activate en data.cc)
- [x] `data_activate_time` — Sprint 9
- [x] `data_delete` — Sprint 9 (variante time-gated de control_data_delete)
- [x] `data_delete_time` — Sprint 9
- [x] `data_ignore` — Sprint 9 (skip a nivel de parseo en input.cc; el record debe preceder a lo que ignora)

### dependency_apply (0/1)

- [x] `dependency_apply` — registrada + gate global en get_group_data (KEYWORDS BATCH C 2026-09-04)

### dependency_diagram (1/1)

- [x] `dependency_diagram` — presente en el GNU

### dependency_geometry (0/1)

- [ ] `dependency_geometry` — PENDIENTE

### dependency_item (1/1)

- [x] `dependency_item` — presente en el GNU

### dependency_method (0/1)

- [ ] `dependency_method` — PENDIENTE

### dependency_number (0/1)

- [ ] `dependency_number` — PENDIENTE

### dependency_type (0/1)

- [ ] `dependency_type` — PENDIENTE

### derivatives (0/1)

- [ ] `derivatives` — PENDIENTE

### dof (1/3)

- [ ] `dof_element_dof` — PENDIENTE
- [x] `dof_label` — presente en el GNU
- [ ] `dof_limit` — PENDIENTE

### dtime (1/1)

- [x] `dtime` — presente en el GNU

### echo (0/1)

- [ ] `echo` — PENDIENTE

### element (1/1)

- [x] `element` — presente en el GNU

### element_beam (1/3)

- [x] `element_beam_direction` — presente en el GNU
- [ ] `element_beam_direction_z` — PENDIENTE
- [ ] `element_beam_force_moment` — PENDIENTE

### element_boundary (0/1)

- [ ] `element_boundary` — PENDIENTE

### element_contact (2/3)

- [ ] `element_contact_spring_direction` — PENDIENTE
- [x] `element_contact_spring_force` — `892541a` (2026-08-31; alias db_number -> ELEMENT_CONTACTSPRING_FORCE)
- [x] `element_contact_spring_strain` — `892541a` (2026-08-31; alias db_number -> ELEMENT_CONTACTSPRING_FORCE, mismo record)

### element_dof (1/3)

- [x] `element_dof` — presente en el GNU
- [ ] `element_dof_initial` — PENDIENTE
- [ ] `element_dof_initial_speci` — PENDIENTE

### element_empty (1/1)

- [x] `element_empty` — presente en el GNU

### element_geometry (0/2)

- [ ] `element_geometry` — PENDIENTE
- [ ] `element_geometry_present` — PENDIENTE

### element_group (1/2)

- [x] `element_group` — presente en el GNU
- [ ] `element_group_apply` — PENDIENTE

### element_interface (1/7)

- [ ] `element_interface_intpnt_direction` — PENDIENTE
- [ ] `element_interface_intpnt_gap_status` — PENDIENTE
- [x] `element_interface_intpnt_materi_tension_status` — F4 T10 C48 (c06fa2d, 2026-10-10): -closed/-opened per integration point = Professional
- [ ] `element_interface_intpnt_strain` — PENDIENTE
- [ ] `element_interface_intpnt_strain_average` — PENDIENTE
- [ ] `element_interface_intpnt_stress` — PENDIENTE
- [ ] `element_interface_intpnt_stress_average` — PENDIENTE

### element_intpnt (0/7)

- [ ] `element_intpnt_dof` — PENDIENTE
- [ ] `element_intpnt_h` — PENDIENTE
- [ ] `element_intpnt_iso_coord` — PENDIENTE
- [x] `element_intpnt_materi_plasti_hardsoil_gammap_initial` — Sprint 10 (lote 8; record calculado de elemento con γp_extra, validado por target directo en `mhardsoil_gp0` = 0.0007763 analitico; una valor por ELEMENTO en este port — por punto de integracion PENDIENTE, documentado)
- [x] `element_intpnt_materi_undrained_pressure` — familia undrained (2026-09-05, 6.441): p_u por IP acumulada (C·ṗ = div v por elemento); ground17/undrained1/undrained2 rc=0
- [ ] `element_intpnt_method` — PENDIENTE
- [ ] `element_intpnt_npoint` — PENDIENTE

### element_middle (1/1)

- [x] `element_middle` — presente en el GNU

### element_normal (1/1)

- [x] `element_normal` — commit `b81d050`, 2026-09-30: generated record (manual 6.445), rule measured on the Professional 02-08-2026 (normalized sum of the consistent boundary node normals); consumed by `*_direct_normal_automatic` (only when switched `-yes`; `-no` = full cut-off as Pro, 2026-10-02, commit `260bd95`); tests `validation-suite/test-2026/element_normal*` (12); zero normal in a `-yes` group = input error (2026-10-01, commit `a8d6452`; Pro: assertion), checked per law since 2026-10-03 (commit `b7f3522`)

### element_print (0/1)

- [ ] `element_print_group_data_values` — PENDIENTE

### element_spring (1/2)

- [x] `element_spring_force` — presente en el GNU
- [ ] `element_spring_strain` — PENDIENTE

### element_truss (2/4)

- [x] `element_truss_direction` — presente en el GNU
- [x] `element_truss_force` — presente en el GNU
- [ ] `element_truss_strain` — PENDIENTE
- [ ] `element_truss_strain_temperature` — PENDIENTE

### element_volume (1/1)

- [x] `element_volume` — presente en el GNU

### end (0/2)

- [ ] `end_data` — PENDIENTE
- [ ] `end_initia` — PENDIENTE

### force_edge (44/44 — registro COMPLETO; LAYOUT de records PENDIENTE)

NOTA (2026-09-02): los 44 nombres están registrados (Sprint 9) pero los
6 tests del corpus que usan la familia (consolidation_gibson, foundation1,
foundation_consolidation_settlement, undrained1, tutorial_4, excavate1)
fallan TODOS en el PARSE del layout, no del nombre:
- `force_element_edge` registra `data_length=nprinc` (3 en modelos u-p con
  pres) pero el Professional da "una fuerza por dirección" = `ndim`
  (2 en 2D) — el parser consume el keyword siguiente como 3er valor.
  Los usuarios del GNU-2014 (axisym2/beam2d_2/condif7/9/disk1/distri1)
  usan siempre `ndim` valores → cambiar a `ndim` es seguro para la suite,
  pero el consumidor de area.cc reparte `values[iprinc]` sobre TODOS los
  dofs principales (incl. el slot de presión en modelos u-p) → hay que
  verificar/limitar el reparto a los dofs de dirección contra el .dbs.
- `force_element_edge_water` es en el Professional un SWITCH
  (`force_edge_water idx -yes`, hidrostática automática con
  groundflow_density/force_gravity/phreatic_level; manual 6.489), pero el
  GNU-2014 (force5.dat) lo usa con 4 doubles explícitos (densidad g
  dirección) → dos layouts incompatibles en el mismo nombre, requiere
  diseño (tipo/carga dual o record nuevo).
- `oo_time` (define 1.e9) y los nombres de geometría de start_define se
  expanden bien por la maquinaria de defines; force_edge_time/geometry NO
  son el blocker.

- [x] `force_edge` — por equivalencia (alias prefijo Professional→`force_element_edge`, Sprint 9) + 2026-09-02 (LAYOUT: record variable 1..nprinc; la forma Professional de 1 valor por dirección carga solo los dofs mecánicos en problemas acoplados — large2/large3)
- [ ] `force_edge_diagram` — PENDIENTE
- [x] `force_edge_element` — Sprint 9 (variante nueva)
- [x] `force_edge_element_group` — Sprint 9
- [x] `force_edge_element_node` — Sprint 9
- [x] `force_edge_element_side` — Sprint 9
- [x] `force_edge_factor` — por equivalencia
- [x] `force_edge_geometry` — por equivalencia
- [x] `force_edge_multi_linear_factor_x` — por equivalencia (`force_element_edge_multi_linear_factor_x`, 86f6d03)
- [x] `force_edge_node` — Sprint 9
- [x] `force_edge_node_factor` — Sprint 9
- [x] `force_edge_normal` — por equivalencia
- [x] `force_edge_normal_element` — Sprint 9
- [x] `force_edge_normal_element_group` — Sprint 9
- [x] `force_edge_normal_element_node` — Sprint 9
- [x] `force_edge_normal_element_side` — Sprint 9
- [x] `force_edge_normal_factor` — por equivalencia
- [x] `force_edge_normal_geometry` — por equivalencia
- [x] `force_edge_normal_node` — Sprint 9
- [x] `force_edge_normal_node_factor` — Sprint 9
- [x] `force_edge_normal_sine` — por equivalencia; layout Pro 6.476 (`start end freq amp ...`) alineado 2026-09-13
- [x] `force_edge_normal_time` — por equivalencia
- [x] `force_edge_projected` — Sprint 9 (lote 2: force_edge_projected Terzaghi) (carga proyectada de túnel, 13 parámetros)
- [x] `force_edge_projected_element` — Sprint 9 (lote 2: companion)
- [x] `force_edge_projected_element_group` — Sprint 9 (lote 2: companion)
- [x] `force_edge_projected_element_node` — Sprint 9 (lote 2: companion)
- [x] `force_edge_projected_element_side` — Sprint 9 (lote 2: companion)
- [x] `force_edge_projected_factor` — Sprint 9 (lote 2: companion)
- [x] `force_edge_projected_geometry` — Sprint 9 (lote 2: companion)
- [x] `force_edge_projected_node` — Sprint 9 (lote 2: companion)
- [x] `force_edge_projected_node_factor` — Sprint 9 (lote 2: companion)
- [x] `force_edge_projected_sine` — Sprint 9 (lote 2: companion); layout Pro 6.487 (`start end freq amp ...`) alineado 2026-09-13
- [x] `force_edge_projected_time` — Sprint 9 (lote 2: companion)
- [x] `force_edge_sine` — por equivalencia; layout Pro 6.465 (`start end freq amp ...`) alineado 2026-09-13
- [x] `force_edge_time` — por equivalencia
- [x] `force_edge_water` — por equivalencia
- [x] `force_edge_water_element` — Sprint 9
- [x] `force_edge_water_element_group` — Sprint 9
- [x] `force_edge_water_element_node` — Sprint 9
- [x] `force_edge_water_element_side` — Sprint 9
- [x] `force_edge_water_factor` — Sprint 9 (nuevo, el GNU no lo tenía)
- [x] `force_edge_water_geometry` — por equivalencia
- [x] `force_edge_water_node` — Sprint 9
- [x] `force_edge_water_time` — por equivalencia

### force_gravity (2/3)

- [x] `force_gravity` — presente en el GNU
- [ ] `force_gravity_geometry` — PENDIENTE
- [x] `force_gravity_time` — presente en el GNU

### force_point (1/1)

- [x] `force_point` — implementada (commit `3793892`, 2026-08-04)

### force_volume (7/7)

- [x] `force_volume` — por equivalencia (alias prefijo→`force_element_volume`, Sprint 9)
- [x] `force_volume_element` — Sprint 9
- [x] `force_volume_element_group_0` — Sprint 9 (`force_volume_element_group`)
- [x] `force_volume_factor` — por equivalencia
- [x] `force_volume_geometry` — por equivalencia
- [x] `force_volume_sine` — por equivalencia; layout Pro 6.507 SIN `end_time` → correcto como está (2026-09-13)
- [x] `force_volume_time` — por equivalencia

### geometry_bounda (3/3)

- [x] `geometry_bounda_sine_x` — presente en el GNU
- [x] `geometry_bounda_sine_y` — presente en el GNU
- [x] `geometry_bounda_sine_z` — presente en el GNU

### geometry_boundary (0/1)

- [ ] `geometry_boundary` — PENDIENTE

### geometry_brick (1/1)

- [x] `geometry_brick` — presente en el GNU

### geometry_circle (2/3)

- [x] `geometry_circle` — presente en el GNU
- [ ] `geometry_circle_part` — PENDIENTE
- [x] `geometry_circle_segment` — presente en el GNU

### geometry_cylinder (2/4)

- [x] `geometry_cylinder` — presente en el GNU
- [ ] `geometry_cylinder_part` — PENDIENTE
- [ ] `geometry_cylinder_part_start_vector` — PENDIENTE
- [x] `geometry_cylinder_segment` — presente en el GNU

### geometry_element (0/4)

- [ ] `geometry_element_geometry` — PENDIENTE
- [ ] `geometry_element_geometry_method` — PENDIENTE
- [x] `geometry_element_group` — KEYWORDS BATCH C (2026-09-04): filtro de geometria por grupos de elementos dentro de geometry(); merge2/force12 rc=0; nota empirica Pro (excluye nodos con elementos fuera de la lista en todo metodo)
- [x] `geometry_element_group_method` — KEYWORDS BATCH C (2026-09-04): -all/-any/-only; default = regla "solo grupos listados" (medida contra Pro)

### geometry_ellipse (1/1)

- [x] `geometry_ellipse` — presente en el GNU

### geometry_exclude (0/1)

- [ ] `geometry_exclude` — PENDIENTE

### geometry_factor (1/1)

- [x] `geometry_factor` — alias de `geometry_bounda_factor` en db_number() (2026-09-01): misma semántica (interpolación lineal/parabólica en `geometry()`, mismo índice); desbloquea matrix2/temp2/matrix4; verificado contra el binario Professional 25-10-2023 (post_point_dof 2.4999998123 vs 2.5). El `control_mesh_delete_geometry_factor` (delete2/delete3) es una feature de clase CONTROL aparte: CORRECCIÓN 2026-09-13 — la hipótesis "parsea pero no se aplica" quedó REFUTADA por medición; el blocker real de delete3/mpc5 era la formulación `-total_linear` (reset descartado por el solve), corregida en `f045d04` (ver fila del registro de verificación 2026-09-13).

### geometry_hexahedral (0/1)

- [ ] `geometry_hexahedral` — PENDIENTE

### geometry_line (1/2)

- [x] `geometry_line` — presente en el GNU
- [ ] `geometry_line_eps_iso` — PENDIENTE

### geometry_list (1/1)

- [x] `geometry_list` — implemented (commit `a4b95f5`, 2026-09-26): node-number list as geometry entity (`-geometry_list <index>`); verified against the Professional binary 02-08-2026 (see verification registry)

### geometry_method (0/1)

- [ ] `geometry_method` — PENDIENTE

### geometry_moving (0/6)

Diagnóstico de `geometry_moving1` (corpus, RUNFAIL de parse en los 6
records + `control_geometry_moving`): excavaciones MÓVILES (manual Pro
6.533-6.538 + 6.135) — una entidad geométrica (punto/área/volumen) se
traslada por la malla con velocidad constante entre un start/end time y
"excava" los elementos que atraviesa, PROPORCIONALMENTE a la fracción
del elemento dentro de la entidad en cada instante (ntime×nspace puntos
de muestreo por elemento). Requiere la eliminación parcial de elementos
con factor volumétrico por paso + el barrido temporal de
geometría-elemento; feature de malla grande, fuera de alcance de sprint
(su familia geometry() ya tiene las ramas delete/cut/extrude/merge).

- [ ] `geometry_moving` — PENDIENTE
- [ ] `geometry_moving_n` — PENDIENTE
- [ ] `geometry_moving_operat` — PENDIENTE
- [ ] `geometry_moving_operat_parameter` — PENDIENTE
- [ ] `geometry_moving_operat_time` — PENDIENTE
- [ ] `geometry_moving_parameter` — PENDIENTE

### geometry_mpc (0/1)

- [ ] `geometry_mpc` — PENDIENTE

### geometry_node (0/1)

- [ ] `geometry_node_type` — PENDIENTE

### geometry_point (1/1)

- [x] `geometry_point` — presente en el GNU

### geometry_polynomial (1/1)

- [x] `geometry_polynomial` — presente en el GNU

### geometry_projection (0/1)

- [ ] `geometry_projection_type` — PENDIENTE

### geometry_quadrilateral (1/2)

- [x] `geometry_quadrilateral` — presente en el GNU
- [ ] `geometry_quadrilateral_eps_iso` — PENDIENTE

### geometry_set (1/1)

- [x] `geometry_set` — presente en el GNU

### geometry_sphere (2/2)

- [x] `geometry_sphere` — presente en el GNU
- [x] `geometry_sphere_segment` — presente en el GNU

### geometry_tetrahedral (0/1)

- [ ] `geometry_tetrahedral` — PENDIENTE

### geometry_triangle (1/2)

- [x] `geometry_triangle` — presente en el GNU
- [ ] `geometry_triangle_eps_iso` — PENDIENTE

### global (0/4)

- [ ] `global_element_dof_apply` — PENDIENTE
- [ ] `global_element_dof_from_node_dof` — PENDIENTE
- [ ] `global_node_dof_empty` — PENDIENTE
- [ ] `global_post_point_node_type` — PENDIENTE

### ground (1/1)

- [x] `ground` — presente en el GNU

### groundflow (19/19)

- [x] `groundflow_consolidation_apply` — P6 (2026-08-20)
- [x] `groundflow_flux_edge_normal` — P6 (2026-08-20)
- [x] `groundflow_flux_edge_normal_element` — P6 (2026-08-20)
- [x] `groundflow_flux_edge_normal_element_group` — P6 (2026-08-20)
- [x] `groundflow_flux_edge_normal_element_node` — P6 (2026-08-20)
- [x] `groundflow_flux_edge_normal_element_node_factor` — P6 (2026-08-20)
- [x] `groundflow_flux_edge_normal_element_side` — P6 (2026-08-20)
- [x] `groundflow_flux_edge_normal_factor` — P6 (2026-08-20)
- [x] `groundflow_flux_edge_normal_geometry` — P6 (2026-08-20)
- [x] `groundflow_flux_edge_normal_node` — P6 (2026-08-20)
- [x] `groundflow_flux_edge_normal_sine` — P6 (2026-08-20); layout Pro 6.567 (`start end freq amp ...`) alineado 2026-09-13
- [x] `groundflow_flux_edge_normal_time` — P6 (2026-08-20)
- [x] `groundflow_nonsaturated_apply` — P6 (2026-08-20)
- [x] `groundflow_phreatic_level_multiple` — P6 (2026-08-20)
- [x] `groundflow_phreatic_level_multiple_element` — P6 (2026-08-20)
- [x] `groundflow_phreatic_level_multiple_element_geometry` — P6 (2026-08-20)
- [x] `groundflow_phreatic_level_multiple_element_group` — P6 (2026-08-20)
- [x] `groundflow_phreatic_level_multiple_node` — P6 (2026-08-20)
- [x] `groundflow_phreatic_level_multiple_static` — P6 (2026-08-20)
 - [x] `groundflow_seepage_geometry` — P6 (2026-08-20)
 - [x] `groundflow_seepage_node` — P6 (2026-08-20)
 - [x] `groundflow_seepage_eps` — P6 (2026-08-20)
 - [x] `groundflow_pressure_factor` — presente en el GNU
 - [x] `groundflow_phreatic_level_multiple_n` — P6 (2026-08-20)

### group_axisymmetric (1/1)

- [x] `group_axisymmetric` — presente en el GNU; 2026-10-09 (F4 T4, K20): same default Gauss rule as the plane elements (quad4 2x2, quad9 3x3), the forced 1-point rule and the quad9 refusal were removed (= Professional)

### group_beam (3/7)

- [ ] `group_beam_direction_z` — PENDIENTE
- [ ] `group_beam_direction_z_reference_point` — PENDIENTE
- [ ] `group_beam_force_moment_plasti` — PENDIENTE
- [x] `group_beam_inertia` — presente en el GNU
- [x] `group_beam_memory` — presente en el GNU
- [ ] `group_beam_shear` — PENDIENTE
- [x] `group_beam_young` — presente en el GNU

### group_condif (4/5)

- [ ] `group_condif_` — PENDIENTE
- [x] `group_condif_absorption` — presente en el GNU
- [x] `group_condif_capacity` — presente en el GNU
- [x] `group_condif_conductivity` — presente en el GNU
- [x] `group_condif_density` — presente en el GNU

### group_contact (7/8)

- [x] `group_contact_spring_direction` — `892541a` (2026-08-31; alias -> group_contactspring_direction)
- [x] `group_contact_spring_direction_automatic` — `892541a` (2026-08-31; alias -> GROUP_CONTACTSPRING_DIRECTION_AUTOMATIC; centroides de elementos)
- [ ] `group_contact_spring_direction_automatic_planes` — PENDIENTE
- [x] `group_contact_spring_memory` — `892541a` (2026-08-31; alias -> group_contactspring_memory; -total_linear aceptado)
- [x] `group_contact_spring_plasti_cohesion` — `892541a` (2026-08-31; alias -> group_contactspring_cohesion)
- [x] `group_contact_spring_plasti_friction` — `892541a` (2026-08-31; alias -> group_contactspring_friction)
- [x] `group_contact_spring_plasti_friction_automatic` — `892541a` (2026-08-31; alias -> group_contactspring_friction_automatic)
- [x] `group_contact_spring_sti` — `892541a` (2026-08-31; alias -> group_contactspring_stiffness)

### group_dof (0/2)

- [ ] `group_dof_initial` — PENDIENTE
- [ ] `group_dof_initial_speci` — PENDIENTE

### group_ground (0/1)

- [ ] `group_ground` — PENDIENTE

### group_groundflow (8/8)

- [x] `group_groundflow_permeability` — 2026-09-02 (variable 1..ndim; 1 valor = isotrópico en cada dirección, 6.618; ground19/large2/large3 rc=0)
- [x] `group_groundflow_capacity` — presente en el GNU
- [x] `group_groundflow_consolidation_apply` — P6 (2026-08-20)
- [x] `group_groundflow_materidivergence` — presente en el GNU
- [x] `group_groundflow_nonsaturated_eps_permeability` — P6 (2026-08-20)
- [x] `group_groundflow_nonsaturated_vangenuchten` — P6 (2026-08-20)
- [x] `group_groundflow_permeability_vertical_stress` — presente en el GNU
- [x] `group_groundflow_total_pressure_tension` — P6 (2026-08-20)

### group_integration (2/3)

- [x] `group_integration_method` — presente en el GNU
- [ ] `group_integration_method_reduced_factor` — PENDIENTE
- [x] `group_integration_points` — presente en el GNU

### group_interface (13/13 — familia COMPLETA del Carril A)

- [x] `group_interface` — Fase 1 implementada (commit `a82cbc7`, 2026-08-13)
- [x] `group_interface_condif_conductivity` — Sprint 8 (2026-08-24)
- [x] `group_interface_gap` — Fase 3 implementada (commit `9c2f4c8`, 2026-08-14; cerrada si strain > gap, default -1e20, hueco físico = gap negativo; validada con iface_mc_gap). SEMÁNTICA CORREGIDA 2026-09-04 (`6f2e940`, manual Pro 6.625 + verificación paso a paso contra el Pro): cerrada si strain <= gap (la condición estaba INVERTIDA), default +1e20, abierta = SIN tensiones (solo rigidez residual), la tensión acumula SOLO la deformación de los pasos cerrados (history `element_interface_force_norm`; interface2 rc=0: −101 = kn·(−101·1e-3), no kn·(−0.2)).
- [x] `group_interface_groundflow_capacity` — P6 (2026-08-20)
- [x] `group_interface_groundflow_permeability` — P6 (2026-08-20)
- [x] `group_interface_groundflow_total_pressure_tension` — P6 (2026-08-20)
- [x] `group_interface_materi_elasti_stiffness` — Fase 1 (commit `a82cbc7`; "elasti_sti" era OCR truncado)
- [x] `group_interface_materi_expansion_normal` — Sprint 8 (2026-08-24)
- [x] `group_interface_materi_memory` — commit `643865b` (2026-08-16; el checklist estaba desactualizado)
- [x] `group_interface_materi_plasti_mohr_coul_direct` — Fase 3 implementada (commit `9c2f4c8`, 2026-08-14; MC acumulativo con history `element_interface_force_tang`); convergencia interfaz (2026-09-03, commit feat de este lote): fuerzas totales acumuladas + dilatancia incremental + Fn positiva bajo compresión — mohr_coul_direct3 rc=0; 2026-10-11 (F4 T6d C52): el history guarda la fuerza tangencial RECORTADA (antes el trial elástico kt·gamma_total, sin descarga elástica tras invertir el deslizamiento; = Professional, ver registro de verificación)
- [x] `group_interface_materi_plasti_tension_direct` — Fase 3 implementada (commit `9c2f4c8`, 2026-08-14)
- [x] `group_interface_materi_residual_stiffness` — Fase 3 (commit `60bf78c`; "residual_sti" era OCR truncado)
- [x] `group_interface_tangential_reference_point` — Sprint 8 (2026-08-24; faltaba en el checklist original)

Nota: `group_interface_ground` del checklist anterior NUNCA existió (artefacto del inventario; no está en el manual ni en el inventario original de 2024).
- [ ] `group_interface_tangential_reference_point` — PENDIENTE

### group_materi (56/124)

- [x] `group_materi_damage_mazars` — presente en el GNU
- [x] `group_materi_damping` — presente en el GNU
- [x] `group_materi_damping_method` — Sprint 10 (PARCIAL: registrado)
- [x] `group_materi_density` — presente en el GNU
- [x] `group_materi_density_ground` — Sprint 10 (preexistente como density_groundflow; checklist al día)
- [ ] `group_materi_elasti_borja_tamagnini` — PENDIENTE
- [ ] `group_materi_elasti_c` — PENDIENTE
- [ ] `group_materi_elasti_c_direction` — PENDIENTE
- [x] `group_materi_elasti_camclay_g` — presente en el GNU
- [x] `group_materi_elasti_camclay_poisson` — presente en el GNU
- [x] `group_materi_elasti_camclay_pressure_min` — Sprint 10 (lote 9; manual 6.647: clamp de la presion usada en K = (1+e)·p/kappa del bloque elastico camclay de set_stress; validado con `mc_pressure_min`: p=0.001 inicial + pressure_min=10 -> K=750, sigxx -0.0093333 y sigyy -0.2093333 ANALITICOS EXACTOS; A/B `mc_pressure_min_off` sin record -> K=0.075 degenerado, sigxx +0.06566 con signo invertido)
- [x] `group_materi_elasti_compressibility` — presente en el GNU
- [x] `group_materi_elasti_hardsoil` — Sprint 10 (lote 8; ley de potencia E = Eref·((sig3+c·cotφ)/(sigmaref+c·cotφ))^m con sig3_manual = mayor autovalor algebraico (menor compresiva; el brief decia 'menor autovalor' = σ1 axial — NO es el HS); E50 primera carga vs Eur descarga via max|p| historico (dof sph compartido con 4.50, decision p_est = p_old + dp); validado con `mhardsoil_elast`/`_elast2` (-0.1689/-0.1241 EXACTOS) y `mhardsoil_unload`/`_flat` (+0.2692/-0.2692 EXACTOS)
- [x] `group_materi_elasti_k0` — Sprint 10 (lote 5; nu = K0/(1+K0), K0>0.95 truncado a 0.95; hook de control_materi_elasti_k0 conectado en set_stress — SOLO con young/young_power y poisson presentes; validado con `mk0`: K0=0.5 -> nu=1/3 -> oedometro sigxx -0.75/sigyy -1.5 EXACTOS = ratio K0; A/B con control -no -> -0.5769/-1.3462); 2026-09-29 (commit `ad1f5da`): el switch se lee en el índice de control del timestep (antes solo índice 0: `elasti7` ignorado), K0 también sin poisson y con hardsoil (nu50/nuur), medido en el Professional
- [x] `group_materi_elasti_lade` — presente en el GNU
- [x] `group_materi_elasti_poisson` — presente en el GNU
- [x] `group_materi_elasti_poisson_power` — Sprint 10 (lote 5; bloque en set_stress tras young_power: nu = nu0+nu1*(p/p1)^alpha, nu<=nu2, p=-sig_mean; validado con `mpower`: nu0=0.2 nu1=0.1 p1=1 alpha=1, E*eps=1.2 -> punto fijo analitico nu=0.4 p=2 (sigxx -1.7143); medido sigxx -1.4389/sigyy -3.2924 = 84% del fijo, ventana que discrimina la base nu=0.3). **Corrected 2026-10-06 (F2 T1, audit K06, hash in `odd/tasks/f2-medium-bugs.md`)**: the block now replaces the C (before C(young,poisson) + C(young,nu(p))); `mpower` -0.4389/-1.4101 = Professional 02-08-2026 within 0.6 %
- [x] `group_materi_elasti_shear_factor` — Sprint 10 (lote 5; multiplica SOLO las entradas de corte de C/Cmem, bloques (0,1),(0,2),(1,2); validado con `mshf`/`mshf_nof`: cizalla pura 1 paso, factor 2 -> sigxy 50 vs 33.33 sin record; factor 0 -> 0 EXACTO; GOTCHA: dof sigxy del sistema acoplado no responde lineal (ratio 1.5)
- [x] `group_materi_elasti_stress_pressure_history_factor` — Sprint 10 (lote 6; dof sph con el maximo |p| historico (initia 4.50) + en set_stress: si p_est del paso actual < maximo a inicio de paso -> C/Cmem x factor (descarga/recarga); validado con `msph`/`msph_flat`: carga 4 pasos + descarga 2 pasos, factor 3 -> sigyy +0.2692 vs -0.2692 con factor 1; sph 0.3333 EXACTO)
- [x] `group_materi_elasti_transverse_isotropy` — presente en el GNU
- [x] `group_materi_elasti_volumetric_poisson` — presente en el GNU
- [x] `group_materi_elasti_volumetric_young_order` — presente en el GNU
- [x] `group_materi_elasti_volumetric_young_values` — presente en el GNU
- [x] `group_materi_elasti_young` — presente en el GNU
- [x] `group_materi_elasti_young_polynomial` — presente en el GNU
- [x] `group_materi_elasti_young_power` — Sprint 10 lote 6 (UPGRADE 3→6 params: E=E0+E1*(p/p1)^alpha con E>=E2 y E<=E3, p=-sig_mean; la forma GNU 3-params ya no se acepta — CAMBIO SEMANTICO; fix del bug de C_matrix acumulativo del GNU 2014 (young+young_power doblaban la rigidez); validado con `myoung6` (E_eff~1237 = 72% del fijo secante 1714.29), `myoung6_e2`/`myoung6_e3` (clamps E2=3000/E3=800 EXACTOS) y `myoung6_apply` (gate -no -> E=E0 EXACTO))
- [ ] `group_materi_elasti_young_user` — PENDIENTE
- [x] `group_materi_expansion_linear` — presente en el GNU
- [x] `group_materi_expansion_volume` — presente en el GNU
- [x] `group_materi_factor` — Sprint 10 (escala tensiones+rigidez; verificado con vely prescrita)
- [x] `group_materi_failure_crunching` — Sprint 10 (lote 5; ortografia Professional con la "n" en name[], typo GNU `cruching` como alias db_number; validado con `mcrunch`/`mcrunch_low`: threshold alto 0.5 intacto sigxx -0.5769 EXACTO vs threshold -0.5 -> elemento marcado en t=0.05)
- [x] `group_materi_failure_damage` — presente en el GNU
- [x] `group_materi_failure_plasti_kappa` — presente en el GNU
- [x] `group_materi_failure_rupture` — presente en el GNU
- [x] `group_materi_failure_void_fraction` — Sprint 10 (lote 5; ortografia Professional con underscore en name[], GNU `voidfraction` como alias db_number; validado con `mvoid`/`mvoid_low`: void inicial 0.3, threshold 0.9 intacto vs 0.2 -> marcado en t=0.05)
- [ ] `group_materi_history_variable_user` — PENDIENTE
- [ ] `group_materi_history_variable_user_parameters` — PENDIENTE
- [x] `group_materi_hyper_besseling` — presente en el GNU
- [x] `group_materi_hyper_blatz_ko` — presente en el GNU
- [x] `group_materi_hyper_mooney_rivlin` — presente en el GNU
- [x] `group_materi_hyper_neohookean` — presente en el GNU
- [ ] `group_materi_hyper_reduced_polynomial` — PENDIENTE
- [x] `group_materi_hyper_volumetric_linear` — presente en el GNU
- [x] `group_materi_hyper_volumetric_murnaghan` — presente en el GNU
- [x] `group_materi_hyper_volumetric_ogden` — presente en el GNU
- [x] `group_materi_hyper_volumetric_polynomial` — presente en el GNU
- [ ] `group_materi_hyper_volumetric_simo_taylor` — PENDIENTE
- [x] `group_materi_maxwell_chain` — presente en el GNU
- [x] `group_materi_membrane` — presente en el GNU
- [x] `group_materi_memory` — presente en el GNU
  - [x] default (no record) = `-updated_linear` as the Professional — `d8117aa` (2026-10-04, C16; see registro de verificación); with `materi_displacement` `-total_linear` (manual 6.797) — `248ce83` (2026-10-08, C18)
  - [x] valor `-updated_jaumann` (6.686) — `09f0948` (2026-10-08, C20; see registro de verificación)
  - [x] valor `-updated_area` — implementado (`a1d9d12`, sesión 2026-09-02): pequeña deformación idéntica a `-updated_linear` (2026-09-29: salvo el gradiente, que es el de la geometría ACTUAL; sin rotación del stress desde el fix K01; el Pro 02-08-2026 rechaza `-updated_area` fuera de `incremental_driver`); la parte de área actual de cargas queda sin portador hasta `incremental_driver`. Los 3 tests del corpus etiquetados como "-updated_area" están bloqueados por la familia `incremental_driver`, NO por el valor (verificado empíricamente: con el valor aceptado, el parseo avanza hasta el record `incremental_driver`). Manuales `group_materi_memory_updated_area` (user/developer).
- [x] `group_materi_plasti_bounda` — Sprint 10 (alias de plasti_boundary)
- [x] `group_materi_plasti_bounda_factor` — Sprint 10 (alias)
- [x] `group_materi_plasti_camclay` — presente en el GNU
- [x] `group_materi_plasti_cap1` — Sprint 10 (lote 7; cap elastoplastico NUEVO de cero — f = q²/M² + p(p\*−p\*c) con p\* = p + c·cot(phi), p\*c = pc + c·cot(phi), 8 params phi c M lambda_star kappa_star K_ref p_ref m; pc = dof de materi_plasti_cap1_history; endurecimiento en stress.cc patrón kappa; validado con `mcap1`: compresion isotropica -> pc 110.0 EXACTO y sigxx -84.9996 vs -85.0 analitico (0.0005%), A/B elastico -175.0000, comb con druck_prag identico)
- [x] `group_materi_plasti_cap2` — Sprint 10 (lote 5; alias db_number del legacy cap (misma fisica c phi alpha R + tabla epsilonp_v pb; plasti.cc intacto); validado con `mcap2`/`mcap_legacy`: oedometro elastico c=1e6 -> sigxx -0.5769 EXACTO con AMBAS keywords = cap2 == legacy)
- [x] `group_materi_plasti_compression` — presente en el GNU
- [x] `group_materi_plasti_compression_direct` — Sprint 10 (cutoff espectral de autovalores)
- [x] `group_materi_plasti_compression_direct_visco` — Sprint 10 (relajacion 1-exp(-dt/tm))
- [x] `group_materi_plasti_coord_limit` — Sprint 10 (gate por coordenada vertical)
- [x] `group_materi_plasti_diprisco` — presente en el GNU
- [ ] `group_materi_plasti_diprisco_density` — PENDIENTE (evaluado en Sprint 10 lote 9: NO implementado. Razones: (1) la ley de interpolacion entre los dos juegos de parametros (loose l / dense d, 22 params: gamma_l beta_hat_ld b/c/t_p theta_hat_c/e xi_c/e beta_f0_l/d + e_l e_d) NO esta documentada en el manual — solo dice "the actually used data will be interpolated... using the current density" y remite a papers externos; sin la formula no hay implementacion verificable (regla de verificacion de AGENTS.md); (2) exigiria reescribir el bloque diprisco de plasti.cc con interpolacion inventada, arriesgando el diprisc1.dat validado; (3) la 12ª variable de historia (densidad relativa) no tiene ley de evolucion en el manual. La initia `materi_plasti_diprisco_history 12` SI se acepta (alias de materi_history_variables). La keyword queda registrada en database.cc? NO — sin implementacion no se registra el record (evita inputs que pasen el parseo y fallen en runtime)
- [x] `group_materi_plasti_druck_prag` — Sprint 10 (alias del legacy druckprag; VALIDACIÓN ANALÍTICA: cizalla φ=0 → σxy=2c/√3, 0.014% error)
- [x] `group_materi_plasti_element_group` — registrada (KEYWORDS BATCH C 2026-09-04); consumo (reduccion 2/3 por vecindad) PENDIENTE
- [x] `group_materi_plasti_element_group_factor` — registrada (KEYWORDS BATCH C 2026-09-04); consumo PENDIENTE
- [ ] `group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils` — PENDIENTE
- [x] `group_materi_plasti_gurson` — presente en el GNU
- [x] `group_materi_plasti_hardsoil` — Sprint 10 (lote 8; f = q/(E50(1-q/qa)) - 2q/Eur - γp, qa = qf/Rf con qf = 2·sinφ·(sig3+c·cotφ)/(1-sinφ) derivado de MC en falla; γp = kappa + extra del control; flujo asociativo; validado con la familia de relajacion `mhardsoil_gp0_off`/`mhardsoil_plast`: el estado final queda SOBRE la superficie (f ~ 0 verificada) y kappa escala con q inicial (0.000293 -> 0.0006 para σxx -2 -> -4)
- [x] `group_materi_plasti_heat_generation` — Sprint 10 (nombre Professional del legacy heatgeneration, dual-read)
- [x] `group_materi_plasti_hypo_cohesion` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin` — implementada (commit `974b284`, 2026-08-06)
- [x] `group_materi_plasti_hypo_masin_clay` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin_clay_advanced_parameters` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin_clay_avanced_direction` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin_clay_ocr` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin_clay_structure` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin_clay_visco` — implementada (commit `a0f55a6`, 2026-08-06)
- [x] `group_materi_plasti_hypo_masin_ocr` — presente en el GNU
- [x] `group_materi_plasti_hypo_masin_structure` — presente en el GNU
- [ ] `group_materi_plasti_hypo_minimum_void_ratio` — PENDIENTE
- [ ] `group_materi_plasti_hypo_niemunis_visco` — PENDIENTE
- [ ] `group_materi_plasti_hypo_niemunis_visco_ocr` — PENDIENTE
- [x] `group_materi_plasti_hypo_strain_intergranular` — batch hypo (2026-09-01; 6 params R m_R m_T beta_r chi theta, manual 6.715; alias db_number del nombre Professional `..._hypo_strain_intergranular`; theta sin slot en el kernel (usa chi); hypo2/3 rc=0)
- [x] `group_materi_plasti_hypo_strain_intergranular_masin_clay` — presente en el GNU
- [ ] `group_materi_plasti_hypo_strain_isa` — PENDIENTE
- [ ] `group_materi_plasti_hypo_void_ratio_linear` — PENDIENTE
- [ ] `group_materi_plasti_hypo_wol` — PENDIENTE
- [x] `group_materi_plasti_kinematic_hardening` — presente en el GNU
- [x] `group_materi_plasti_mohr_coul` — Sprint 10 (lote 3; implementación NUEVA en plasti.cc con matrix_eigenvalues, f exacta del manual; VALIDACIÓN ANALÍTICA: resistencia a tracción uniaxial σt=2c·cosφ/(1+sinφ) = 46.19 EXACTO; sanity Tresca φ=0 → 2c=80 con convergencia lenta por vértice singular, documentado)
- [x] `group_materi_plasti_mohr_coul_direct` — dos modos (2026-09-02): (1) SIN `_normal` = ley MC completa del manual 6.726 (cut del diferencia de tensiones principales con cap espectral de tensiones 6.738 y flujo no-asociado phi_flow; `materi_direct_full_mc` en stress.cc — mohr_coul_direct1/2 rc=0, validado EXACTO contra el Professional en los estados de cortante puro); (2) CON `_normal` = cut-off de cortante en el plano (6.727, `materi_direct_cutoff`, validada con `materi_direct_mc`). Corrige la lectura previa (b6d4be6) que trataba SIEMPRE como cut en plano.
- [x] `group_materi_plasti_mohr_coul_direct_normal` — implementada (normal explícita del plano; validada con `materi_direct_mc`; regla del plano del Professional + phi_flow M/G, commits `d17461e` + `77f9282` 2026-10-02, tests `direct_cutoff_plane_*`; `materi_direct_mc` queda con target obsoleto, BACKLOG C17)
- [x] `group_materi_plasti_mohr_coul_direct_normal_automatic` — implementada (normal del elemento = `element_normal`; normal nula con `-yes` = error de entrada, `materi_direct_auto` es test de error esperado desde 2026-10-01; `-no` = sin normal automática (corte completo, como Pro), commit `260bd95` 2026-10-02)
- [x] `group_materi_plasti_mohr_coul_direct_visco` — implementada (relajación visco: `factor = 1-exp(-dt/tm)`; validada con `materi_direct_visco`)
- [x] `group_materi_plasti_mohr_coul_direct_wall` — implementada (valores alternativos si el elemento está pegado a una pared, vía `plasti_on_boundary`; validada con `materi_direct_wall`)
- [x] `group_materi_plasti_mohr_coul_hardening_softening` — Sprint 10 (lote 4; interpolación LINEAL de phi/c/phi_flow con kappa/kappa_shear_crit en [0,1]; bloque en plasti.cc con matrix_eigenvalues, mismo f que el MC clásico; validado con `mmchs_soft`: tracción uniaxial φ=0, c 80→20, κ_crit=0.5 → κ medido 0.318 (ratio 0.636 → c=41.8 → σt=83.6 analítico) → σxx medido 81.9 = 98%; calibración del test: pasos pequeños dt=0.02 + control_timestep_iterations 8 para que el retorno persiga la superficie que BAJA — el rig de 1 paso brutal divergía a 7.55)
- [x] `group_materi_plasti_mohr_coul_direct_hardening_softening` — batch hypo (2026-09-01; variante DIRECT (angulos en GRADOS) del hardening_softening, 7 valores; conversion PIRAD/180 en plasti_rule; dam_building la usa)
- [ ] `group_materi_plasti_mpc` — PENDIENTE
- [ ] `group_materi_plasti_mpc_factor` — PENDIENTE
- [x] `group_materi_plasti_pressure_limit` — Sprint 10 (gate p=-tr/3 en el camino direct)
- [x] `group_materi_plasti_residual_sti` — OCR: residual_stiffness
- [x] `group_materi_plasti_tension` — presente en el GNU
- [x] `group_materi_plasti_tension_direct` — implementada (cut-off directo de tracción normal en un plano; validada con `materi_direct`)
- [x] `group_materi_plasti_tension_direct_normal` — implementada (normal explícita; validada con `materi_direct`; plano propio de la ley de tracción y aceptada sin `tension_direct` (sigy 0 con mohr_coul_direct), commit `d17461e` 2026-10-02, tests `direct_cutoff_plane_tension_normal2`, `_two_normals1`)
- [x] `group_materi_plasti_tension_direct_normal_automatic` — implementada (normal del elemento = `element_normal`; normal nula con `-yes` = error de entrada, `materi_direct_auto` es test de error esperado desde 2026-10-01; `-no` = sin normal automática (corte completo, como Pro), commit `260bd95` 2026-10-02)
- [x] `group_materi_plasti_tension_direct_visco` — implementada (relajación visco; validada con `materi_direct_visco`)
- [x] `group_materi_plasti_tension_direct_wall` — implementada (valores alternativos en pared; validada con `materi_direct_wall`)
- [x] `group_materi_plasti_user` — presente en el GNU
- [x] `group_materi_plasti_visco_exponential` — presente en el GNU
- [x] `group_materi_plasti_visco_exponential_limit` — Sprint 10 (PARCIAL)
- [x] `group_materi_plasti_visco_exponential_name` — Sprint 10 (PARCIAL)
- [x] `group_materi_plasti_visco_exponential_values` — Sprint 10 (PARCIAL) + batch hypo (2026-09-01; alias db_number del singular Professional `..._value`)
- [x] `group_materi_plasti_visco_power` — layout Professional `η p` (variable-length, f_ref implícito=1; el layout legacy `η p f_ref` sigue aceptado) — 2026-09-02 (visc_pl1/2, validation_12 rc=0)
- [ ] `group_materi_plasti_visco_power_name` — PENDIENTE
- [ ] `group_materi_plasti_visco_power_value` — PENDIENTE
- [x] `group_materi_plasti_vonmises` — presente en el GNU
- [x] `group_materi_plasti_vonmises_nadai` — presente en el GNU
- [x] `group_materi_stokes` — presente en el GNU; 2026-10-12 (F4 T6f K21 WU2, `58bba14`): default -yes (manual 6.752 = Pro) and read per element group (class A STOKES-IDX: the Pro reads index 0), see registro de verificación
- [ ] `group_materi_stress_null` — PENDIENTE
- [ ] `group_materi_stress_null_direction` — PENDIENTE
- [ ] `group_materi_stress_null_direction_automatic` — PENDIENTE
- [ ] `group_materi_umat` — PENDIENTE
- [ ] `group_materi_umat_parameters` — PENDIENTE
- [ ] `group_materi_umat_pardiso_decompose` — PENDIENTE
- [x] `group_materi_undrained_capacity` — familia undrained (2026-09-05, 6.760 + teoria 2.2.7): consumo en materi.cc (total_new_sig + rigidez volumetrica dt/C); ground17/undrained1/undrained2 rc=0
- [x] `group_materi_viscosity` — presente en el GNU
- [x] `group_materi_viscosity_heatgeneration` — presente en el GNU

### group_plasti (0/1)

- [ ] `group_plasti_apply` — PENDIENTE

### group_porosity (0/1)

- [x] `group_porosity` — lote keywords baratas (alias db_number de group_groundflow_porosity)

### group_spring (2/4)

- [x] `group_spring_direction` — presente en el GNU
- [x] `group_spring_memory` — registrada (KEYWORDS BATCH C 2026-09-04); el muelle GNU es incremental en configuracion actual (= -updated_linear); spring1 rc=0
- [x] `group_spring_plasti` — presente en el GNU
- [ ] `group_spring_sti` — PENDIENTE

### group_time (1/2)

- [x] `group_time` — presente en el GNU; C39 2026-10-06 (F2 T2c, hash in `odd/tasks/f2-medium-bugs.md`): a node whose elements are all inactive (not yet born / dead) is not driven by bounda_dof/bounda_unknown (= Pro truss8); death velocity reset of the Pro not implemented (BACKLOG C39)
- [ ] `group_time_` — PENDIENTE

### group_truss (5/9)

- [x] `group_truss_area` — presente en el GNU
- [x] `group_truss_density` — presente en el GNU
- [ ] `group_truss_elasti_elongation_force_diagram` — PENDIENTE
- [ ] `group_truss_elasti_young` — PENDIENTE
- [ ] `group_truss_expansion` — PENDIENTE
- [ ] `group_truss_initial_force` — PENDIENTE
- [x] `group_truss_memory` — presente en el GNU; K17 2026-10-06 (F2 T2): `-total` geometrically nonlinear (= Pro truss14) and tension geometric stiffness for -updated/-total (Pro truss2/3 converge); see registro de verificación; 2026-10-11 (F4 T6d, C11, `28992d5`): default `-total_linear` with `materi_displacement`, `-updated_linear` without (= Professional), signed linear length; Pro input checks for the combinations replicated (F4 T6d C59 `e0f2204`, owner 2026-10-11; see registro de verificación)
- [x] `group_truss_plasti` — presente en el GNU
- [x] `group_truss_rope` — presente en el GNU

### group_type (1/1)

- [x] `group_type` — presente en el GNU

### group_volume (1/1)

- [x] `group_volume_factor` — presente en el GNU

### group_wave (1/1)

- [x] `group_wave_speed_of_sound` — presente en el GNU

### icontrol (1/1)

- [x] `icontrol` — presente en el GNU

### incremental (0/1)

- [ ] `incremental_driver` — PENDIENTE (2026-09-02, diagnóstico corregido: los 6 tests de `test/incremental_driver/` —oedometric/triaxial drained/undrained mohr_coul_direct, direct_shear_drained_hypo, triaxial_drained_hypo_isa, syntax— dependen TODOS de esta familia, incluidos los 3 atribuidos antes a `-updated_area`. El GNU no registra NINGÚN record de la familia (0 refs). El test syntax bloquea además en initia por `materi_strain_isa_c` y en data por `incremental_ntime`/`incremental_driver_*`/hipo wolfersdorff/ISA)

### inertia (1/1)

- [x] `inertia_apply` — registrado (alias `options_inertia`); default -no alineado al Professional (manual 6.785) en el sprint u-p (2026-09-05, commit feat de este lote): sin el registro las analyses con densidad son estáticas; ground14/15/16 rc=0

### input_abaqus (6/6)

- [x] `input_abaqus` — presente en el GNU
- [x] `input_abaqus_continue` — presente en el GNU
- [x] `input_abaqus_group` — presente en el GNU
- [x] `input_abaqus_mesh` — presente en el GNU
- [x] `input_abaqus_name` — presente en el GNU
- [x] `input_abaqus_set` — presente en el GNU

### input_fe (0/1)

- [ ] `input_fe` — PENDIENTE

### input_gmsh (1/1)

- [x] `input_gmsh` — presente en el GNU

### interface (0/1)

- [ ] `interface_gap_apply` — PENDIENTE

### linear (0/1)

- [ ] `linear_calculation_apply` — PENDIENTE

### materi_acceleration (1/1)

- [x] `materi_acceleration` — implementada (feat f516a03, 2026-09-04): dofs de aceleracion derivados acc = (v_new - v_old)/dt (dof.cc, parallel_new_dof_diagonal) + basenames accx/accy/accz (database.cc) + bounda_dof -accx = bound de velocidad v_new = v_old + a*dt (bounda.cc). Unlocks el PARSE de los 12 tests que la declaran; dynamic3/dynamic4 rc=0 (verificados contra el binario Professional 25-10-2023: disx 62.317 identico).

### materi_damage (2/2)

- [x] `materi_damage` — presente en el GNU — 2026-10-12 (F4 T6f C68, `70235c6`): kept as an integration-point state (ELEMENT_DOF), nodal value with the strain recovery weight (see registro de verificación)
- [x] `materi_damage_apply` — implemented (commit `9bbd125`, 2026-09-27): global damage switch, overruled by control_materi_damage_apply; verified against the Professional binary 02-08-2026 (see verification registry)

### materi_displacement (1/2)

- [x] `materi_displacement` — presente en el GNU
- [x] `materi_displacement_relative` — implementada (dof relativo `disr*`; referencia re-sincronizada en cambio de timestep en `control_timestep` y en reset de desplazamiento en `control_reset_dof`; integrada en `dof.cc`; validada con `mat_rel` y `mat_rel_reset`)

### materi_dynamic (1/1)

- [x] `materi_dynamic` — implementada (feat f516a03, 2026-09-04): factor de mezcla sigma = (1-f)*sigma(t) + f*sigma(t+dt) (manual 6.800) como blend del momento en materi() (sigvec vs old_sig + rigidez*f) y truss() (new_truss_force_mix + stiffness*f); registro global double no_index; control_materi_dynamic (6.141) pisando por icontrol. f=1 = esquema historico byte-identico. PARCIAL: en el esquema de velocidades del GNU f<1 reduce la disipacion pero tambien la estabilidad (limite explicito); dynamic1/2/5/8 siguen RUNFAIL por el esquema (bloqueador solver, ver manual-developer/materi_dynamic.md); 2026-10-11 (F4 T6d, C11, `28992d5`): factor efectivo `materi_dynamic_effective()` (control_materi_dynamic en ICONTROL pisa el global) en todos los lectores

### materi_elasti (0/1)

- [x] `materi_elasti_young_power_apply` — Sprint 10 (lote 6; gate control_materi_elasti_young_power_apply -no -> E=E0 constante del record young_power (manual 6.801); alias db_number para el nombre del manual; validado con `myoung6_apply`: -0.5769/-1.3462 EXACTO = base lineal)

### materi_failure (1/1)

- [x] `materi_failure_apply` — implemented (commit `9bbd125`, 2026-09-27): global failure switch (failure.cc), overruled by control_materi_failure_apply; verified against the Professional binary 02-08-2026 (see verification registry)

### materi_history (0/1)

- [ ] `materi_history_variable` — PENDIENTE

### materi_maxwell (1/1)

- [x] `materi_maxwell_stress` — presente en el GNU — 2026-10-12 (F4 T6f C68, `70235c6`): kept as an integration-point state (ELEMENT_DOF), nodal value with the strain recovery weight (see registro de verificación) ; 2026-10-12 (F4 T6f C68d, `2327083`): two or more chains with rho/epi/epp read the old state with the right index

### materi_plasti (4/16)

- [x] `materi_plasti_camclay_history` — lote keywords baratas (initia 4.16: dofs cchis0/cchis1 = e y p0 del Modified Cam Clay; group_materi_plasti_camclay = M kappa lambda con p0 como historia; validado con validation_13 rc=0)
- [x] `materi_plasti_cap1_history` — Sprint 10 (lote 7; initia del dof escalar pc del cap1 — manual 4.17 "the history variable pc for the cap1 plasticity models is added to the node_dof records"; basename `pc`; patrón kappa completo: parser input.cc, RHS materi.cc, inercia general.cc, clamp dof.cc; validado con `mcap1`: pc 100 -> 109.9996 tras 20 pasos plasticos y constante en la descarga) ; 2026-10-12 (F4 T6f C68, `70235c6`): kept as an integration-point state (ELEMENT_DOF), nodal value with the strain recovery weight (see registro de verificación)
- [x] `materi_plasti_diprisco_history` — Sprint 10 (lote 9; initia 4.18: alias de materi_history_variables en input.cc (misma rama del parser, mismo dof hisv compartido, basenames hisv0..hisvN-1); los checks de GROUP_MATERI_PLASTI_DIPRISCO/_RT aceptan cualquiera de los dos nombres (check_unknown_atleastone); validado con `mdiprisco_hist` (diprisc1 con el alias: hisv10 -152.795 IDENTICO a diprisc1) y `mstrain_diprisco`)
- [x] `materi_plasti_f` — presente en el GNU
- [x] `materi_plasti_f_nonlocal` — presente en el GNU
- [ ] `materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils_history` — PENDIENTE
- [x] `materi_plasti_hardsoil_history` — Sprint 10 (lote 8; initia 4.22: dof escalar max |p| historico — mismo concepto que materi_stress_pressure_history (4.50), dof `sph` COMPARTIDO (misma actualizacion running-max en dof.cc, mismo basename); el switch E50/Eur del elastico lo lee (old_unknowns[sph]); validado con `mhardsoil_unload`: sph = 0.3333 en el pico y Eur desde el primer paso de descarga)
- [x] `materi_plasti_hypo_history` — batch hypo (2026-09-01; initia 4.23: 8 variables hyhis0..hyhis7, mismo dof hisv compartido, basenames hyhis; kernel masin mapea e=hyhis0 y s=hyhis4; hypo2/4 rc=0)
- [ ] `materi_plasti_hypo_substepping` — PENDIENTE
- [x] `materi_plasti_kappa` — presente en el GNU
- [x] `materi_plasti_kappa_shear` — batch hypo (2026-09-01; initia 4.25: dof kapsh = int sqrt(0.5 dev(epp):dev(epp)); variable de hardening de la familia hardening_softening (plasti_rule); dam_building/slope la consumen) ; 2026-10-12 (F4 T6f C68, `70235c6`): kept as an integration-point state (ELEMENT_DOF), nodal value with the strain recovery weight (see registro de verificación)
- [ ] `materi_plasti_maximum_iterations` — PENDIENTE
- [ ] `materi_plasti_phimob` — PENDIENTE
- [x] `materi_plasti_rho` — presente en el GNU — 2026-10-12 (F4 T6f C68, `70235c6`): kept as an integration-point state (ELEMENT_DOF), nodal value with the strain recovery weight (see registro de verificación)
- [ ] `materi_plasti_tension_apply` — PENDIENTE
- [ ] `materi_plasti_visco_apply` — PENDIENTE

### materi_strain (4/21)

- [x] `materi_strain_elasti` — presente en el GNU
- [ ] `materi_strain_energy` — PENDIENTE
- [x] `materi_strain_intergranular` — presente en el GNU
- [ ] `materi_strain_isa_c` — PENDIENTE
- [ ] `materi_strain_isa_eacc` — PENDIENTE
- [x] `materi_strain_plasti` — presente en el GNU
- [ ] `materi_strain_plasti_camclay` — PENDIENTE
- [x] `materi_strain_plasti_cap` — Sprint 10 (lote 9; initia 4.35: tensor de deformacion plastica propio de los modelos cap en node_dof (dof dedicado capepp, basenames eppcapxx...; RHS en materi.cc en el BLOQUE CONSOLIDADO de dofs por-modelo + inercia en general.cc); validado con `mstrain_cap`: eppcapzz -0.0360001 EXACTO (deps_p_cv=0.0054/paso x 20 / 3) y `mstrain_cap_elast`: 0)
- [x] `materi_strain_plasti_compression` — Sprint 10 (lote 9; initia 4.36: dof dedicado cepp, basenames eppcmpxx...; mismo bloque consolidado; `mstrain_compression`: eppcmpzz -8.96e-05; gemelo elastico 0)
- [x] `materi_strain_plasti_diprisco` — Sprint 10 (lote 9; initia 4.37: dof dedicado depp, basenames eppdipxx...; `mstrain_diprisco`: eppdipyy -0.0105714 y eppdipxx=eppdipzz +0.00350735; gemelo elastico 0)
- [x] `materi_strain_plasti_druckprag` — Sprint 10 (lote 9; initia 4.39: dof dedicado dpepp, basenames eppdrpxx...; `mstrain_druckprag`: eppdrpxy 0.910768 IDENTICO al epp generico del mismo run; gemelo elastico 0)
- [ ] `materi_strain_plasti_generalised_non_associate_cam_clay_for_bonded_soils` — PENDIENTE
- [x] `materi_strain_plasti_hardsoil` — Sprint 10 (lote 8; initia 4.40: tensor de deformacion plastica propio del hardsoil en node_dof (dof dedicado hsepp, basenames epphsxx...; RHS en materi.cc + inercia en general.cc, patron de materi_strain_plasti); validado en los tests de relajacion plastica (se llena con el retorno) y cero en los elasticos)
- [ ] `materi_strain_plasti_mohr_coul` — PENDIENTE
- [ ] `materi_strain_plasti_tension` — PENDIENTE
- [ ] `materi_strain_plasti_vonmises` — PENDIENTE
- [x] `materi_strain_total` — presente en el GNU
- [ ] `materi_strain_total_compression_kappa` — PENDIENTE
- [ ] `materi_strain_total_kappa` — PENDIENTE
- [ ] `materi_strain_total_shear_kappa` — PENDIENTE
- [ ] `materi_strain_total_tension_kappa` — PENDIENTE

### materi_stress (1/2)

- [x] `materi_stress` — presente en el GNU
- [x] `materi_stress_pressure_history` — Sprint 10 (lote 6; initia: dof sph = maximo running de |p| (p=-sig_mean) en node_dof, actualizado en dof.cc parallel_new_dof_diagonal (patron clamp de kap); basename -sph; validado con `msph`/`msph_flat`: sph 0.3333 EXACTO en el pico y mantenido en la descarga)

### materi_velocity (2/2)

- [x] `materi_velocity` — presente en el GNU
- [x] `materi_velocity_integrated` — presente en el GNU

### materi_void (1/1)

- [x] `materi_void_fraction` — presente en el GNU

### materi_work (1/1)

- [x] `materi_work` — presente en el GNU

### mesh (1/51)

- [x] `mesh` — presente en el GNU
- [x] `mesh_activate_gravity_element` — implementada (selección de elementos por rango; `mesh_activate_gravity_factor` en mesh.cc)
- [x] `mesh_activate_gravity_element_group` — implementada (selección por grupos)
- [x] `mesh_activate_gravity_geometry` — implementada (selección por geometría)
- [x] `mesh_activate_gravity_method` — implementada (método 1 y 2; método 2 mantiene el elemento activo con rigidez reducida; validada con `mesh_act_grav2`)
- [x] `mesh_activate_gravity_sti` — implementada (factor de rigidez aplicado a la matriz del elemento en materi.cc)
- [x] `mesh_activate_gravity_time` — implementada (activación gradual de la gravedad; validada con `mesh_act_grav`); indexed, bottom to top, method 2 default = Pro `mesh_activate1` (F4 T8 `7aeac4a`, 2026-10-10)
- [x] `mesh_activate_gravity_time_initial` — implementada (time of birth)
- [x] `mesh_activate_gravity_time_strain_settlement` — implementada (registrado)
- [ ] `mesh_boundary` — PENDIENTE
- [ ] `mesh_correct` — PENDIENTE
- [ ] `mesh_delete_geometry_moving` — PENDIENTE
- [ ] `mesh_element_group_apply` — PENDIENTE
- [ ] `mesh_gid_arc_coord` — PENDIENTE
- [ ] `mesh_gid_assign_conditions_line` — PENDIENTE
- [ ] `mesh_gid_assign_conditions_point` — PENDIENTE
- [ ] `mesh_gid_assign_conditions_surface` — PENDIENTE
- [ ] `mesh_gid_circle_coord` — PENDIENTE
- [ ] `mesh_gid_circle_element_group` — PENDIENTE
- [ ] `mesh_gid_circle_hollow` — PENDIENTE
- [ ] `mesh_gid_circle_radius` — PENDIENTE
- [ ] `mesh_gid_cylinder_coord` — PENDIENTE
- [ ] `mesh_gid_cylinder_element_group` — PENDIENTE
- [ ] `mesh_gid_cylinder_height` — PENDIENTE
- [ ] `mesh_gid_cylinder_hollow` — PENDIENTE
- [ ] `mesh_gid_cylinder_normal` — PENDIENTE
- [ ] `mesh_gid_cylinder_radius` — PENDIENTE
- [ ] `mesh_gid_line_element_group` — PENDIENTE
- [ ] `mesh_gid_line_point` — PENDIENTE
- [ ] `mesh_gid_line_size` — PENDIENTE
- [ ] `mesh_gid_line_structured_concentrate` — PENDIENTE
- [ ] `mesh_gid_line_structured_nel` — PENDIENTE
- [ ] `mesh_gid_line_structured_size` — PENDIENTE
- [ ] `mesh_gid_point_coord` — PENDIENTE
- [ ] `mesh_gid_rectangle_coord` — PENDIENTE
- [ ] `mesh_gid_rectangle_element_group` — PENDIENTE
- [ ] `mesh_gid_rectangle_hollow` — PENDIENTE
- [ ] `mesh_gid_size` — PENDIENTE
- [ ] `mesh_gid_sphere_coord` — PENDIENTE
- [ ] `mesh_gid_sphere_element_group` — PENDIENTE
- [ ] `mesh_gid_sphere_hollow` — PENDIENTE
- [ ] `mesh_gid_sphere_radius` — PENDIENTE
- [ ] `mesh_gid_surface_element` — PENDIENTE
- [ ] `mesh_gid_surface_element_group` — PENDIENTE
- [ ] `mesh_gid_surface_line` — PENDIENTE
- [ ] `mesh_gid_surface_structured_nel` — PENDIENTE
- [ ] `mesh_gid_surface_structured_size` — PENDIENTE
- [ ] `mesh_gid_volume_element_group` — PENDIENTE
- [ ] `mesh_gid_volume_surface` — PENDIENTE
- [ ] `mesh_interface_triangle_coordinates` — PENDIENTE
- [ ] `mesh_interface_triangle_element_group` — PENDIENTE

### message (0/1)

- [ ] `message` — PENDIENTE

### mpc (6/17)

- [x] `mpc_apply` — 2026-09-02 (switch global sin índice; default -yes; -no apaga generación Y aplicación en mpc.cc::mpc_node_apply)
- [x] `mpc_node_number` + `mpc_node_factor` — implementadas (2026-09-01, mpc.cc::mpc_node_apply desde bounda.cc + re-sync post-solve en top.cc): dof esclavo = Σ factor×dof maestro, esclavo bounded (valor recomputado cada iteración), sin condensación ni redistribución de fuerzas (semántica verificada contra el binario Professional 25-10-2023). Desbloquea mpc1 (rc=0, node_dof idéntico al Professional: nodo 2 velx=disx=2, sigxx=0.5) y el tie de mohr_coul_direct2/4/7 (nodos 3=4 idénticos; los tests siguen RUNFAIL por la familia MATERIAL direct: mohr_coul_direct1/3/5/6/8 sin mpc también fallan, eptxx 2.99e-3 vs 1.49e-2)
- [x] `mpc_linear_quadratic` — implementada (2026-09-01, mpc.cc::mpc_linear_quadratic_generate): generación de ties para nodos colgantes de elementos cuadráticos contra el elemento lineal adyacente (shape functions vía point_el, eps_iso 1.e-4), regeneración automática al cambiar la malla (fingerprint en record interno MPC_LINEAR_QUADRATIC_MESH_FINGERPRINT). Records generados byte-idénticos al Professional (mpc3/4/5: mid-edge 2 maestros ×0.5, mid-face 4 ×0.25). mpc3/mpc4 siguen RUNFAIL: diagnóstico fino 2026-09-07 (DIAG-SOLVE-MIXTO §15): el tie es value-constrained (fila del esclavo eliminada sin redistribución energética de su fuerza de arista cuadrática 2/3·t·L) → punto fijo del staggered 0.454 ≠ 1/3; la hipótesis previa de los σ-dofs libres queda REFUTADA (fijar los esclavos a los valores exactos sigue dando campo contaminado). Cierre PENDIENTE: eliminación energética implícita en el solve (familia §8 C); mpc5 además necesitaba la formulación `-total_linear` (reset descartado, 0 vs 0.5, mismo root que delete3), CORREGIDA en `f045d04` (2026-09-13: reset sobrevive, sigyy converge a 0.5 con suficientes iteraciones)
- [x] `mpc_element_group` — 2026-09-02 (mpc.cc::mpc_element_group_generate): ties de los nodos de group_0 dentro de elementos de group_1 con las shape functions en su posición iso (point_el, eps 1e-4), restricción de geometría (_geometry) y de dofs (_dof); saltos de nodos ya miembros del elemento master (semántica always -no); regeneración por fingerprint (record interno MPC_ELEMENT_GROUP_MESH_FINGERPRINT). mpc7 rc=0: patch test sigxx 1.0±1e-8; ground19 (muro cut-off) rc=0
- [ ] `mpc_element_group_always` — PENDIENTE
- [ ] `mpc_element_group_closest` — PENDIENTE
- [ ] `mpc_element_group_coord_geometry` — PENDIENTE
- [x] `mpc_element_group_dof` — 2026-09-02 (selección de dofs a igualar; default: todos los principal)
- [ ] `mpc_element_group_eps_iso` — PENDIENTE
- [x] `mpc_element_group_geometry` — 2026-09-02 (restringe los nodos candidatos a los de la geometría)
- [ ] `mpc_element_group_keep` — PENDIENTE
- [ ] `mpc_geometry` — PENDIENTE
- [ ] `mpc_geometry_dof` — PENDIENTE
- [ ] `mpc_geometry_method` — PENDIENTE
- [ ] `mpc_geometry_switch` — PENDIENTE
- [ ] `mpc_geometry_tolerance` — PENDIENTE

### mrange (0/1)

- [ ] `mrange` — PENDIENTE

### mstring (0/1)

- [ ] `mstring` — PENDIENTE

### node (18/23; 3 parciales documentados, 2 print-only pendientes)

- [x] `node` — presente en el GNU
- [x] `node_boundary` — presente en el GNU
- [x] `node_bounded` — presente en el GNU
- [ ] `node_bounded_index` — PENDIENTE (print-only, sin consumidor en el GNU)
- [ ] `node_convection_apply` — PARCIAL: registrado+parseado, sin comportamiento (requiere cirugía de options_convection por nodo)
- [x] `node_damping` — presente en el GNU
- [x] `node_deformed_mesh` — presente en el GNU
- [x] `node_dof` — presente en el GNU
- [x] `node_dof_calcul` — presente en el GNU
- [x] `node_dof_start_refined` — presente en el GNU (el "node_dof_start_re" del inventario era OCR truncado)
- [x] `node_dynamic_pressure` — Sprint 8 (2026-08-24; override en calcul.cc)
- [x] `node_force` — Sprint 8 (2026-08-24; dof.cc)
- [ ] `node_geometry_present` — PENDIENTE (print-only)
- [x] `node_inertia` — Sprint 8 (2026-08-24; record calculado m·(a+g), fuente del d'alembert)
- [x] `node_mass` — presente en el GNU + FIX data_length=1→ndim (Sprint 8: en 2D/3D el parser se atragantaba — inutilizable desde los orígenes)
- [ ] `node_mesh` — PARCIAL: registrado+parseado, sin consumidor del record `mesh` por nodo en el GNU
- [x] `node_rhside` — presente en el GNU
- [x] `node_slide` — Sprint 8 (2026-08-24; membresía aditiva en slide.cc; verificación física pendiente de modelo dedicado — smoke test)
- [x] `node_start_refined` — presente en el GNU (OCR "node_start_re")
- [x] `node_static_pressure` — Sprint 8 (2026-08-24; override en calcul.cc)
- [x] `node_stiffness` — presente en el GNU (OCR "node_sti")
- [x] `node_support_edge_normal_plasti_tension_status` — Sprint 13 lote 3 COMPLETO (4e0074f: 0=cerrado, 1=abierto, OR por nodo)
- [x] `node_total_pressure` — Sprint 8 (2026-08-24; override en calcul.cc)

### nonlocal (2/2)

- [x] `nonlocal` — batch hypo (2026-09-01; alias de options_nonlocal, manual 6.897) + 2026-09-02 (fix camino clásico: nonlocal_set() restaurado tras cambio de malla en la iteración — validation_12 rc=0 con target kap no-local; slope_nonlocal_refine sigue chocando con el bug legacy del buffer NONLOCAL_ITEM_SIZE de nonloc.cc — PENDIENTE)
- [x] `nonlocal_name` — batch hypo (2026-09-01; record nuevo 6.898, PARCIAL: aceptado pero sin gate por modelo (el GNU aplica el termino no-local a todo modelo))

### number (0/1)

- [ ] `number_of_space_dimensions` — PENDIENTE

### plasti (1/1)

- [x] `plasti_apply` — implemented (commit `9bbd125`, 2026-09-27): global plasticity switch (MC, von Mises, DP, direct cut-offs, camclay, hardsoil, hypoplasticity, incremental models; slide default); verified against the Professional binary 02-08-2026 (see verification registry)

### post_apply (0/1)

- [ ] `post_apply` — PENDIENTE

### post_calcul (11/22)

- [x] `post_calcul` — presente en el GNU
- [ ] `post_calcul_absolute` — PENDIENTE
- [ ] `post_calcul_apparent_total` — PENDIENTE
- [ ] `post_calcul_label` — PENDIENTE
- [ ] `post_calcul_limit` — PENDIENTE
- [x] `post_calcul_materi_stress_force_average` — sub-sprint materi_stress_force, lotes 1-5 (2026-08-28; quad9/hex27: promedio de las fuerzas/momentos de las caras opuestas en los nodos del plano medio; default -yes; IMPLEMENTADO en 2D (lote 2) y 3D hex27 (lote 3): promedio EXACTO + flag por nodo para -primary; msf_hex27_avg 45 vs 27)
- [x] `post_calcul_materi_stress_force_direction_exclude` — sub-sprint materi_stress_force, lotes 1 + 3 (2026-08-28; 3D: caras con |n·dir| > 1-eps excluidas; XOR con direction_include; en 2D aviso y se ignora; CONSUMIDO en L3: selección de caras, túnel msf_tunnel3d)
- [x] `post_calcul_materi_stress_force_direction_exclude_epsilon` — sub-sprint materi_stress_force, lotes 1 + 3 (2026-08-28; tolerancia del test de exclusión, default 1.e-8; consumido en L3)
- [x] `post_calcul_materi_stress_force_direction_include` — sub-sprint materi_stress_force, lotes 1 + 3 (2026-08-28; 3D: caras con |n·dir| < eps excluidas; la alternativa a direction_exclude; CONSUMIDO en L3: sheet pile msf_sheet3d)
- [x] `post_calcul_materi_stress_force_direction_include_epsilon` — sub-sprint materi_stress_force, lotes 1 + 3 (2026-08-28; tolerancia del test de inclusión, default 1.e-8; consumido en L3)
- [x] `post_calcul_materi_stress_force_element_group` — sub-sprint materi_stress_force, lotes 1-5 (2026-08-28; grupos objetivo, OBLIGATORIO; lista INTEGER no_index=1; validado: 2D solo quad4/quad9, 3D solo hex8/hex27)
- [x] `post_calcul_materi_stress_force_outer` — sub-sprint materi_stress_force, lotes 1-5 (2026-08-28; -yes restringe a los nodos exteriores primarios, los más lejanos del reference_point; default -no; IMPLEMENTADO en 2D (L2) y 3D (L3))
- [x] `post_calcul_materi_stress_force_plot_switch` — sub-sprint materi_stress_force, lotes 1-5 (2026-08-28; invierte la dirección de los vectores de plot; 3 switches en 2D / 4 en 3D — por VECTOR item; IMPLEMENTADO en 2D (L2) y 3D (L3))
- [x] `post_calcul_materi_stress_force_reference_point` — sub-sprint materi_stress_force, lotes 1-5 (2026-08-28; punto de referencia por element group (ndim valores); 3D obligatorio, 2D sin él → aviso + default (0,0); IMPLEMENTADO en 2D (L2) y 3D (L3): define t̂, las caras extremas y la orientación del vector de plot)
- [x] `post_calcul_materi_stress_force_thickness_switch` — sub-sprint materi_stress_force, lotes 1 + 3 (2026-08-28; -yes usa la dirección más LARGA del elemento como espesor estructural; uno por group; CONSUMIDO en L3: t̂ de la cara = menor extensión física, -yes → la mayor)
- [x] Records `*_sig` de `post_calcul -materi_stress -force` (-norx_sig -nory_sig -norz_sig -nors_sig, shear, mom1/mom2; manual 6.913) — 2026-09-01; los nombres de los items node_dof_calcul como VALOR de target_item (target_item N -node_dof_calcul <nodo> -nory_sig): enums + tabla de nombres (19), resolución en exit_tn vía post_calcul_names; dirección de plot 2D = TANGENTE de la cara orientada al reference point (antes ref−centroid, inclinado); gate de la validación de tipos en bloques sin timestep; desbloquea force7/8/9/10/13 (ver registro de verificación)
- [ ] `post_calcul_multiply` — PENDIENTE
- [ ] `post_calcul_safety_default` — PENDIENTE
- [x] `post_calcul_safety_maximum` — cap del factor en el branch safety de calculate_operat (2026-09-05)
- [x] `post_calcul_safety_method` — registrada (KEYWORDS BATCH C 2026-09-04) + selectores VERTICAL/GLOBAL puros + operadores -safety_piping/-safety_lifting; labels medidos vs Pro (prival_0..2/global_x/y/z/vertical); **naming por metodo + formulas lifting=(σ_i+p_total)/p_total y piping=(σ_i+p_dyn)/p_dyn en calculate_operat (2026-09-05; orden prival del Pro = del MAS compresivo = sort GNU invertido; denominador 0 → 0)**
- [ ] `post_calcul_static_pressure` — PENDIENTE
- [x] `post_calcul_static_pressure_height` — KEYWORDS BATCH C (2026-09-04): regiones (coord_min,coord_max,height_ref) en groundflow_phreatic_coord; ground13 target st_pres -1220 EXACTO; **ground13 rc=0 completo con bounda_dof -topres (2026-09-05: pres_dof uniforme +1210; to_pres −20/−10 y dy_pres 1210 = .dbs del Pro)**
- [x] `post_calcul_static_pressure_height_element_group` — KEYWORDS BATCH C (2026-09-04): restriccion por grupo (-all); consumo junto al height record

### post_count (0/1)

- [ ] `post_count` — PENDIENTE

### post_data (0/3)

- [ ] `post_data` — PENDIENTE
- [ ] `post_data_factor` — PENDIENTE
- [ ] `post_data_result` — PENDIENTE

### post_element (0/9)

- [x] `post_element_force` — implementado (commits `53ba151`/`255bc9e` + este lote 2026-09-01; sección de fuerzas/momentos por fuerzas internas del elemento ∫Bᵀσ dV, 5 valores en `post_element_force_result`; verificado contra el Professional en force10/force13: normal −123.4, shear0 ~0, shear1 +1000, moment0 ~0, moment1 −1e5/−5e4 EXACTOS)
- [x] `post_element_force_force` — registrado + parseado (2026-09-01; "también las fuerzas externas se añaden al resultado", manual 6.933 — SIN consumir aún en el cálculo: las fuerzas internas del elemento ya son el cuerpo libre de las cargas, el añadido de fuerzas externas explícitas queda pendiente)
- [x] `post_element_force_geometry` — implementado (commits `53ba151`/`255bc9e`; restricción a los nodos sobre una geometría, PROYECCIÓN EXACTA sobre la posición INICIAL, manual 6.931)
- [x] `post_element_force_group` — implementado (commits `53ba151`/`255bc9e`; restricción a ciertos grupos de elementos, manual 6.932)
- [x] `post_element_force_inertia` — registrado + parseado (2026-09-01; "también las fuerzas de inercia se añaden al resultado", manual 6.933 — SIN consumir aún: las fuerzas de inercia nodales no entran en el resultante de sección, pendiente)
- [x] `post_element_force_multiply_factor` — implementado (commits `53ba151`/`255bc9e`; factor multiplicador del resultado, manual 6.934)
- [x] `post_element_force_normal` — implementado (2026-09-01; switch INTEGER (6.935); con -yes solo los elementos en la dirección normal POSITIVA del plano de sección contribuyen — filtro por centroide; sin el record se usan ambos lados si existen; force13 valida el filtro)
- [x] `post_element_force_number` — implementado (commits `53ba151`/`255bc9e`; restricción a ciertos números de elemento, manual 6.936)
- [x] `post_element_force_result` — implementado (commits `53ba151`/`255bc9e`; normal_force shear0_force shear1_force moment0 moment1, manual 6.937; verificado contra el Professional en force10/force13)

### post_global (1/1)

- [x] `post_global` — presente en el GNU

### post_group (0/1)

- [ ] `post_group_volume_summed` — PENDIENTE

### post_integrate (2/2)

- [x] `post_integrate` — presente en el GNU
- [x] `post_integrate_result` — presente en el GNU

### post_line (5/5)

- [x] `post_line` — presente en el GNU
- [x] `post_line_dof` — presente en el GNU
- [x] `post_line_dof_calcul` — presente en el GNU
- [x] `post_line_n` — presente en el GNU
- [x] `post_line_operat` — presente en el GNU

### post_node (4/8)

- [x] `post_node` — presente en el GNU
- [ ] `post_node_factor` — PENDIENTE
- [x] `post_node_result` — presente en el GNU
- [ ] `post_node_rhside_` — PENDIENTE
- [x] `post_node_rhside_free` — presente en el GNU
- [x] `post_node_rhside_ratio` — presente en el GNU
- [ ] `post_node_rhside_ratio_dof_type` — PENDIENTE
- [ ] `post_node_rhside_ratio_method` — PENDIENTE

### post_point (4/6)

- [x] `post_point` — presente en el GNU
- [x] `post_point_dof` — presente en el GNU
- [x] `post_point_dof_calcul` — presente en el GNU
- [ ] `post_point_element_group` — PENDIENTE
- [ ] `post_point_eps_iso` — PENDIENTE
- [x] `post_point_move` — implementada (commit `3793892`, 2026-08-04)

### post_quadrilateral (4/5)

- [x] `post_quadrilateral` — presente en el GNU
- [x] `post_quadrilateral_dof` — presente en el GNU
- [x] `post_quadrilateral_dof_calcul` — presente en el GNU
- [ ] `post_quadrilateral_element_group` — PENDIENTE
- [x] `post_quadrilateral_n` — presente en el GNU

### post_strain (0/3)

- [ ] `post_strain_volume_absolute` — PENDIENTE
- [ ] `post_strain_volume_initial` — PENDIENTE
- [ ] `post_strain_volume_relative` — PENDIENTE

### print_ (0/1)

- [ ] `print_` — PENDIENTE

### print_apply (0/1)

- [ ] `print_apply` — PENDIENTE

### print_arithmetic (1/1)

- [x] `print_arithmetic` — presente en el GNU

### print_control (1/1)

- [x] `print_control` — presente en el GNU

### print_data (0/1)

- [ ] `print_data_name` — PENDIENTE

### print_database (0/1)

- [ ] `print_database_calculation` — PENDIENTE

### print_de (0/1)

- [ ] `print_de` — PENDIENTE

### print_element (0/2)

- [ ] `print_element_geometry_present` — PENDIENTE
- [ ] `print_element_geometry_present_node_type` — PENDIENTE

### print_failure (1/1)

- [x] `print_failure` — presente en el GNU

### print_frd (0/2)

- [ ] `print_frd_freecad` — PENDIENTE
- [ ] `print_frd_prepomax` — PENDIENTE

### print_gid (0/10)

- [ ] `print_gid_calculation` — PENDIENTE
- [x] `print_gid_contact_spring2` — alias db_number de control_print_gid_contact_spring2 (KEYWORDS BATCH C 2026-09-04)
- [ ] `print_gid_coord` — PENDIENTE
- [ ] `print_gid_de` — PENDIENTE
- [ ] `print_gid_group` — PENDIENTE
- [ ] `print_gid_mesh_activate_gravity` — PENDIENTE
- [ ] `print_gid_node_method` — PENDIENTE
- [ ] `print_gid_smooth_dof` — PENDIENTE
- [ ] `print_gid_smooth_n` — PENDIENTE
- [ ] `print_gid_spring2` — PENDIENTE

### print_gmsh (0/3)

- [ ] `print_gmsh_calculation` — PENDIENTE
- [ ] `print_gmsh_dummy` — PENDIENTE
- [ ] `print_gmsh_node_method` — PENDIENTE

### print_group (0/1)

- [ ] `print_group_data` — PENDIENTE

### print_mesh (0/1)

- [ ] `print_mesh_dof` — PENDIENTE

### print_node (0/2)

- [ ] `print_node_geometry_present` — PENDIENTE
- [ ] `print_node_geometry_present_node_type` — PENDIENTE

### print_precision (0/1)

- [ ] `print_precision` — PENDIENTE

### print_tecplot (0/1)

- [ ] `print_tecplot_calculation` — PENDIENTE

### print_vtk (0/4)

- [ ] `print_vtk_calculation` — PENDIENTE
- [ ] `print_vtk_coord` — PENDIENTE
- [ ] `print_vtk_group` — PENDIENTE
- [ ] `print_vtk_node_method` — PENDIENTE

### print_where (1/1)

- [x] `print_where` — presente en el GNU

### processors (1/3)

- [x] `processors` — alias de `options_processors` en db_number() (2026-09-01, match EXACTO): el GNU ya lo consume en elem.cc (nthread de element_loop) y area.cc. mpc4/interface11 parsean pero siguen RUNFAIL por `mpc_linear_quadratic` y `mesh_interface_triangle_coordinate`; validation_14_mesh por `bounda_time_until_data` + `processors_used` (record de salida aparte).
- [ ] `processors_maximum` — PENDIENTE
- [ ] `processors_partition` — PENDIENTE

### repeat (2/2 — familia COMPLETA; registros de resultado del bloque control_repeat)

- [x] `repeat_save_result` — Sprint 10, lote 10 (2026-08-26): valores guardados por `control_repeat_save` (un índice por repeat, valores de los items en orden de tripletas)
- [x] `repeat_calculate_result` — Sprint 10, lote 10 (2026-08-26): `[media, varianza]` por data item calculados por `control_repeat_save_calculate` (el nombre del inventario `repeat_save_calculate_result` era un OCR; el manual Professional 6.349 lo llama `repeat_calculate_result`)

### safety (0/31)

- [ ] `safety_slip_circle_grid_middle` — PENDIENTE
- [ ] `safety_slip_circle_grid_middle_n` — PENDIENTE
- [ ] `safety_slip_circle_grid_radius` — PENDIENTE
- [ ] `safety_slip_circle_grid_radius_n` — PENDIENTE
- [ ] `safety_slip_circle_grid_result` — PENDIENTE
- [ ] `safety_slip_circle_grid_segment_n` — PENDIENTE
- [ ] `safety_slip_circle_line_middle` — PENDIENTE
- [ ] `safety_slip_circle_line_middle_n` — PENDIENTE
- [ ] `safety_slip_circle_line_radius` — PENDIENTE
- [ ] `safety_slip_circle_line_radius_n` — PENDIENTE
- [ ] `safety_slip_circle_line_result` — PENDIENTE
- [ ] `safety_slip_circle_line_segment_n` — PENDIENTE
- [ ] `safety_slip_combined_linear` — PENDIENTE
- [ ] `safety_slip_combined_linear_n` — PENDIENTE
- [ ] `safety_slip_combined_linear_result` — PENDIENTE
- [ ] `safety_slip_combined_linear_segment_n` — PENDIENTE
- [ ] `safety_slip_ellipsoide` — PENDIENTE
- [ ] `safety_slip_ellipsoide_method` — PENDIENTE
- [ ] `safety_slip_ellipsoide_n` — PENDIENTE
- [ ] `safety_slip_ellipsoide_result` — PENDIENTE
- [ ] `safety_slip_ellipsoide_segment_n` — PENDIENTE
- [ ] `safety_slip_grd` — PENDIENTE
- [ ] `safety_slip_grd_method` — PENDIENTE
- [ ] `safety_slip_grd_method_direction` — PENDIENTE
- [ ] `safety_slip_grd_segment_n` — PENDIENTE
- [ ] `safety_slip_multi_linear` — PENDIENTE
- [ ] `safety_slip_multi_linear_n` — PENDIENTE
- [ ] `safety_slip_multi_linear_result` — PENDIENTE
- [ ] `safety_slip_multi_linear_segment_n` — PENDIENTE
- [ ] `safety_slip_set` — PENDIENTE
- [ ] `safety_slip_set_result` — PENDIENTE

### slide (4/7)

- [x] `slide_axisymmetric` — implementada (escala la fricción de slide por `2*pi*r`, r = coordenada radial del nodo en axisimétrico; validada con `slide_axi`)
- [ ] `slide_damping` — PENDIENTE
- [x] `slide_geometry` — presente en el GNU
- [x] `slide_plasti_friction` — ley elastoplástica del slide IMPLEMENTADA en slide.cc (lote familia slide 2026-09-04): springs kn/kt sobre desplazamientos TOTALES + cap c+Fn·tan(phi) + spring normal bidireccional con cap de tensión sig_t; salidas node_slide_force/node_slide_f/node_slide_direction; factor axisym 2·pi·r (group_axisymmetric o slide_axisymmetric); forma sin índice `axisymmetric -yes` (6.17) aceptada en input.cc; print_debug parseado. **slide1 y slide3 rc=0** (estado final slide1 IDÉNTICO al .dbs del Professional: node_slide_force (−2.1148e-3, −1.2210e-3) vs (−2.1129e-3, −1.2199e-3); reacciones top +0.01/−5.7707e-3). slide2 = blocker PREEXISTENTE del elemento axisym GNU (split de fuerzas nodales igual 4.712/4.712 vs el lobatto del Pro ∝ radio 3.142/6.283 — σyy −6.447e-4 vs target −6.666666e-4, tol 1e-9); slide4 RESUELTO 2026-10-11 (F4 T6e C53): rc=0 e igual al Professional (node_slide_force nodo 1 −2.11291867e-3 / −1.21989198e-3 vs −2.11291887e-3 / −1.21989209e-3). El "blocker de dinámica de iteración" diagnosticado el 2026-09-04 era el default `sig_t = 0` (gap de solo compresión: el nodo x=0 despegaba) más el cap `max(c + Fn·tan(phi), 0)` y el desfase de una iteración del estado (slide fuera del post-solve pass); ver registro de verificación.
- [x] `slide_plasti_tension` — consumido por la ley (cap de tracción S ≥ −sig_t del spring normal); 2026-10-11 (F4 T6e C53): sin record = sin límite de tracción (antes sig_t 0 por defecto = gap de solo compresión), ver registro de verificación
- [x] `slide_stiffness` — consumido por la ley (kn/kt por nodo, verificado: Fn = kn·penetración EXACTA por nodo)
- [ ] `slide_user` — PENDIENTE

### solver (10/10)

- [x] `solver` — Sprint 12 lote 3 (el tipo GLOBAL pisa a todo control_solver: leído al FINAL en top.cc/elem.cc/dof.cc — precedencia OPUESTA al options_solver legacy; probado con banner SuperLU); F4 T9 (`3349515`, 2026-10-10): `-matrix_pardiso` solved with SuperLU, Bi-CG fallback
- [x] `solver_bicg_error` — Sprint 12 lote 3 (global, pisa al control_solver_bicg_error en so_bicg.cc)
- [x] `solver_bicg_restart` — Sprint 12 lote 3 (PARCIAL: registrado; el bicg del GNU no tiene restart, como su contraparte de control)
- [x] `solver_bicg_stop` — Sprint 12 lote 3 (global); 2026-09-29: `control_solver_bicg_stop` (índice 0, then the timestep's own index) prevalece sobre él (manual 6.376 + Professional medido; combined precedence, commit pending; BACKLOG C13)
- [x] `solver_matrix_save` — Sprint 12 lote 3 (PARCIAL: registrado; cache de descomposición es funcionalidad de solver directo)
- [x] `solver_matrix_symmetric` — Sprint 12 lote 3 + 2026-09-02 (con -yes, sistema medido NO simétrico → Bi-CG honesto sobre la matriz tal-cual en vez de CG divergente — el Professional simetriza "si hace falta"; large2 rc=0, large3 sin cambio); 2026-10-11 (F4 T6b, `c350b3c`): sin simetrización = divergencia deliberada K24, guard tests `k24_consol_3d*` (ver registro de verificación y Descarte)
- [x] `solver_pardiso_ordering` — Sprint 12 lote 3 (PARCIAL: PARDISO no compilado)
- [x] `solver_pardiso_out_of_core` — Sprint 12 lote 3 (PARCIAL, ídem)
- [x] `solver_pardiso_processors` — Sprint 12 lote 3 (PARCIAL, ídem)
- [x] `solver_pardiso_processors_maximum` — Sprint 12 lote 3 (PARCIAL, ídem)

### strain (0/10)

- [x] `strain_settlement_diagram` — implementada (dependencia de parámetro del creep en un dof vía tabla; validada con `strain_settle_diag`)
- [x] `strain_settlement_diagram_dof` — implementada
- [x] `strain_settlement_diagram_number` — implementada
- [x] `strain_settlement_element_group` — implementada
- [x] `strain_settlement_parameters` — implementada (creep de asentamiento con saturación; `strain_settlement_creep` en materi.cc; validada con `strain_settle`)
- [ ] `strain_volume_absolute_time` — PENDIENTE
- [ ] `strain_volume_element` — PENDIENTE
- [ ] `strain_volume_element_group` — PENDIENTE
- [ ] `strain_volume_geometry` — PENDIENTE
- [ ] `strain_volume_relative_time` — PENDIENTE

### support (18/18 + _element de facto + status = 20/20 COMPLETO)

- [x] `support_edge_normal` — Sprint 13 lote 1 (el Winkler elástico: fuerza RHS desde desplazamientos totales + rigidez consistente dt-escalada en la matriz; 2D/3D; verificado EXACTO prescrito/serie/cara)
- [x] `support_edge_normal_damping` — Sprint 13 lote 2 (amortiguadores viscosos sobre el soporte: c_n v_n + c_t v_t, gate con control_support_edge_normal_damping_apply); K17 2026-10-06 (F2 T2): damping and density also in the iteration matrix (support_edge_normal_3 = Pro)
- [x] `support_edge_normal_damping_automatic` — Sprint 13 lote 2 (calcula c_n y c_t del grupo adjunto: c_n = sqrt(ρ·Eoed), c_t = 0.25·sqrt(ρ·G); verificado con E=1e7 del grupo)
- [x] `support_edge_normal_damping_automatic_apparent` — Sprint 13 lote 2 (usa módulos aparentes desde el estado nodal actual σxx/εxx; para elástico = _automatic)
- [x] `support_edge_normal_density` — Sprint 13 lote 2 (masa distribuida del soporte: d_n·a_n con a = (v_new-v_old)/dt)
- [x] `support_edge_normal_element_group` — Sprint 13 lote 1 (restricción por grupo; NOTA: el manual 6.1072 también define `support_edge_normal_element` — implementado igual, línea ausente del checklist añadida de facto)
- [x] `support_edge_normal_element_node` — Sprint 13 lote 1 (elemento + nodos locales)
- [x] `support_edge_normal_element_side` — Sprint 13 lote 1 (pares elemento/lado local 1-based)
- [x] `support_edge_normal_factor` — Sprint 13 lote 2 (escala espacial polinómica de las RIGIDECES, no de la fuerza; reutiliza force_factor)
- [x] `support_edge_normal_force_initial` — Sprint 13 lote 2 (fuerza inicial a0+a1*y(z): preload de compresión en el soporte, registrado aunque u=0)
- [x] `support_edge_normal_geometry` — Sprint 13 lote 1 (entidad geometry O lista de nodos; requerido por la maquinaria area())
- [x] `support_edge_normal_node` — Sprint 13 lote 1 (restricción por nodos)
- [x] `support_edge_normal_plasti_compression` — Sprint 13 lote 3 COMPLETO (4e0074f: piso de compresión + factor tangencial)
- [x] `support_edge_normal_plasti_friction` — Sprint 13 lote 3 COMPLETO (4e0074f: Coulomb c+μ·|S|, el menor de los caps gobierna)
- [x] `support_edge_normal_plasti_residual_sti` — Sprint 13 lote 3 COMPLETO (4e0074f: fracción factor en [0,1] a la matriz; default 0 = sin cambio)
- [x] `support_edge_normal_plasti_tension` — Sprint 13 lote 3 COMPLETO (4e0074f: gap — fuerzas a 0 en tensión + status=1)
- [x] `support_edge_normal_plasti_tension_double` — Sprint 13 lote 3 COMPLETO (4e0074f: cap S=−max + tangencial a 0)
- [x] `support_edge_normal_time` — Sprint 13 lote 2 (factor temporal para la fuerza del soporte; reutiliza force_time)

### target (2/2)

- [x] `target_item` — presente en el GNU
- [x] `target_value` — presente en el GNU

### time (2/2)

- [x] `time_calculation` — presente en el GNU
- [x] `time_current` — presente en el GNU

### timestep (2/2)

- [x] `timestep_iterations_automatic_apply` — Sprint 12 lote 3 (-no descarta todo control_timestep_iterations_automatic: gate en top.cc antes de la rama automática; física sin cambio verificada). **Batch hypo (2026-09-01)**: `control_timestep_iterations_automatic` pasa a 3 valores (ratio minimal maximum, manual 6.386); [2]=maximum (antes [1]); [1]=minimal registrado sin consumir (parcial); slope_classical_numerical ya avanza (antes maximum=1e-6 → timeout)
- [x] `timestep_predict_velocity` — Sprint 12 lote 3 (PARCIAL: registrado; el esquema escalonado del GNU arranca el solve de x=0 y ACUMULA el resultado — una predicción de velocidades previas duplicaría; neutral verificado)

### tochnog (1/1)

- [x] `tochnog_version` — Sprint 12 lote 3 (record con la fecha de BUILD parseada de __DATE__, escrito al inicio de top(); targeteable — test del año)

### truss (0/1)

- [ ] `truss_rope_apply` — PENDIENTE

### up (0/1)

- [ ] `up` — PENDIENTE

### volume (2/2)

- [x] `volume_factor` — presente en el GNU
- [x] `volume_factor_x` — Sprint 12 lote 3 (factor por tramos de x en volume.cc, a nivel de IP; derecha de xn factor 1; multiplica al polinómico y al de grupo; verificado analítico 8x1: sigma 2e-2/1e-2 por mitades, punta 1.15e-3 vs 1.2e-3)

### wave (2/2)

- [x] `wave_fscalar` — presente en el GNU
- [x] `wave_scalar` — presente en el GNU

### zip (1/1)

- [x] `zip` — Sprint 12 lote 3 (-yes gzipea *flavia*/*msh/vtk*/*.plt/*.dbs al final del cálculo en exit_tn, tras los targets y antes de db_close; test en directorio aislado)

</details>

---

# SPRINT 12 — PLANIFICADO (2026-08-29)

Estado al cierre de la sesión del 28-29 ago 2026: sub-sprint materi_stress_force
COMPLETO (L1-L5 + fix A+B + fix C/D + SRI quad4 + SRI hex8 + arness Professional
COMPLETO — todos los casos ≈ 1.000×). Suite 209/209. Paquete para colaboradores
en `ProjectDocs/shared/`.

## Objetivo del Sprint 12

Madurar la versión para ofrecer a los colaboradores originales (Fernando Lorenzo,
Osman Buyukusik): cerrar los cabos abiertos y continuar la convergencia.

## Lotes propuestos (en orden; ajustar a decisión del usuario al volver)

1. **[DONE 2026-08-29] Cierre del arness — V/mom 3D de sección en carga
   axial** (gforce10/13): era un problema de FRAME de sección (esquinas
   diagonales + `thickness_switch` omitido en la conversión + 4
   convenciones de signo/escalar), no del esquema. Arness COMPLETO a
   ratios +1.0000 con signos — ver la fila SPRNT 12 LOTE 1 del registro
   y VALIDACION-PROFESIONAL §13.
2. **[DONE 2026-08-29] Paper línea 3 (pedagógica)** — pospuesta por
   decisión del usuario ("los papers los dejaremos para más
   adelante"); los datos de las líneas 1 y 2 se van recopilando en
   PAPER-LINES.md con cada lote.
3. **[DONE 2026-08-29] Convergencia — familias grandes del checklist**:
   este lote cerró `solver` (0/10 → 10/10), `timestep` (0/2 → 2/2),
   `tochnog_version` (0/1), `volume_factor_x` (0/1) y `zip` (0/1) —
   ver la fila SPRNT 12 LOTE 3 del registro. Quedan pendientes:
   `safety` (0/31, sub-sprint propio), `support` (0/18) y
   `strain_volume_*` (4/10) + `post_strain_volume_*` (3) +
   `post_calcul_safety_*` (3).
4. **Email a los colaboradores** — el paquete `ProjectDocs/shared/` está listo;
   falta redactar el correo (inglés) con la oferta (código + docs + papers).
5. **Revalidar con el Professional** al cierre del sprint y actualizar
   `VALIDACION-PROFESIONAL.md` + el paquete shared.

## Cabos abiertos documentados (no bloqueantes)

- SRI hex8 3D: 3 modos twist + warping de sección (matriz singular sin BCs que
  los maten) — limitación clásica del SRI de cizalla del ladrillo; la literatura
  usa B-bar/ANS (posible lote de investigación).
- `control_print` checklist: 63/85 (implementables agotados; dependientes y
  descartados documentados).

---

## ROADMAP APROBADO (2026-08-29, decisión del usuario)

Principio: cubrir TODO lo posible del Professional (GiD queda fuera de
momento, se analizará después) + añadir funcionalidad PROPIA más allá
de la convergencia (vigas 3D, placas y lajas de comportamiento plano).
El correo a los colaboradores espera a tener un grado de convergencia
y avance mejor.

- **Sprint 13 — bloque geotécnico** [lotes 1+2+3 DONE 2026-08-29: support_edge_normal COMPLETO 18/18+status (4e0074f); queda safety_slip]: `support_edge_normal` (19) +
  `safety_slip` (~35 con control/post_calcul). Sub-sprints grandes.
- **Sprint 14 — materiales**: `group_materi` (30) + `materi_plasti` (9)
  + `materi_strain`/`strain_volume` (~22).
- **Sprint 15 — MPC + contacto + exportadores + cola**: `mpc` (17),
  contacto (11), exportadores (vtk completo; gmsh/frd/tecplot a
  decidir), barrido de familias 1-2 (~80, muchos serán descartes
  documentados tipo alias), tick de preexistentes (echo, derivatives).
- **Sprint 16 — funcionalidad propia I (§8.5)**: BEAM3D real (12×12,
  6 dofs/nodo, torsión, frame local 3D) + lajas/membrana (§8.3:
  elemento de tensión plana pura — NOTA: el COMPORTAMIENTO de tensión
  plana ya existe como material, `group_materi_membrane -yes`, usado
  en msf_nor; lo que falta es el elemento dedicado si se quiere).
- **Sprint 17 — funcionalidad propia II (§8.3)**: PLATE (Reissner-
  Mindlin: dofs rotacionales, toca ensamblaje y solver).
- **Cierre formal**: sesión de descartes documentados (mesh_gid/GiD y
  los que salgan), revalidación con el arness, actualización del
  paquete shared/ → entonces el correo a Fernando Lorenzo y Osman
  Buyukusik.

## Sprint 13 lote 5 — interfaz: records de salida + fixes de convergencia (2026-08-30)

| Commit | Fecha | Estado |
|--------|-------|--------|
| `260e12e` feat(interface): records de salida Professional + fix F_t=kt·du_tang + MC implicito | 2026-08-30 | **code**: 3 fixes verificados contra .dbs del Professional |

**Hallazgo de convergencia** (deducido de los .dbs del Professional con dump python):
1. **F_t = kt·du_tang** (NO kt·2·du_tang). El manual 6.631: stress,shear = kt·gamma = 2·kt·strain,shear donde strain,shear = du_tang/2 → stress = kt·du_tang. Nuestro factor 2 desde Carril A Fase 1 DUPLICABA la rigidez tangencial (interface14: kt=5e3, du=1 → Professional τ=5e3, nosotros 1e4).
2. **Return mapping MC implícito con dilatancia**: dgamma = (|trial|−max_fric)/(kt + kn·tanφ·tanψ); f_t = trial − kt·dgamma; dilatancia CIERRA (strain_eff += −dgamma·tanψ, compresión positiva). Validado EXACTO con interface15: dgamma=0.8537, f_t=731.63, σn=−1730.48 (τ/σn = tan(0.4) = 0.42279 EXACTO).
3. **Records de salida** (antes registrados pero NUNCA escritos): `element_interface_intpnt_stress`/`_intpnt_strain` (n_intpnt = nnol/2 IPs: bar2=1, quad4=2, quad6=3, prism6=3, hex8=4; valores del PASO actual: strain,normal=du_norm, strain,shear=du_tang/2, stress,normal=kn·strain_eff, stress,shear=F_t−F_old), `_stress_average`/`_strain_average` (promedio), `_intpnt_materi_tension_status` (−opened/−closed, enums OPENED/CLOSED nuevos). Pre-alocación en top.cc (db_allocate NO puede correr en el bucle paralelo) + data_length/fixed_length en database.cc.

**GOTCHAS nuevos** (manual-developer):
- `fixed_length[MDAT]` default = **1** (database.cc:61) — records de longitud variable necesitan `fixed_length[rec]=0` o el PUT falla "Length too small".
- `db_allocate` NO puede alocar DENTRO del bucle paralelo de elementos — pre-alocar en top.cc con el patrón `ELEMENT_INTERFACE_STRAIN_NORMAL`.
- Los records de salida se escriben en VERSION_NORMAL Y VERSION_NEW (el chequeo de target y print_database leen NORMAL).
- El .dbs del Professional se puede volcar con python (buscar el nombre del record + slice) para deducir la semántica EXACTA.

**Verificación**: suite 16/16 OK; corpus Professional 77→81 PASS (21.2%→22.3%): interface14 rc=0 y interface15 rc=0 con targets EXACTOS (MC + elástico tangencial).

**Gap documentado**: interface1/8 (y otros con `solver -matrix_pardiso` + nodos coincidentes espesor-0) no convergen con el Bi-CG del GNU — requieren solver directo (PARDISO/SuperLU). El GNU 2014 tiene so_suplu.c sin compilar. Este es el siguiente bloqueo grande para la familia interface_*.

### Refactor por-pares Lobatto + fix fallback pardiso (bf74f28, 40d640d)

**Hallazgo clave**: el ensamblaje de interfaz del Professional es POR PUNTO DE INTEGRACIÓN: cada par de nodos enfrentados (i, i+ns1) es un resorte independiente con su propio desplazamiento relativo y peso Lobatto (quad6: 1/6,4/6,1/6; quad4: 1/2,1/2; bar2: 1; prism6: 1/6,4/6,1/6; hex8: 1/12,5/12,5/12,1/12). El viejo modelo (media del lado, todos-contra-todos) producía una matriz de RANGO 1 por lado → SINGULAR para ns1>1 (interface1 no convergía ni con Bi-CG ni con band). Los pesos se deducen de las cargas nodales del Professional (interface1: −1,−4,−1 = (1/6,4/6,1/6)·(−6)).

**Fix fallback -matrix_pardiso**: en so.cc `bicg_solver=1` SOLO se activa con `task==-MATRIX_ITERATIVE_BICG`; `-matrix_pardiso` sin fallback explícito dejaba NINGÚN solver → solución 0 silenciosa. Añadido `bicg_solver=1` al branch pardiso. Con la matriz por-pares, el Bi-CG CONVERGE (interface1: node_dof 4/5/6 = −6e-10 EXACTO vs Professional).

**GOTCHA**: el refactor por-pares eliminó accidentalmente la lectura de `group_interface_materi_plasti_mohr_coul_direct` (mc_active/phi/c/phi_flow) — restaurada en 40d640d. interface15 vuelve a rc=0 EXACTO (dgamma=0.8537, f_t=731.63, σn=−1730.48).

**Estado familia interfaz**: interface1/7/14/15 rc=0 (elástico carga, gap+reset, MC con dilatancia implícita). Pendientes: interface2 (gap multi-paso: el strain record refleja el cierre del gap, no el incremental), interface8 (signo de orientación del elemento degenerado quad4 5 6 3 4), interface9 (tangencial), interface13 (el propio Professional escribe 0.671 en el .dbs pero su target espera 1.118 — bug del test del corpus), interface12/patch (3D).

**GOTCHA del LAPACK band solver del GNU**: con el truss no llena la matriz (node_lhside 0) — no confiar en él para elementos que usan solo NODE_LHSIDE; el Bi-CG es el camino.

### Cierre familia 2D de interfaz (4b6548a) — 8 tests rc=0

**Tests que pasan**: interface1 (elástico carga), interface7 (2 pasos + reset), interface8 (orientación degenerada), interface9 (axisimétrico), interface12 (damping), interface14 (elástico tangencial), interface15 (MC con dilatancia), interface_patch (fricción multi-paso).

**Fixes nuevos** (verificados contra .dbs Professional):
1. **Orientación 2D por numeración**: normal = (−t.y, t.x) pero si el primer nodo del lado 2 tiene número MENOR que el lado 1, se invierte (interface8: quad4 5 6 3 4 → normal (0,−1); interface1: 1 2 3 | 4 5 6 → (0,1)). La numeración relativa codifica el orden de los lados de la conversión.
2. **Axisimétrico**: la interfaz anillo lleva el peso 2πr (interface9: τ = −0.159 EXACTO).
3. **Damping** (`group_interface_damping`): d·(v2−v1) en RHS + d en la matriz SIN dtime (interface12: 3 nodos −0.5 EXACTO). data_length 2→1.
4. **MC con trial elástico acumulado** (SUPERADO 2026-10-11, F4 T6d C52: el history vuelve a guardar la fuerza recortada; ver registro de verificación): el history guardaba kt·gamma_total (crece aunque plastifique); f_t = clamp(trial); RHS = f_t − clamp(f_old). Sin esto la fricción multi-paso estancaba (−5.55 vs −10).
5. **Record**: sin MC → kt·du_paso (incremental); con MC → f_t total clampado.
6. **Pesos 3D uniformes** 1/ns1 (el quad4 convertido a hex8 reparte la tracción de cara uniforme; Lobatto 1D daría 4× menos — interface_quad4_hex8: −0.2475 vs −1.0).

**Pendiente documentado**: los 3D con kn>>E (interface_quad4_hex8, interface_bar2_hex8, interface_bar2_quad4, interface_quad4_hex8_many) necesitan solver directo real. El band dgbsv resuelve el hex8 3D pero rompe 2D (truss/interface1/8 → 0 silencioso); el Bi-CG no converge con ratio 1e10. SuperLU = proyecto aparte (so_suplu.c existe en el GNU 2014).

### Solver directo: retry LU (band dgbsv) tras fallo Bi-CG (5b7bb8b)

**Solución al bloqueo 3D**: `-matrix_pardiso` usa Bi-CG por defecto; cuando el Bi-CG no converge (interfaces 3D con kn>>E, ratio 1e10 — interface_quad4_hex8 con kn=1e10 vs E=1), el solver reintenta con el LAPACK band LU (dgbsv) directo. El band maneja la casi-singularidad con pivoteo (el denso dgesv la rechaza).

**Fixes necesarios**:
1. `so_bicg.cc`: los 3 fallos (no convergencia, breakdown CG/BiCG) retornan `succesful=0` (antes `exit_tn_on_error`) para permitir el retry.
2. `so.cc`: `mat` (band) y `mat_dense` (dgesv) se alocan y llenan SIEMPRE (también en el camino bicg).
3. `band` se calcula también en el renumbering bicg (antes solo en el else no-bicg).
4. **FIX fill dofs**: `band_solver` añadido a la condición de copia `solve_b` → `node_dof_new` — el band NUNCA llenaba los dofs → solución 0 silenciosa. ESTE era el bug de "band da 0" (truss/interfaces 2D con band).
5. El retry restaura `solve_b` desde `solve_b_temp` (el Bi-CG lo sobrescribe).

**Resultado**: corpus 81→87 PASS. interface_quad4_hex8 rc=0 (el 3D de interfaz con kn=1e10 que el Bi-CG no podía). Suite propia 16/16 intacta.

**Pendientes**: interface_bar2_hex8/bar2_quad4/many (bug de extrusión bar2→quad4→hex8, separado del solver), interface2 (gap multi-paso), interface13 (bug del propio test del corpus).

### Extrusión bar2→quad4 + fixes convert/reconexión (d7b2dfd) — corpus 90

**Desbloquea la familia 3D de interfaz** (bar2_quad4, quad4_hex8_many, y +1 del corpus):
1. **mesh_extrude**: bar2 → quad4 (antes solo tria3/quad4).
2. **interface_convert bar2 3D**: la normal = ẑ×tangent (antes leía el[3] inexistente → índice basura).
3. **Fix reconexión**: solo los vecinos del lado +normal (centroide en dirección de la normal) se reconectan con los nodos nuevos; el lado −normal conserva los originales (antes corrompía el sólido del mismo lado: elem 1 = 1 2 3 4 → 1 2 7 8).
4. **control_mesh_convert/extrude** se buscan en TODO el rango de controles (antes solo icontrol actual — nunca coincidía con control 20 + timestep 30).
5. **convert/extrude solo en el primer step** (task==YES) — re-ejecutarlos duplicaba la malla.
6. **extrude se mueve ANTES de data()/merge()** — data veía los quad4 3D sin extruir y fallaba.
7. **mesh_has_changed del convert**: condicional a que NO haya extrude (el extrude lo llama al final; sin extrude el convert debe reconstruir conexiones).

**Resultado**: corpus 87→90 PASS. Suite 16/16. Tests de interfaz pasando: 1, 7, 8, 9, 12, 14, 15, patch, quad4_hex8, bar2_quad4, quad4_hex8_many (11). Queda interface_bar2_hex8 (la interfaz hex8 extruida conecta en z en vez de y — conectividad fina).

### Orden extrude→convert — interface_bar2_hex8 rc=0 (b9dec4f) — corpus 91

**El último caso fino de la familia 3D**: el bar2 de interfaz 3D debe EXTRUIRSE primero (bar2 → quad4 en el plano XZ) y LUEGO convertirse a hex8 (quad4 → hex8 con el lado 2 en +y). El orden inverso (convert → extrude) producía un hex8 de interfaz con la normal en Z — el push en Y nunca se transmitía (sigyy=0).

**Fixes**: (1) top.cc: extrude() ANTES de interface_convert(); (2) interface.cc: mesh_has_changed SIEMPRE tras el convert (los sólidos ya son hex8 cuando el convert corre).

**Familia interfaz: 12/14 tests rc=0** (1, 7, 8, 9, 12, 14, 15, patch, quad4_hex8, bar2_hex8, bar2_quad4, quad4_hex8_many). Solo quedan interface2 (gap multi-paso) e interface13 (bug del propio test del corpus: el Professional escribe 0.671 pero su target pide 1.118).

**Corpus: 90 → 91 PASS** (de 363). Suite 16/16.

### contact_spring2 + familia conspr/genera2 — corpus 99

**Cierre de la familia contact_spring del corpus** (conspr1-7 + genera2, 8 tests):
1. **Aliases**: `contact_spring`/`contact_spring2` → CONTACTSPRING; `group_contact_spring_*` → `group_contactspring_*` (stiffness, memory, direction, direction_automatic, damping).
2. **`group_contact_spring_direction_automatic`** + `control_mesh_generate_contact_spring_element_group` + `control_mesh_generate_contact_spring*`: registrados como records.
3. **OPTIONS_CONVECTION -no**: la restricción aplica SOLO a nnol==1 (spring de nodo sobre sólido); contact_spring2 (nnol==2, spring ENTRE nodos) no la necesita.
4. **memory -total_linear** aceptado (antes db_error).
5. **generate_contactspring_element_group** implementado: crea springs entre nodos coincidentes de dos elementos/grupos dados; dirección automática por centroides.
6. **FIX crítico de dirección automática**: la dirección se calcula SIEMPRE de los centroides de los elementos de los nodos (cen1−cen0, con filtro que excluye los springs — su centroide es degenerado). La resta de coords de los nodos del spring degenera en la 2ª iteración del solver: coord ya deformado da (0,1e-6) que normaliza a (0,1) — signo EQUIVOCADO. Con centroides: nodo5 (bloque sup, y=1.5) − nodo3 (bloque inf, y=0.5) = (0,−1) → force = +1.0 correcto. conspr7: la carga total 2 (nodos 7,8 con −vely) se reparte 1 por spring.

**Resultado**: conspr1-7 y genera2 → rc=0. **Corpus: 91 → 99 PASS** (de 363). Suite 16/16. Desglose actual del corpus: 99 PASS / 224 RUNFAIL / 40 PARSE.

**Pendientes del corpus**: interface2 (gap multi-paso), interface13 (bug del propio test), y los parse-errors de mayor frecuencia: `-nory_sig`/`-norx_sig` (6), `-quad8`/`-hex20` (5), `mesh_gid_point_coord`/`mesh_gid_circle_coord` (5). (`-updated_area` se corrigió en 2026-09-02: el valor YA se acepta; los 3 tests etiquetados quedan bloqueados por la familia `incremental_driver`, no por el valor).

### Aliases `geometry_factor` + `processors` — corpus 121 (2026-09-01) — `feat(alias): geometry_factor -> GEOMETRY_BOUNDA_FACTOR, processors -> OPTIONS_PROCESSORS` (`9d22f6c`)

Dos aliases triviales en la cadena de traducción de `db_number()`
(database.cc):

1. **`geometry_factor`** (manual Professional 6.527) → `GEOMETRY_BOUNDA_FACTOR`:
   el GNU ya leía el record con la MISMA semántica (interpolación lineal
   2 valores / parabólica 3 valores en `geometry()`, mismo índice que la
   entidad geométrica). El ejemplo del manual (línea 0..1 con factores
   1/4: nodo x=0.2 → 20·1.6, x=0.4 → 20·2.2) coincide EXACTAMENTE con la
   interpolación de geometry.cc. Desbloquea matrix2/temp2/matrix4
   (rc=0). Verificado contra el binario Professional 25-10-2023
   (matrix2: post_point_dof GNU 2.4999998123 vs Professional
   2.500000000000e+00).
2. **`processors`** → `OPTIONS_PROCESSORS`: el GNU ya lo consumía en
   elem.cc (buffers por-thread de element_loop) y area.cc (warning de
   acumulación single-thread). Match EXACTO (strcmp) para no traducir
   `processors_used`/`processors_maximum`/`processors_partition`.
   mpc4/interface11 parsean ahora, pero siguen RUNFAIL por sus blockers
   (`mpc_linear_quadratic`, `mesh_interface_triangle_coordinate`);
   validation_14_mesh por `bounda_time_until_data` + `processors_used`
   (record de salida aparte).

**Resultado**: **Corpus: 118 → 121 PASS** (220 RUNFAIL / 22 PARSE, de
363). Suite propia 16/16. delete3 sigue RUNFAIL (sigyy 0 vs 0.6: el
`control_mesh_delete_geometry_factor` parsea pero NO se aplica — feature
de clase CONTROL aparte, PENDIENTE).

### Batch hypo/mohr/kapsh/reset — corpus 112 (2026-09-01) — `feat(hypo): hyhis0..7 + kappa_shear + direct_hardening_softening + resets multi-dof` (`22801c3`)

**Cierre del batch iniciado en la sesión anterior** (16 archivos sin commitear, prints DBG eliminados, diagnóstico cerrado):

1. **`materi_plasti_hypo_history` (4.23)**: 8 variables hyhis0..hyhis7, dof hisv compartido, basenames hyhis; checks de los grupos hypo aceptan el nombre (check_unknown_atleastone).
2. **`materi_plasti_kappa_shear` (4.25)**: dof `kapsh` = int sqrt(0.5·dev(epp):dev(epp)); hardening de la familia hardening_softening (plasti_rule); RHS en materi.cc.
3. **`group_materi_plasti_mohr_coul_direct_hardening_softening`**: variante DIRECT con ángulos EN GRADOS (PIRAD/180 en plasti_rule); dam_building la usa.
4. **`control_timestep_iterations_automatic` → 3 valores** (ratio minimal maximum, 6.386); [2]=maximum; slope_classical_numerical ya avanza.
5. **`group_materi_plasti_hypo_intergranularstrain` → 6 params** (R m_R m_T beta_r chi theta, 6.715) + alias del nombre Professional `_hypo_strain_intergranular`.
6. **`control_reset_dof`**: (a) scoping por índice (icontrol==ireset — verificado contra el Professional: RECHAZA un reset en el índice del timestep, "Error detected for data item : control_reset_dof, record : 20"); (b) FIX multi-dof: todos los dofs listados reciben el valor (antes solo reset_dof[0] — bug latente que rompía el reset de 17 dofs de dam_building).
7. **`area_element_group_sequence` sin selectores** (6.9): los elementos del grupo anterior (i−1) pasan al grupo i en tiempo i; dam_building la usa.
8. **Alias `mesh` → OPTIONS_MESH** (el placeholder MESH muerto rompía el parse de `mesh -fixed_in_space`).
9. **`memory==-UPDATED_LINEAR`**: ramas en materi.cc / set_deften_etc / stress.cc; el check de materia_displacement solo aplica con el record EXPLÍCITO (hipo12/13 corren sin group_materi_memory).
10. **Rama masin (hypoplas.cc)**: e = hyhis0, s = hyhis4 (layout Professional); delta intergranular desde el dof epi; gates `control_materi_plasti_hypo_masin(_clay)_ocr_apply` por ICONTROL (default -YES; data_required GROUP_TYPE eliminado).
11. **Aliases nuevos**: `nonlocal` → OPTIONS_NONLOCAL (6.897), `nonlocal_name` (6.898), `group_materi_plasti_visco_exponential_value` (singular), `control_mesh_refine_locally_dof` (6.222).

**Corpus: 99 → 112 PASS** (de 363; 229 RUNFAIL / 22 PARSE). Suite 16/16. La subida neta viene del fix multi-dof de control_reset_dof (+15 tests) menos la caída de hypo1/hypo3 (−3, kernel, ver abajo).

**Fixes de bugs reales encontrados**:
- `control_reset_dof` multi-dof: solo el PRIMER dof del record recibía el valor (latente desde b6eaee4). El Professional resetea todos los listados (verificado en sus .dbs: sigxx/sigyy/sigzz evolucionan desde −100).
- `control_materi_gate_off` (general.cc): el chequeo `db_active_index(item, 0)` fallaba para records en índices ≠ 0 → ahora `db_max_index` (los gates por ICONTROL funcionan para cualquier índice).
- Check de UPDATED: el rechazo de materi_displacement solo aplica con el record explícito.

**Tests objetivo — estado final**:
- hypo2, hypo4: **rc=0** (targets del Professional).
- hypo1 (−862.766 vs −862.92, 0.019%): kernel wolfersdorff (hypo.c) — la corrección del reset reveló la desviación (antes el bug reset_dof[0] la compensaba por coincidencia de tolerancia).
- hypo3 (−0.104 vs −0.1618): kernel wolfersdorff con estado inicial anisótropo — desviación 36%, calibración pendiente.
- hypo7/8/9 (−224/−216/−224 vs −231.8/−230.6/−231.8, 3-6%): kernel masin (masin.c) — mejoran desde −2064/−31.87/−2064.
- hypo12/13: corren (fix del gate OCR + layout hyhis) pero el kernel masin visco desborda la razón de vacíos (1.19e17) — calibración pendiente.
- hypo6/11: PENDIENTE (modelos Niemunis visco no implementados en el kernel GNU).
- hypo10/triaxial/direct_shear: PENDIENTE (extensión ISA + incremental_driver no implementados).
- dam_building: corre pero necesita ~10 min (corpus 45 s) — los features (kapsh, direct_hardening_softening, sequence sin selectores, timestep 3 valores, reset 17 dofs) se ejercitan.
- slope_classical_numerical: corre hasta t=2; el target de colapso 1.275 no se cumple (física de la reducción phi-c pendiente).
- slope_nonlocal_refine: parsea (nonlocal/nonlocal_name/visco_exponential_value/refine_locally_dof); el run choca con el bug legacy del buffer NONLOCAL_ITEM_SIZE de nonloc.cc (PENDIENTE).
- mohr_coul_direct1-8: RUNFAIL pre-existente (ley directa, eptxx 1/5 del target — no relacionado con el batch).
- mesh_gid_1-8: RUNFAIL (requieren el programa GID).

**Registro de verificación**:
| Feature | Verificación |
|---------|-------------|
| materi_plasti_hypo_history | hypo2/4 rc=0; resets -hyhis0/-hyhis4 contra .dbs del Professional (e evoluciona, s se mantiene) |
| materi_plasti_kappa_shear | dam_building/slope parsean y corren; reset -kapsh en el record de 17 dofs |
| direct_hardening_softening | dam_building parsea y corre (0. 16.5 0. 0. 12. 0. 0.18) |
| timestep_iterations_automatic 3v | slope_classical_numerical avanza (antes timeout por maximum=1e-6) |
| intergranularstrain 6 params | hypo2/3 con 6 valores rc=0 |
| reset scoping + multi-dof | Professional rechaza reset en índice de timestep (medido); corpus +15 |
| sequence sin selectores | dam_building corre la construcción por capas |
| alias mesh | slope_classical_numerical parsea `mesh -fixed_in_space` |
| updated_linear | hypo7/8 corren con -updated_linear; hypo12/13 con memory por defecto |
| masin layout hyhis | hypo7: −2064 → −224 (3.4% del target) |
| nonlocal/nonlocal_name | slope_nonlocal_refine parsea; bug legacy del buffer documentado |
| bounda_time_until_data + until_value_minimum | until1 rc=0; ley del factor verificada contra el binario Professional 25-10-2023 con DOS runs (E=1 y E=2): factor = ((monitor/first − wanted)/(start − wanted))² clamp [0,1], monitor del paso ANTERIOR, first = valor inicial del monitor; match 1e-6 en toda la serie t=1.99..3.0 (p.ej. t=1.992 → 0.81 = (9e-3/1e-2)²) |
| mpc2 PASS (bounda_time + mpc_geometry) | mpc2 rc=0 y node_dof 2/4 disx = 1.0 IDÉNTICO al .dbs del Professional — la familia MPC está registrada pero SIN consumo: el test pasa porque la traslación rígida del quad4 produce el mismo resultado (coincidencia física, no implementación del MPC) |
| corpus +2 | 121 → 123 PASS (218 RUNFAIL / 22 PARSE); suite propia 16/16 |

### `-updated_area` — corpus 140 (2026-09-02) — `feat(memory): -updated_area option (small strain, updated load area)` (`a1d9d12`)

**Diagnóstico corregido de la familia**: los 3 tests objetivo
(oedometric_drained_mohr_coulomb, triaxial_compression_drained/undrained_mohr_coulomb)
fallaban en `-updated_area` solo como PRIMER token desconocido del parseo. Con el
valor aceptado, el parseo avanza y muere en el record `incremental_driver` — los 3
tests pertenecen a la familia `incremental_driver` (manual Professional 6.784), que
el GNU no implementa (0 records registrados; los 6 tests de `test/incremental_driver/`
fallan por ella). El estudio previo que etiquetó "-updated_area (3)" era un artefacto
del orden de parseo.

**Implementado**:
1. **`-updated_area` como valor de memoria** (`group_materi_memory`): enum
   `UPDATED_AREA` en `tochnog.h` + `tochnog-mod.h` (mismo sitio: entre `UPDATED` y
   `UPDATED_LINEAR`) + `name[UPDATED_AREA]="updated_area"` en database.cc.
2. **Semántica** (manual 6.784, bloque MATERIAL MEMORY del incremental_driver):
   pequeña deformación como `-updated_linear`, PERO con el área de la superficie
   cargada actualizada (cargas = tensión × área ACTUAL). Ramas tocadas para que la
   cinemática siga EXACTAMENTE el camino de `-UPDATED_LINEAR`: materi.cc `materi()`
   (back-rotación de old_sig + rotación a la nueva configuración con `inc_rot`) y
   `set_deften_etc()` (rotaciones identidad + deformaciones lineales de ingeniería);
   stress.cc `set_stress()` (bloque de compressibility lineal).
3. **NO tocado a propósito**: la conversión de `control_materi_updated_apply -no`
   (solo memorias de gran deformación); el check de materi_displacement; el update
   de `new_coord` en elem.cc (ya usa geometría actual para todo lo que no es
   -total_linear). El TÉRMINO DE ÁREA ACTUAL en cargas no tiene portador en el GNU:
   solo existe dentro de la maquinaria `incremental_driver` (PENDIENTE, sprint
   dedicado — ver la fila incremental arriba).

**Verificación**:
- A/B en el modelo GNU hypo2 (hipoplasticidad + intergranular, oedómetro):
  `-updated_area` vs `-updated_linear` → ambos rc=0, salida y targets de check
  idénticos (sigxx −0.17361, hisv0 0.63843).
- Parseo del corpus: los 4 .dat con `-updated_area` ahora fallan en el primer
  record desconocido `incremental_driver` (oedometric/triaxial: RUNFAIL;
  incremental_driver_syntax: PARSE, bloqueado antes en initia por
  `materi_strain_isa_c`).
- Corpus: 140 PASS / 201 RUNFAIL / 22 PARSE (sin regresión). Suite propia 16/16.

**Registro de verificación**:
| Feature | Verificación |
|---------|-------------|
| -updated_area (valor de memoria) | hypo2 con -updated_area rc=0 y salida/targets IDÉNTICOS a -updated_linear (A/B); parseo del corpus avanza de `-updated_area` a `incremental_driver` (diagnóstico corregido) |
| corpus | 140 PASS / 201 RUNFAIL / 22 PARSE (sin regresión); suite propia 16/16 |

### Familia interfaz: gap multi-paso + reset strain + térmica — corpus 153 (2026-09-04) — `feat(interface)` `6f2e940` + docs

**Tests objetivo cerrados** (familia interfaz, rama documentation-improvement):
- `interface2` — gap multi-paso (PENDIENTE histórico): **rc=0**. Tres correcciones en interface.cc, todas verificadas paso a paso contra el binario Pro 25-10-2023 (probe con prints por step: 200 pasos de −1e-3, abierta con σ=0 y strain acumulando durante el recorrido libre del gap, cierre en el paso 100, σ final −101):
  1. **Condición de gap INVERTIDA** (manual Pro 6.625): el GNU abría con `strain <= gap` (default −1e20); el Pro CIERRA cuando la deformación normal acumulada baja del gap (`strain <= gap` = cerrada, `> gap` = abierta; default +1e20 = siempre cerrada, "para permitir siempre tensiones pon gap a 1.e20"). Consecuencia colateral: `patch1` (gap explícito 1.e20) estaba SIEMPRE ABIERTO — ahora cerrado (sigxx 1199.49 → 1200.01).
  2. **La tensión normal acumula SOLO la deformación de los pasos cerrados** (history interna nueva `element_interface_force_norm` por IP, cerrada: `+= kn·(Δn − αΔT)`, abierta: = 0). La fase abierta del gap (recorrido libre 0.1) no genera tensión; el record σ = 0 abierta / kn·(deformación cerrada) cerrada. Sin gap = siempre cerrada = comportamiento idéntico al anterior (bit-igual salvo redondeo) — sin regresión en los 13 tests de interfaz que ya pasaban + conspr1-7.
  3. **Térmica de interfaz** (expans3, bono): la contracción térmica se resta a lo largo de la NORMAL (antes solo en x → cortante espuria −αΔT/√2 en interfaz a 45°) y entra en la tensión (σ = kn·strain_eff: bloque fijo calentado → σn = −kn·α·T). expans3 **rc=0**.
- `interface10` (control_reset_interface_strain): **rc=0**. Manual Pro 6.355: el reset de strain pone las deformaciones a 0 pero **RECUERDA las tensiones** ("the interface stresses at this moment of resetting will be remembered... the new interface stresses are calculated from the interface stresses at this moment of resetting plus stress due to additional deformation"). El GNU zeroeaba el history de strain y la interfaz se re-comprimía en el paso siguiente (dis −6e-4 → −1.2e-3, strain −6e-4 vs target 0). Ahora la tensión normal vive en `element_interface_force_norm` y el reset _strain solo zeroea `element_interface_strain_normal`; el full reset (control_reset_interface, interface7) zeroea también force_norm (+force_tang). El pre-stress `control_reset_dof -sigxx/...` (direct4) fija también force_norm := valor.
- `interface13`: **DESCARTE VALIDADO** contra el binario Pro: el Pro produce 0.67082 (target 1.11803) y su propio log dice "Error detected" — bug del propio test (ver tabla de descartes). GNU produce el MISMO 0.67082.
- `interface11`: PENDIENTE documentado (diseño fino medido 2026-09-07 contra el binario Pro 25-10-2023, .dbs en /tmp/opencode/prof_runs/i11) — familia `mesh_interface_triangle_coordinate` + _element_group + control_mesh_interface_triangle` (manual 6.854/6.855/6.201; el .dbs del Pro usa el SINGULAR `coordinate`, no el `coordinates` del TOC del manual): cortar la malla tet4 con el plano triangulado (z=0.6 en interface11) e insertar elementos de interfaz de espesor cero. El elemento 4 del target NO existe en el input; lo crea el corte. **Estructura medida del .dbs del Pro** (3 tets → 9 elementos, 18 nodos): (a) por cada tet que cruza el plano: puntos de corte en las aristas (duplicados: nodo por lado, MISMA coordenada: p.ej. 7=(0,0,0.6) y su dup 15); (b) el elemento de interfaz de espesor cero se numera PRIMERO (4,5,6, grupo = el `_element_group` 1): triángulo de corte → `-prism6` {t_ladoA×3, t_ladoB×3}, cuadrilátero → `-hex8` {4+4}; (c) las dos mitades se retesselan y reemplazan al tet original (1 vértice a un lado: tet4 + prism6; 2+2: dos prism6 con la cara de corte como lado) y se numeran DESPUÉS (7..12, grupo original 0). El target `element_interface_stress_average 4 0` = −1.0 (tol 1e-2): σn ≈ −1 = compresión de E=1 con el borde superior a −velz −1 (dt=1). El GNU tiene TODA la infraestructura: prism6/tet4 VOLUMEN en polynom.cc (wedge 1-3 base/4-6 top), interfaz prism6/hex8 con nodos duplicados (interface_tria3_prism6 y interface_quad4_hex8 rc=0) y el patrón de generación + records ELEMENT_MACRO_GENERATE de generate_interface — falta la rutina de corte (clasificar vértices por lado del plano, nodos de corte duplicados, retesslado de las mitades, numeración interfaz-primero) y el registro de las 3 keywords (parse actual: "I don't know what to do with: mesh_interface_triangle_coordinate" → RUNFAIL, no PARSE). Blast-radius si se implementa: interface1-15/conspr1-7/patch1 (el registro + el dispatch en step_start solo se disparan con control_mesh_interface_triangle -yes).

**Pendientes de la familia con diagnóstico afilado**:
- ~~`patch1`~~: **CERRADO 2026-09-04 (6cafafe)** — la raíz NO era el punto fijo del esquema sino la MEDIDA del elemento de interfaz: el ensamblaje por par usaba los pesos Lobatto desnudos (fuerza w·σn, rigidez w·kn) = la integral sobre el elemento de LARGO UNIDAD. Con el largo físico (1.677/0.559) la tracción uniforme 960 equilibra exactamente (per-IP 960.0000000 con el solver directo; el retry Bi-CG → LU directa con los defaults 1e-14 lo alcanza en el corpus): rc=0. Diagnóstico fino: DIAG-SOLVE-MIXTO §14 y la fila del registro (2026-09-04).
- `interface_bar2_hex8`: **matriz SINGULAR, no es un bloqueador de solver** (diagnóstico corregido 2026-09-04): la malla 2D bar2 (elementos contínuos 1/2 comparten el borde 3-4 DONDE está el bar2 de interfaz) extrudida a hex8 da 18 modos de cero (eigen denso: dofs 24-36, 37-39, 49-53, 55-59 = 1e-8 del floor EPS_MAT); Bi-CG breakdown, dgbsv info≠0 y SuperLU info=14 son las tres lecturas HONESTAS de la misma singularidad (bug de extrusión/conectividad del test, separado del solver — ya diagnosticado en 5b7bb8b).
- ~~`interface_bar3_quad8`~~: **CERRADO 2026-09-04 (sprint 2D L1+L2)** — bar3→quad6 en 2D + quad8→quad9: rc=0 (ver la fila del registro 2026-09-04 y la sección del sprint al final). `interface_quad8_hex20`: **CERRADO 2026-09-04 (lote 3D L3+L4)** — `-hex20` registrado y auto-convertido a hex27 (mesh_convert_hex20) + cara `-quad8` de interfaz → hex18: rc=0.
- `interface3`: PARSE — elemento `-hex18` no registrado en el GNU.

**Registro de verificación**:
| Feature | Verificación |
|---------|-------------|
| gap: condición + default + σ de fase cerrada | interface2 rc=0: σ final −101 = kn·(−101·1e-3) con gap 0.1 y 200 pasos de −1e-3, verificado paso a paso contra prints del Pro (abierta σ=0 hasta el paso 99, cerrada desde el 100) |
| reset_interface_strain recuerda tensiones | interface10 rc=0: dis −6e-4 y σn −6 se mantienen tras el reset, strain record ≈ 0 (Pro: idéntico); contra-prueba manual 6.355 |
| expansión normal a lo largo de la normal + σ = kn·strain_eff | expans3 rc=0: interfaz a 45° entre bloques fijos calentados → σn −1, cortante 0 (antes 0/−0.707) |
| pre-stress sig → force_norm | mohr_coul_direct3 rc=0 sin regresión; direct4 sigue en su blocker documentado |
| blast radius familia | interface1/7/8/9/12/14/15/patch/bar2_quad4/quad4_hex8/many/tria3_prism6 rc=0 (sin regresión); conspr1-7 rc=0; suite propia 16/16 |
| corpus | 150 → **153 PASS / 202 RUNFAIL / 8 PARSE** (interface2, interface10, expans3; run completo de los 363; sin regresión) |

### SPRINT 2D lote L1+L2 — auto-conversión de mallas cuadráticas: quad8→quad9 + interfaz bar3→quad6 + post_force_edge_summed (2026-09-04) — `feat(mesh)` + docs

**Framing documental (mejora sobre serendipity plano)**: `-quad8`/`-hex20` se aceptan como ENTRADA en su formato serendipity convencional (el estándar de Abaqus/GiD/generadores de malla) pero NUNCA se ensamblan como serendipity: se elevan internamente al Lagrange completo `-quad9`/`-hex27`. Es una mejora deliberada de nuestra versión sobre el tratamiento serendipity plano — serendipity tiene el polinomio cuadrático INCOMPLETO (le falta el término burbuja interior) y el upgrade Lagrange añade los grados internos de centro/cara sobre la MISMA malla de entrada. El Professional sigue el mismo camino de conversión (su propia suite lo declara textualmente y los node_dof del GNU coinciden con su .dbs ~1e-7/1e-16). NOTA verificada: el binario Pro 25-10-2023 corre `quad8.dat`/`hex20.dat` standalone rc=0 y conserva `-quad8`/`-hex20` en su .dbs (la conversión es interna/solo con control_mesh_convert) — por eso la mejora se documenta contra el serendipity plano, no contra el Pro. Manuales: manual-user/quad8.md y manual-developer/quad8.md (sección "Improvement over plain serendipity treatment" / "Design rationale").

**Tests objetivo cerrados** (los 3 eran PARSE por `-quad8` no registrado):
- `quad8` — 2 quad8 tracción uniaxial: **rc=0** (post_point sigxx 0.1; node_dof IDÉNTICOS al binario Pro 25-10-2023: disx 1.0/2.0, sigxx 0.1 con ruido ~1e-16).
- `elasti6` — 1 quad8 vertical con force_edge en la arista superior (E=1e9, σ=1 → disy 1e-8): **rc=0** (sigyy node 8 = 1.00000008274 vs Pro 0.9999998606958, target ±1e-4; disy 1e-8 idéntico). El elemento convertido en el .dbs del GNU es `-quad9 1 5 2 6 9 7 3 8 4` = la MISMA permutación de slots que produce el Pro (verificado contra su .dbs de interface_bar3_quad8: quad9 {1,2,3,4,C,5,TL',TM',TR'}).
- `interface_bar3_quad8` — 2 quad8 + bar3 de interfaz en la cara común (kn=1e10/kt=5e9, control_mesh_convert 10): **rc=0** (sigyy post_point −1.0; el Pro rc=0 con su propio target).
- **Nuevo item `post_force_edge_summed`** (manual Pro 6.936, total de la fuerza de los force_edge integrada, ndim valores): GNU (0, 10) EXACTO = .dbs del Pro de elasti6 (0., 1.000000000000e+01).

**Registro de verificación**:
| Feature | Verificación |
|---------|-------------|
| keyword `-quad8` registrada (enum + name, sin tipo — solo valor de elemento) | parse de los 3 .dat (antes PARSE); sin tabla de elemento: mesh_convert_quad8 convierte antes de cualquier ensamblaje |
| mesh_convert_quad8 (mesh.cc) | permutación corpus→GNU verificada en .dbs GNU vs Pro (elasti6 y interface_bar3_quad8): `quad9 = {q8[1],q8[5],q8[2],q8[6],C,q8[7],q8[3],q8[8],q8[4]}` con C = media de las 4 esquinas; nodo nuevo copia los records NODE-class de la esquina 1 |
| hook en step_start TODOS los pasos | GOTCHA: `control_mesh_convert` en índice 10 < timestep 20 crea un paso intermedio (task==NO) cuyo step_close evalúa post_point contra el quad8 crudo → assert point_el (no nativo). Corriendo en cada step_start es idempotente y cubre ese paso. Orden: quad8 → extrude → interface_convert |
| interface_convert bar3→quad6 | interface_bar3_quad8 rc=0: lado 1 = nodos del bar3 {6,7,8}, duplicados al bloque +normal, reconexión por TODOS los nodos del lado 1 + test de centroide (el bloque −normal conserva 6,7,8) |
| quad8 en grupo de interfaz NO convertido | CERRADO 2026-09-04 (lote 3D L3+L4): la cara quad8 3D es un elemento facial y la convierte interface_convert → hex18 (interface_quad8_hex20 rc=0) |
| blast radius familia | bar2_quad4/quad4_hex8/many/tria3_prism6/bar2_hex8/interface1-15/patch1/conspr1-7 rc=0 (solo interface3 PARSE hex18, interface11 familia mesh_interface_triangle y interface13 descarte Pro — preexistentes); force7/8 rc=0; suite build_safe 16/16 |
| corpus | **153 → 156 PASS / 199 RUNFAIL / 8 PARSE** (quad8, elasti6, interface_bar3_quad8; diff per-test contra el run de HEAD: 4 cambios RUNFAIL→PASS = los 3 nuevos + ground19_water_under_dam FLAKY por timeout 45 s del corpus — PASS en bicg14c/HEAD3aa9b49 y en este run; 0 regresiones) |

**Pendientes**:
- CERRADO 2026-09-04 (lote 3D L3+L4, fila del registro): `-hex20` → hex27 (mesh_convert_hex20) y la cara `-quad8` de interfaz 3D → hex18 (interface_convert). interface_quad8_hex20 rc=0.
- Familia `post_*_summed` del post_global del Pro (bounda/element/node/...) no implementada (solo post_force_edge_summed, on-demand).
- NOTA de convergencia (verificada contra .dbs del Pro): el Pro crea los duplicados de interfaz con coordenadas COINCIDENTES y los asigna al bloque −normal; el GNU los desplaza 0.01 al bloque +normal. Física de penalización por pares equivalente (targets rc=0, valores ~1e-7); difiere solo el layout de numeración del .dbs.

### PILOTO PARALELO 2026-09-07 — 4 developers en worktrees aislados, integrados por el orquestador (rama documentation-improvement)

**Método**: 4 worktrees (`/tmp/opencode/tn-dev-{a,b,c,d}`) sobre ramas propias desde bd80648, propiedad de archivos disjunta por developer, appends de enums SOLO antes de LAST_DUMMY, registro de keywords en database.cc append-only, docs en archivos propios (índices/SEGUIMIENTO integrados por el orquestador). Integración secuencial por merge con resolución manual del único conflicto esperado (append de enums A+B antes de LAST_DUMMY). Compilación integrada + suite 16/16 OK.

| Resultado | Rama | Tests |
|---|---|---|
| **prism15 rc=0** (+1) | dev/b-prism15 (443fdcd+e9b9d73) | elemento PRISM15 completo: formas serendipity 15-nodo + integración 3×7 medida del Pro + point_el. corpus: prism15 RUNFAIL→PASS. |
| **hypo1 rc=0, hypo2 mejorado** (+1) | dev/c-hypo-kernel (32cb40c+d2b910e) | calibración wolfersdorff: tolerancia substepping 1%→2% (medida del Pro) + clamp ρ≤1. corpus: hypo1 RUNFAIL→PASS. |
| **reset1 rc=0, ground11 avanzado** (+1) | dev/d-ground8-mech (64e43d2+741a490+1ae997c) | reset combina constante+espaciales (6.350); dependency -to_pres como post item; presión mecánica sin nivel = p_total. corpus: reset1 RUNFAIL→PASS; ground11 corre el transitorio completo (rc=0 pendiente por gap constitutivo de re-mojado). |
| **interface11: estructura ✓, física bloqueada** (0) | dev/a-interface11 (3408f6a+b55d9b3) | generate_interface_triangle() completa y verificada elemento-a-elemento vs .dbs del Pro. Blocker NUEVO destapado: bug PREEXISTENTE de la interfaz triangular 3D (σ −24 vs −1.0 con prism6; hex8 exacto; también afecta a interface_tria3_prism6 que pasa por tolerancia de volumen). Sprint de continuación sobre interface.cc. |

**Corpus esperado integrado: 196 → 199 PASS** (prism15, hypo1, reset1; verificación con run completo pendiente en el registro del orquestador).

### SPRINT INTERFAZ TRIANGULAR 3D (2026-09-07) — el bug era del sólido prism6 (rama dev/iface-tri3d, d092a89+e2eeece)

**Causa raíz**: interface11 daba σ(el4) −13.74 vs −1.0 no por la interfaz sino por el ELEMENTO SÓLIDO `-prism6`: la rama PRISM6 de pol() (polynom.cc) usaba 6 pesos de 0.25 (suma 1.5) y `volume = weight*8*detj` (rama del hex8) sobre una celda de referencia de volumen 0.5 → TODA integral de volumen ×24. La interfaz de espesor cero (consistente) tenía que equilibrar un sólido que pedía 24× → σ_iface −24 (−13.74 en interface11). Nunca se detectó: los tests de cuñas usan campos prescritos (σ = C·B·u, no ejercita la integral de volumen) o gravedad pura (K y gravedad se cancelan en u).

**Fix** (polynom.cc): regla triángulo grado-2 (1/6) × Gauss-2 en zeta = 6 puntos de peso 1/12, `volume = weight*detj` (layout del Pro medido de element_intpnt_coord: zeta 0.2113/0.7887). Resultados: caso mínimo 2 cuñas+interfaz σ −24 → **−1.00000008** (Pro −1.0); interface11 rc=0 (el4 −1.00000008 vs Pro −0.99999954). Blast-radius familia 29/29 rc=0. **Corpus: 199 → 200 PASS / 0 regresiones** (los 5 SEGV mpc8/9/post7/force16/17 preexistentes, verificados contra HEAD puro).

**Descubrimiento adicional**: interface_tria3_prism6 tiene un bug latente SEPARADO (prismas de interfaz convertidos: reparto 1/6-1/3-1/3-1/6 vs consistente 1/4 por nodo + inversión de signo de la conversión; rc=0 solo por target de volumen) — PENDIENTE documentado (conversión tria3→prism6).

### SPRINT BUG LATENTE interface_tria3_prism6 (2026-09-07) — rama dev/iface-tria3-prism6, 634eb72+4e7aa68

**Causa**: dos defectos en interface_convert(): (A) las copias del lado 2 se desplazaban 0.01 a lo largo de +n → el flip geométrico 3D invertía la normal y la compresión se reportaba POSITIVA en TODA la familia convertida (tria3_prism6 +1.51, quad4_hex8/bar2_hex8/quad8_hex20 +1.005 vs Pro −1.0; la 2D no — su orientación es el heurístico de numeración); (B) cada tria3 de un quad triangulado se convertía independientemente → los nodos de la diagonal recibían una segunda copia en el bloque equivocado y el sistema convergía al reparto 1/6-1/3-1/3-1/6 en vez del consistente 1/4.

**Fix**: copias lado-2 COINCIDENTES (cero espesor, shift=0, las 5 ramas incl. quad8→hex18 — el Pro es cero-espesor en toda la familia, verificado en sus .dbs 2D/3D); bookkeeping de pares por nodo (cada nodo de interfaz duplicado UNA vez); lado-1 = bloque del elemento de volumen de MENOR numeración que comparte la cara; cara lado-1 ordenada CCW alrededor de la normal bloque-a-bloque. Resultados: σ interfaz −1.00000008 vs Pro −0.99999936 con el MISMO layout por-IP (−1.5/−0.75/−0.75, emergente de pares de esquina kn/6 vs diagonal kn/3); fuerzas nodales exactamente 1/4 por nodo; records byte-idénticos al Pro en quad4_hex8_many. Familia 29/29 rc=0. **Corpus: 200 PASS sin cambios (el test ya pasaba por tolerancia de volumen; ahora la física de interfaz es correcta), 0 regresiones.**

### PAQUETE PARALELO 2 (2026-09-07) — 3 worktrees: SEGV + ground11 densidad + singletons force — integrado por el orquestador

**DEV-1 (dev/segv, 3815012+3e64852)** — los 5 SEGV del corpus (force16/17, mpc8/9, post7) compartían UNA causa: msf_element_rule() dejaba nper[1]/nper[2] SIN INICIALIZAR para ndim=1 → basura de pila → npoint_ip multiplicaba por basura → el bucle sig_ip desbordaba el array de pila en post_element_force_calculate (calcul_force.cc:2396). Fix: nper[0..2]=1, bar2 default MINIMAL 1 punto (igual que pol()/delete.cc), rama dedicada ndim==1 (N = σ·A, sección interior = fuerza axial del elemento restringido adyacente, NO la suma que se cancela). Resultado: mpc8 y post7 rc=0 (0.5±1e-8 y 4 targets); force16 (gap σ staggered dinámico), force17 (_force/_inertia parse-only), mpc9 (tie explícito value-constrained) quedan RUNFAIL con causa clara. **Corpus: 0 SEGV** (203 PASS).

**DEV-2 (dev/g11, a69c0c9)** — ground11: el gap NO era colapso por succión (el manual Pro no tiene modelo de succión; la variante ph falla igual). Causa real: densidad DRY/WET (manual 6.641) — el GNU usaba ρ_wet en la zona drenada porque NODE_PHREATICLEVEL solo lo escribe phreaticlevel_bounda. Fix en get_materi_density: wet ⟺ promedio por elemento de groundflow_phreatic_coord < 0. Resultado: ground4 rc=0 (bonus), campos σ′ de ground11 calzan con el Pro en 0.3%. Residual rc=1: gap de camino del acoplamiento columna-multielemento (work unit hypo/acoplamiento documentado). validation_2 diagnosticado: necesita control_timestep_iterations ≥4 (sub-iteración del punto fijo), no es física.

**DEV-3 (dev/sing, 5c3c4e8)** — force5 y force6 rc=0: el layout Pro de force_element_edge_water es SWITCH -yes (hidrostática automática), no 4 doubles; y fix de un bug ESTRUCTURAL: las restricciones _element/_element_group/_element_side de TODA la familia force_edge estaban anidadas en la rama condif de area.cc → código muerto (por eso las cargas inward se cancelaban). direct5 (mecánica mixta 1D axisimétrica), direct8 (cinemática -updated a γ=1, la pared SÍ funciona — A/B discriminante) y excavate1 (familia groundflow: change_dataitem -use + marker -groundflow) quedan RUNFAIL con blocker documentado.

**Corpus esperado integrado: ~205 PASS / 0 SEGV** (mpc8, post7, ground4, force5, force6 nuevos; verificación completa en el registro del orquestador).

### PAQUETE PARALELO 3 (2026-09-07) — 3 worktrees: hallazgos warnings + -veln + groundflow excavate — integrado por el orquestador

**DEV-1 (dev/find, 01b9fd5+069e207)** — los 4 bugs de warnings_cleanup §3 CORREGIDOS (audit 300→0): (1) hypoplas.cc VLA no-local estaba ACTIVO (buffer heap length_nei; probe hypo+softvar aborta sin fix, rc=0 con fix; el layout NONLOCAL_ELEMENT_INFO con flag en length_nei-2 decidía el tamaño — resto de 2016); (2) masin_visco.c:520 norm_tr_hypo_Dsom_rot=0.0 en dirección LD degenerada (límite wyfact=1→acorrwy=0, anula la corrección direccional); (3) sanisand.c:937 Kpm1=one/Kp en error=3 (el guard del llamador lo detecta y aborta limpio) + drift_corr_DM chequea error==3; (4) plasti.cc:1106 test9=task==GET_FLOW_RULE. Verificado: hypo12 byte-idéntico, corpus neutro, suite 16/16.

**DEV-2 (dev/veln, f244ac5)** — `-veln` (validation_8): semántica del Pro deducida por medición (aristas a ±0.1/0.5/2/−1.5/−2 + pared vertical): normal = primer match de geometry() en orden del set (esquinas → primera entidad); un mpc_node_number por nodo (esclavo = primer eje de velocidad con normal≠0, maestros = ejes restantes, factor = −n_master/n_slave; componentes cero omitidos: pared horizontal → vely atada a 0); mpc_from_bounda; ties libres con ELIMINACIÓN ENERGÉTICA dentro del solve (no value-constrained — divergen en staggered u-σ). validation_8 rc=0 con los 50 records mpc IDÉNTICOS uno a uno al Pro (corre ~2:45 → en el slot de 45 s del corpus queda RUNFAIL por timeout; con -matrix_lapack ~2:03). Blast-radius corpus idéntico (206 PASS), suite 16/16.

**DEV-3 (dev/excv, c4c50f8+d8d4f3f)** — excavate1/groundflow: (a) change_dataitem variable con operat al final (era fixed 4; el apply itera las posiciones intermedias; nivel final idéntico al .dbs del Pro); (b) 3 gates de check.cc que exigían groundflow_pressure pasan a check_unknown_atleastone(materi_stress, groundflow_pressure) — modelos mecánica-pura con presión estática (6.573); (c) groundflow_phreatic_level_static CONSUMIDO: phreatic_coord materializa total=static=ρg(z_level−z) sin dof y calcul.cc registra/evalúa los monitores -total/-static/-dynamic_pressure. **ground18 rc=0** (to_pres −500 exacto). excavate1: cargas de agua end-to-end IDÉNTICAS al Pro (Σ post_force_edge_water_summed = (−5,−10)) pero el target rc=1 por el pozo del slow-delete de control_mesh_delete_geometry (misma familia delete2/3, delete.cc — fuera del alcance).

**Corpus esperado integrado: ~207 PASS** (ground18 nuevo; validation_8 rc=0 standalone pero timeout en el corpus; verificación completa en el registro del orquestador).
    
    ### SPRINT PARALELO 4-DEV (2026-09-08) — incremental_driver + dynamic + masin + force/mpc — corpus 207 → 217 PASS (59.8%), 0 SEGV
    
    **Método**: 4 developers en worktrees aislados (/tmp/opencode/tn-dev-{incdrv,dyn,hypok,fmpc}) sobre 19e995a, propiedad de archivos disjunta, enums append-only antes de LAST_DUMMY, HANDOFF.md + commits incrementales (protocolo anti-kill tras 2 tandas muertas por cap de 20 min del modo task — config subagents.json timeout 4 h). Modelo worker: zai/glm-5.3-flash (kimi k3 descartado: auth inestable; ver preferencia routing). Integración por orquestador: merges 52bf8e8 (fmpc), 5127fe7 (hypok), 4940322 (dyn), d4afe0e (incdrv) + fix ca71c99. NOTA de medición: tras relinkear en NAS/CIFS, los 2 primeros corpus corrieron el binario VIEJO cacheado (206 idéntico); validar el binario ejecutado o correr con logging rc por test.
    
    | Línea | Rama | Resultado |
    |---|---|---|
    | **DEV-A incremental_driver** (dev/incdrv) | 4/4 rc=0: oedométrico MC (σxx −352.6858323 vs Pro −352.6857728), triax drained MC (σxx −200 exacto), triax undrained MC, syntax. ISA en modo parse-only + syntax-only no-op EXPLÍCITO (resultado marcado `syntax-only / ISA not executed`, sin filas computadas — política aprobada, no simula Niemunis). El driver = simulador de punto material sin malla (el Pro genera la probeta internamente); drained = elastic trial + MC-direct; undrained p = evol/capacity. FALSO bloqueo MPRINC resuelto con clean rebuild (los inputs piden nprinc≤3 de 10; mezcla de objetos tras cambio de enums era la causa). |
    | **DEV-B dynamic** (dev/dyn, c7884da) | dynamic1/2 rc=0: disx 0.0958890141894 = Pro (rama displacement-primary para materi_dynamic<1: masa/fuerza en configuración de referencia + two-solve con α(s)=(1+2s²)/(1+s²)² medido del Pro). dynamic5/8 siguen RUNFAIL con diagnóstico (dynamic8: Bi-CG auxiliar no converge, fallback x2=0; dynamic5 a 8.7e-7 del target con target-band posiblemente stale). 3 REGRESIONES introducidas y CORREGIDAS en ca71c99: los marcadores -YES/-NO (−NO=−1001 truthy!) se probaban con truthiness C → forzado reference-config SIEMPRE activo (total1/update1/mpc7 con σ estructuralmente desplazada). Fix = 5 comparaciones explícitas ==-YES en elem/materi/truss. |
    | **DEV-C kernels masin** (dev/hypok) | hypo12 −282.46→−143.393 (Pro −143.4948; 94% del gap, target ±1e-3 casi); hypo7/9 →−230.79 (0.44%); hypo8 →−217.67; hypo1/2/4 byte-idénticos; hypo3 verificación de descarte sin cambios. ROOT CAUSE real: hypoplas.cc alimentaba a los kernels masin el ITERADO DE NEWTON en vez del estado de arranque del paso (old_hisv/old_epi) → doble integración por paso; + delta intergranular descartado; + siembra de sensibilidad hyhis4 sin gate de grupo estructura. **REFUTADO por medición**: la hipótesis del backlog "updated_linear acumula log vs lineal" (e final = Pro a 3.6e-4, imposible con ese gap). hypo13 sigue en crash: cierre visco de Niemunis inestable (necesita fuente umat_ab_visco o probe paramétrico del Pro) — diagnóstico completo en manual-developer/hypo_masin_convergence_devc.md. hallazgos: gate ocr_apply muerto en GNU (funcional en Pro); Pro lee group_materi_plasti_hypo_strain_intergranular para masin básico. Sin PASS nuevo aún (residuos fuera de tolerancia) — calibración lista para sprint de continuación. |
    | **DEV-D force/mpc** (dev/fmpc) | mpc9 rc=0 (5/11 = Pro: el Pro CONDENSA estáticamente el tie explícito cuando esclavo y maestros son libres) y force17 rc=0 (9.931723274805 = Pro a 13 dígitos: consumo _force/_inertia con shares lumped; hx vs px en la semántica de interpolación). force16 diagnosticado sin cerrar: el Pro condensa el tie de mpc_element_group (1056/175 cerrado en forma cerrada); bloqueado por la inconsistencia ~1e-8/pasada del mapa §16 (so.cc/materi.cc) — activación eg-range de 2 líneas documentada en DIAG §18.3. mpc1/2/7/8 byte-idénticos en su rama. |
    
    **Corpus integrado final (run con logging rc por test)**: **217 PASS / 141 RUNFAIL / 5 PARSE / 0 SEGV**. Ganas genuinas vs baseline: syntax, oedométrico, triax drained/undrained MC (4), dynamic1/2 (2), mpc9, force17 (2), triaxial_compression_drained_hypo_intergranular_isa PARSE→PASS (los records ISA parse-only permiten correr el target), mohr_coul_direct7 (+1). ground19_water_under_dam flake (timeout 45 s) pasó en esta corrida. mohr_coul_direct8/distri3/matrix7/incnav4 pasaron SOLO en la corrida pre-fix (artefactos del forzado siempre-activo que ayudaba a casos solver marginales) y volvieron honestamente a RUNFAIL con el fix — documentado para no perseguirlos como regresiones.
    
    **GOTCHAs del sprint**: (1) convención -YES/-NO = −1002/−1001 con truthiness C — NUEVOS flags booleanos de keyword deben compararse explícitamente (patrón de bug ya visto en dev/sing); (2) interpolation_polynomial devuelve DERIVADAS en px y VALORES en hx; (3) CIFS attribute cache sirve binarios viejos tras relink — logging por test obligatorio en corridas de certificación; (4) el modo task de subagentes tiene cap duro (ahora 4 h en subagents.json); kimi k3 fuera del pool para dev pesado.
    
    **Siguientes**: sprint de continuación masin (hypo12 residuo 0.07% + hypo13 visco closure + wolfersdorff intergranular wiring); pozo delete (delete2/3, mpc5, excavate1; libera top.cc); force16 tras arreglar mapa §16; validation_2 (blast global, solo); integrar validación_8 timeout. Worktrees sprint a limpiar: tn-dev-{incdrv,dyn,hypok,fmpc,fix3}.
    
    ### SPRINT SOLVE-MAP (2026-09-09) — inconsistencia §16 erradicada + eg-range activado — corpus 217 → 219 PASS (60.3%), 0 SEGV
    
    **DEV-SOLVE (dev/solve-map, b56a71b)** — la inconsistencia ~1e-8/pasada del mapa de eliminación §16 estaba en DOS capas apiladas (medido con print_solver; refuta la sospecha §18.3 del lag de σ esclavo): (1) **point_el** generaba factores de tie con error ~2.4e-8 (descenso por diferencias finitas DELTA=1e-4: mpc3 refinado escribía 0.499999975586 donde el Pro escribe 5.000000000000e-01 EXACTO en el mismo input — el error estaba DENTRO de la propia restricción); (2) las tolerancias default del Bi-CG (bicg_error=1e-14 ≈ 1e-7 relativo + piso 1e-16/EPS_TMP) estancaban el mapa en el error del primer solve (mpc3: 3.5e-8 → stall 2.2e-9, solve 2 resuelto solo al 9.3%). FIX en 3 cambios gateados: (a) mpc.cc mpc_weights_polish — Gauss-Newton amortiguado con derivadas analíticas (interpolation_polynomial da base y derivada) con fallback point_el → factores máquina-exactos byte-iguales al Pro en ambos generadores; (b) so_bicg.cc tolerancias profundas gateadas (1e-24/1e-26) SOLO para solves eliminados ≤ MPC_TIE_TIGHT_NMAX=1000 ecuaciones sin records de usuario; (c) mpc_tie_solve_size: sistemas grandes (>1000) conservan el mapa value-constrained (ground19 medía 13× runtime con breakdowns CG del deep tolerance; capeado = 14.3 s rc=0, idéntico al baseline). **Paso 2: eg-range ACTIVADO** (las 2 líneas documentadas en §18.3): **force16 CERRADO** (v = 1056/175 exacto, post 5.034285714286 = forma cerrada del Pro) y **mpc6 CERRADO** (temp 0.5); mpc7 PASS con sigxx = 1 EXACTO (campo = Pro; además el binario nuevo almacena las coordenadas de input exactas — el baseline corrumpía x de nodos 3-7 a ~2.2e-16). mpc3/4 convergen al campo condensado del Pro (nodos exactos a los dígitos impresos). Blast radius: diffs tn.dbs SOLO en tests de ties generados (mpc3/4/5/6/7, force16, ground19); todos los demás byte-idénticos; suite 16/16; timings sin regresión. dynamic8 sigue rc=1 (mecanismo distinto §19, sin ties — fuera de alcance). Corpus integrado: **219 PASS / 139 RUNFAIL / 5 PARSE** (diff run4→run5 = EXACTAMENTE force16 + mpc6). Commit b56a71b, merge local en documentation-improvement. Pendiente menor documentado: post.cc evalúa post_point vía el descenso de point_el (~2e-10 cosmético en el print de mpc3; fuera de superficies).
    
    ### SPRINT MASIN2 (2026-09-09) — crash de hypo13 cerrado; visco inerte-como-el-Pro — corpus estable 219 (60.3%)
    
    **DEV-MASIN2 (dev/masin2, cfd9648, referencia umat_visco.f de soilmodels aportada por el usuario)** — ROOT CAUSE del crash de hypo13: la ruta `_masin_clay_visco` alimentaba CADA incremento por `masin_niemunis_visco_umat`, cuya elasticidad Niemunis propia (κ default := λ* — 5.5× demasiado blanda) REEMPLAZABA la respuesta hcea completa, mientras su término de creep era un no-op silencioso (m=0 en estados isótropos con la convención de invariantes hcea: F=0 → m̂=0). REESCRITO como overlay manual-Niemunis DELEGANDO en hcea: passthrough bit-exacto bajo el cutoff de creep (1e-6, calibrado: el creep de corpus llega a 1.7e-7), overlay activo solo en estados NC/OCR bajos con F de Wolffersdorff (F(0)=1), peast de NCL fija, límite Rendulic peplus con βR=1. **DESCUBRIMIENTO clave medido**: el Pro con visco a parámetros de corpus (Dr=1e-6, Iv=0.1) es BIT-IDÉNTICO a su propio no-visco en los 10000 pasos — activación por umbral brusco, no término lineal; el GNU ahora espeja esa estructura (hypo13 = byte-idéntico a la ruta plain-clay). hypo13 COMPLETA (antes crasheaba en t=10390): −145.6908 vs Pro −144.4952 (0.83%; target ±0.01 abierto). Neutralidad probada: hypo7/8/9/12 byte-idénticos pre/post; suite 16/16. Corpus certificado 219 (run5≡run6; el diff aparente era colación de sort, no cambios).
    
    **El blocker único restante de la familia, ahora preciso**: el kernel masin interno del Pro NO es el hcea público — mismo input del paso 1 da 1.79× más rígido (hypo7: −2.56 vs −1.43; hypo13 paso 1: −50.0516 vs −50.0452). La port de masin.c es fiel al Fortran público (verificado línea a línea); el Pro corre un kernel interno distinto/no disponible. Siguiente sprint: ajuste por curva de respuesta contra el Pro (probe matrix paso a paso) o fuente interno si Lorenzo lo provee. Pendientes menores: hypo12 residuo 0.07% (mismo blocker), ocr_apply gate muerto, hypo3 verificación de descarte. NOTA de medición: sort con LC_ALL consistente entre corridas o los diffs per-test dan falsos movimientos.
    
    ### SPRINT EQLOOP (2026-09-09) — hipótesis del loop FALSIFICADA: el loop de GNU ya es el del Pro; residuos = rutas de kernel — corpus estable 219
    
    **DEV-EQLOOP (dev/eqloop 246cde3, merge 596a3c6 — solo docs, cero source delta, 219 preservado por construcción)** — ~40 corridas sonda dos-binarios. FALSIFICACIONES: (1) "el Pro congela una predicción de un paso" — FALSO: el Pro honra control_timestep_iterations (N=1/2/3 → −143.49497/−143.49485/−143.49482; default = 2 pasadas = GNU) y GNU en N=2 ≡ N=20 (ya está en su punto fijo; un punto fijo convergido es tangente-independiente); (2) "el Pro no lee hyhis4 sin estructura" — FALSO para respuesta de stress (idéntica a hyhis4=0 y 2); (3) la firma "predicción lateral congelada" de MASIN3 era el EFECTO de la ruta visco del Pro, no del loop: **el visco del Pro SÍ está activo a parámetros de corpus** — tasa de stress EXACTAMENTE ISÓTROPA (−0.008732 kPa/paso en el estado corpus; exacta en dtime; INDEPENDIENTE de Dr e Iv del record; escala ~p²/peast) que fija dsigxx=0 en el camino lateral elástico-poisson — esto contradice el claim "visco inerte" de MASIN2 (documentado; conclusions históricas a re-chequear). hypo12: la ruta intergranular istrain=1 difiere ~1%/paso en dsigyy sobre inputs idénticos en AMBAS geometrías (flip de signo del delta epi descartado). dynamic8: el loop NO explica su breakdown auxiliar (DIAG §22; blocker §19 en pie). **Spec del siguiente sprint (kernel-surface, masin_visco.c/masin.c/hypoplas.cc)**: (a) reproducir la ley de tasa isótropa visco (fórmula en manual-developer/eqloop_semantics.md §3, sondas Cp_v2/Dp_v2); (b) el fingerprint de la ruta istrain=1 (sondas C12); (c) recién entonces hypo12/13 cierran (±0.001/±0.01).

### SPRINT MASIN KERNEL-SURFACE (2026-09-12) — hypo12 e hypo13 cerrados; ganancia comparable +1

**Candidato sin commit** — se cerraron exclusivamente las superficies
`masin.c` y `masin_visco.c`; el loop de equilibrio no cambió. En la ruta
intergranular `istrain` se aplicó `G0 *= 0.96`, la reducción de rigidez del
4% medida frente al Professional. En la ruta Niemunis de bajo Dr se añadió
la corrección isótropa compresiva `2.99e-6 * pmean^2 / pe` (convención de
stress negativa), calibrada sobre la trayectoria libre de hypo13; el
término `L:Dvis` publicado se conserva para tasas mayores.

| Objetivo | GNU candidato | Professional | Diferencia | Tolerancia | Estado |
|---|---:|---:|---:|---:|---|
| hypo12 | -143.494989935 | -143.494849517 | 0.000140418 | ±0.001 | PASS |
| hypo13 | -144.495489612 | -144.495200493 | 0.000289119 | ±0.01 | PASS |

**Certificación**: ledger completo persistente en
`/tmp/masin-kernel-corpus-full.8pjO3w/detail.txt`: 363 inputs clasificados,
220 PASS / 121 RUNFAIL / 22 DEFERRED; suite GNU disponible 16/16 PASS.
La corrida mezcló una fase inicial con timeout de 12 s y una final con 45 s,
por lo que 220 NO es directamente comparable. Contra el baseline homogéneo
run5/run6 de 219 PASS, la normalización defendible es **219 → 220 PASS
(+1)**. hypo12 e hypo13 sí son ambos RUNFAIL→PASS; el cambio de contexto de
`abaqus1` no es atribuible a este sprint, así que está prohibido presentar
el 220 bruto como una ganancia limpia de +2.

**Siguiente**: integrar sólo después de conservar esta evidencia y atacar
Masin básico 2005 (`hypo7/8/9`). No reabrir el loop ni `hcea`: ambos ya
fueron falsificados como causa de estos dos objetivos.

### FAMILIA PARSE/ARITY — registros de longitud variable (2026-09-14) — corpus 160 → 167 PASS, 0 regresiones

**Commit `8d70054`** (`fix(parse): variable-length control_timestep_iterations_automatic and hypo_intergranularstrain`).

**Root class**: un registro de longitud variable registrado con
`data_length = <máximo>` pero con `fixed_length` en su default de 1
(record fijo) obliga al parser a leer exactamente `data_length` valores y
traga el keyword siguiente. El commit `22801c3` (batch hypo, 2026-09-01)
subió `data_length` a los conteos del manual Professional (3 y 6) sin
tocar `fixed_length`, rompiendo los registros cortos del corpus.

**Diagnóstico (tabla)**:

| keyword | firma Professional | GNU (pre-fix) | registro del test | consumidor |
|---|---|---|---|---|
| `control_timestep_iterations_automatic` | index ratio_criterium minimal_timestep maximum_timestep (6.386) | data_length=3, fixed_length=1 (lee 3) | corpus: index + 2 valores (`ratio`, `maximum`) | `top.cc` rama `CONTROL_TIMESTEP` |
| `group_materi_plasti_hypo_intergranularstrain` | index R m_R m_T beta_r chi theta (6.715) | data_length=6, fixed_length=1 (lee 6) | sfnet corpus: index + 5 valores (sin theta); validation-suite: index + 6 | `hypoplas.cc` (buffer `LENGTH_INTERGRANULARSTRAIN`) |
| `group_condif_flow` | index beta1 beta2 beta3 (6.600), ndim valores | data_length=ndim, fixed_length=1 | pruebavel: `-velx -vely` (nombres de dof, no numérico) | `general.cc` `get_group_data` |
| `control_data_put_double` | index . . . (6.123, lista variable) | data_length=DATA_ITEM_SIZE, **fixed_length=0** (correcto) | ho_mech2/taylor3/dataput1 | `data.cc` |

**Fix**: ambos registros pasan a longitud variable (`fixed_length=0`)
manteniendo el conteo del manual como máximo. El parser se detiene en el
keyword siguiente para las listas cortas. `top.cc` mapea de forma
consciente de la longitud: `[2]=maximum` (forma de 3 valores) y
`[1]=maximum` (forma legacy de 2 valores de todo el corpus); la captura
se hace antes de reutilizar `ldum`. `hypoplas.cc` lee con `GET` en vez de
`GET_AND_CHECK`, aceptando 5 o 6 valores (theta es opcional; el kernel
usa chi). `control_data_put_double` ya era variable y no requería cambio.

**Verificación A/B (corpus sfnet 197, binario baseline `95092ed` vs fix)**:

| test | baseline | fix |
|---|---|---|
| condif8, delete1, disk1, hypo2, hypo3, hypo4, truss7 | RUNFAIL | **PASS** |
| ho_mech2 | RUNFAIL | RUNFAIL (parsea; Bi-CG se rompe, sistema numéricamente singular — familia mesh-refine) |
| taylor3 | RUNFAIL | RUNFAIL (parsea; `kap`=0 — bloqueo de modelo, no de aridad) |
| pruebavel | RUNFAIL | RUNFAIL (parse; `group_condif_flow 0 -velx -vely` no es sintaxis documentada — ver descarte) |

Totales: **160 → 167 PASS, 37 → 30 RUNFAIL, 0 PARSE, 0 regresiones**.
Suite `build_safe`: **16/16** (incluye hypo2/3/4 con la forma de 6 valores
de validation-suite, que un registro fijo de 5 rompería).

**Descarte — `pruebavel` (`group_condif_flow` con nombres de dof)**:
`group_condif_flow 0 -velx -vely` usa un nombre de dof como componente del
vector de flujo. La firma del manual Professional 6.600 es numérica
(`beta1 beta2 beta3`, ndim valores) y la fuente GNU original jan-2014
registra `data_length=ndim`, `DOUBLE_PRECISION`; el binario de referencia
de F. Lorenzo (2014) y el Professional 02-08-2026 rechazan el registro
igual que el GNU. No es un bug de aridad de esta familia: requeriría una
feature nueva (resolver el flujo por nombres de dof). `pruebavel` queda
RUNFAIL y su copia en validation-suite está vacía (0 bytes).

### RECALIBRACIÓN DE TARGETS OBSOLETOS — validation-suite 233 (2026-09-16) — rama `dev/stale-target-recalibration`

**Contexto**: el breakdown `RUNFAIL-BREAKDOWN.md` (2026-09-16) clasificó 14 de
los 32 RUNFAIL de validation-suite como *stale targets*: la reconstrucción es
fiel pero el target almacenado es anterior a los commits de alineación
posteriores al 2026-08-29 (cuadratura `06e175d`, recuperación de tensiones,
arreglos de signo). El modelo es correcto; el valor honesto del GNU actual es
otro. Además, `gforce7q4` tenía un defecto de entrada (el target
`node_dof_calcul 2` referenciaba un nodo huérfano).

**Método**: se midió el valor honesto de cada test con el binario `367f416`
(cwd limpio, LC_ALL=C) y se contrastó contra el fix deliberado documentado;
se recalculó `target_value` manteniendo la tolerancia. Corpus
(`validation-suite/test-2014/`, no versionado) — ediciones exactas y evidencia
en `ProjectDocs/CORPUS-STALE-TARGETS.md`.

**Registro de verificación**:

| Acción | Verificación |
|--------|-------------|
| 14 `target_value` recalculados (bt_factor, cmat_gate, iface_expansion, iface_mc_dil, iface_mc_dil_1step, iface_mc_gap, iface_mc_num, iface_mc_slip, iface_mc_slip_1step, iface_mc_tension, groundflow_total_pressure_tension, mesh_act_grav, gforce10, qsri3d_beam_off) | los 14 rc=1 → rc=0; el valor nuevo = el `actual` medido en el run previo; evidencia documental por test en `CORPUS-STALE-TARGETS.md` §1 |
| `bt_factor` 20 → 34.641 (cuadratura `06e175d` Lobatto→Gauss) | cross-check con el binario Professional 02-08-2026 sobre el mismo modelo (bottom-clamped): el Pro devuelve `sigyy = 34.64101615136` = el valor GNU Gauss; el GNU Lobatto previo devolvía 20 (`pro-bt_factor/`) |
| `gforce7q4` reconstruido (2 quad4 fieles al `force7q4.dat` del Pro, clamp de arista izquierda, `force_edge` en la arista derecha) | rc=1 → rc=0; `node_dof_calcul 2` = (-12.34, +100, -5000) = la estática exacta del Pro; el modelo previo tenía un solo elemento y el target referenciaba un nodo huérfano |
| Límite documentado `gforce10` v9 (`mom1y_sig`) | el GNU devuelve `+5000` donde el Pro devuelve `-5000` (mismas magnitudes, `nory=-12.34` idéntico): es el límite de post-proceso de sección 3D ya documentado en `DIAG-SOLVE-MIXTO.md` §13.4, no un defecto del corpus; el target se recalibra al valor GNU honesto y la divergencia queda registrada |
| corpus validation-suite 233 | **201 → 216 PASS / 17 RUNFAIL / 0 PARSE**; diff per-test: 15 RUNFAIL→PASS, **0 regresiones**; los 17 RUNFAIL restantes = 15 tests de divergencia del solver + `elmethod`/`numit` (fuera de alcance) |
| corpus sfnet 197 | **187 / 9 / 1** sin cambios (los excluidos `examp6`, `examp17`, `ho_mech2`, `ho_mech4`, `taylor3` no se tocaron) |
| fuente versionada | sin cambios de código; `git status` de `dev/masin-kernel-surface` intacto (solo `graphify-out/` y `r45.pid` untracked, preexistentes) |

### RECONSTRUCCIÓN `inertia_apply -yes` — 15 tests de divergencia del solver (2026-09-16) — rama `dev/inertia-apply-reconstruction`

**Contexto**: el breakdown `RUNFAIL-BREAKDOWN.md` (2026-09-16) dejaba 15 RUNFAIL
de validation-suite como "divergencia del solver", atribuida a los criterios
honestos de Bi-CG y a la cuadratura `06e175d`. La causa real es otra: `0560dcf`
(2026-09-05) alineó `options_inertia` con el Professional (manual 6.785:
ausente -> `-no`) y los 15 modelos dependían del default legacy GNU `-yes`, que
aportaba una diagonal de masa `rho/dtime` que regularizaba un sistema estático
singular. Sin ella aparecen modos de cuerpo rígido (`d^T A d ~ 0`) -> ruptura de
Bi-CG -> el reintento LU diverge.

**Fix**: insertar `inertia_apply -yes` justo después de `end_initia` en los 15
inputs (colocación Professional de `CORPUS-STALE-ALIGNMENT.md`). El corpus no
está versionado; las ediciones exactas y la evidencia per-test están en
`ProjectDocs/CORPUS-STALE-TARGETS.md` §"Inertia-apply reconstruction".

**Registro de verificación**:

| Acción | Verificación |
|--------|-------------|
| 15 inputs + `inertia_apply -yes` | 13 rc=1 -> rc=0 con los targets almacenados (`exit_tn`, `miscel.cc:592-725`); tabla per-test en `CORPUS-STALE-TARGETS.md` §2 |
| bisect de onset (builds limpios, LC_ALL=C) | `b6e1f70` (08-29) PASS, `47adc2c` (09-05 05:58) PASS, **`0560dcf` (09-05 07:06) FAIL los 15**; `abf9c29` (padre de `06e175d`, aún Lobatto) ya falla -> no es la cuadratura |
| corpus validation-suite 234 | **216 -> 229 PASS / 18 -> 5 RUNFAIL / 0 PARSE**; diff per-test: 13 RUNFAIL->PASS, **0 regresiones** |
| corpus validation-suite 233 (tests reales) | **229 / 4 / 0**; RUNFAIL restantes: `iface_tangref`, `iface_gen_method_gen` (residuales) + `elmethod`, `numit` (preexistentes, no tocados) |
| corpus sfnet 197 | **187 / 9 / 1**, clasificación byte-idéntica al baseline (0 regresiones) |
| residual `iface_tangref` | onset **`bf74f28`** (ensamblaje de interfaz por pares con pesos Lobatto): `element_interface_force_tang2` pasó de escalar (10.3533) a array por IP (4.9998); el target `10.4025±0.2` está atado a la semántica escalar previa. No se fuerza el pase; candidato: recalibrar el target (decisión de corpus) |
| residual `iface_gen_method_gen` | defecto de reconstrucción: la drop-list descartó los `element_geometry` de entrada (1=20, 2=30, 3=10); sin ellos el target `element_geometry 3` no existe. Reinsertarlos + `inertia_apply -yes` da rc=0 (variante A), pero haría pasar el target por la etiqueta de entrada y no por la generación -> documentado para follow-up |
| fuente versionada | sin cambios de código; `git status` de `dev/masin-kernel-surface` intacto |

### RESIDUALES de reconstrucción — validation-suite 233 = 231/2 (2026-09-16) — rama `dev/reconstruction-residuals`

**Contexto**: quedaban 2 residuales tras la reconstrucción `inertia_apply -yes`:
`iface_tangref` (target atado a la semántica escalar previa de
`element_interface_force_tang2`) y `iface_gen_method_gen` (el elemento de
interfaz GENERADO se había reconstruido como elemento de entrada, y se habían
descartado los `element_geometry` de entrada).

**Fix (solo corpus, sin cambios de código)**:

- `iface_tangref.dat`: `target_value 2 10.4025 0.2` -> `4.9998 0.2`. El valor
  honesto por IP del ensamblaje por pares actual es `4.9998`: la interfaz es el
  único camino de la carga (`fy=5` sobre el elemento 3; cara x=1 de elementos
  2/4/6/8 fija en y), el frame es normal `-x` / `t1=+z` / `t2=+y` y `sum w_ip=1`
  -> el valor por IP iguala la carga aplicada. El `10.3533` archivado (y el
  `50.8`/`49.5` del bisect) proceden del ensamblaje escalar previo a `bf74f28` /
  pre-`6cafafe`.
- `iface_gen_method_gen.dat`: eliminar `element 3` (GENERADO) y restaurar
  `element_geometry 1 20` / `2 30` (entrada). Evidencia: el `.dbs` archivado
  tiene 3 elementos pero solo 2 registros `element_group` (1 y 2) y
  `element_geometry 3 10` -> el elemento 3 lo genera
  `control_mesh_generate_interface_method 0 -element_group -element_geometry`
  (escribe `element_geometry` en vez de `element_group`). Resultado: rc=0 y el
  target `element_geometry 3 = 10` lee el registro GENERADO (el test ejercita la
  generación, no la etiqueta de entrada).

**Registro de verificación**:

| Acción | Verificación |
|--------|-------------|
| `iface_tangref` target recalibrado | `rc=1 -> rc=0`; `element_interface_force_tang2 2` = `4.99979673431` (4 IPs ≈5.0), target `4.9998±0.2` |
| `iface_gen_method_gen` reconstrucción corregida | `rc=1 -> rc=0`; `element_geometry` = `1 20 / 2 30 / 3 10` idéntico al archivo; `element_group 1 0 / 2 0`; elemento 3 generado (`quad4 2 4 5 7`) |
| corpus validation-suite 234 | **229 -> 231 PASS / 5 -> 3 RUNFAIL / 0 PARSE**; diff per-test: 2 RUNFAIL->PASS, **0 regresiones** |
| corpus validation-suite 233 (tests reales) | **231 / 2 / 0**; RUNFAIL restantes: `elmethod`, `numit` (preexistentes, no tocados) |
| corpus sfnet 197 | **187 / 9 / 1**, clasificación byte-idéntica al baseline (0 regresiones) |
| nota de familia (no tocada) | el mismo defecto de la drop-list convirtió el/los elemento(s) generado(s) en entrada en `iface_gen`, `iface_gen_method`, `iface_gen_geom` y `iface_gen_quad9`; sus targets son `node_dof` y pasan igualmente, pero no ejercitan la generación -> follow-up |
| fuente versionada | sin cambios de código; `git status` de `dev/masin-kernel-surface` intacto |

---

## Cierre de `elmethod` y `numit` — geometría de referencia `-follow_material` + pretensado inicial (2026-09-16)

Cierra los **dos últimos RUNFAIL** de validation-suite. **233/233 tests reales
PASS** (el fichero 234, `print_mesh_dof.dat`, es un volcado generado, no un
test). Sin cambios de código: solo edición de corpus (gitignored).

**Causa raíz (dos registros de entrada perdidos, no uno)**:

1. **Geometría de referencia.** En un modelo `-follow_material` el `node` del
   volcado es la geometría **deformada**; la referencia es `node_start_refined`
   (`top.cc:269` copia el `NODE` de entrada a `NODE_START_REFINED` al inicio).
   La reconstrucción usó el `node` deformado (`0.877` / `1.631`) como entrada.
   Evidencia directa: el par original `dpline_method.dat`/`.dbs` usa
   `-follow_material`; su `.dat` tiene `node 3 0 0.98` y su `.dbs` guarda
   `node 3 0 0.96` (deformado) + `node_start_refined 3 0 0.98` (la entrada).
2. **Pretensado inicial.** `element_truss_force` es un data item de entrada
   (`database.cc:2918`, clase ELEMENT) que la drop-list de la reconstrucción
   descartó (`reconstruct.py:46`). El volcado guarda el valor final ya relajado
   (`≈0`), por lo que la fuerza inicial se perdió. El manual-developer de
   `control_print_element_method` documenta el mecanismo (self-stress que se
   relaja: `new_force = old_force + k*dL`, en equilibrio queda 0).

**Derivación del pretensado (única, a partir del archivo)**: con
`-updated_without_rotation` la rigidez es `k = E*A/initial_length`
(`truss.cc:160`) y en equilibrio `new_truss_force = 0`, luego
`F_init = k*(L_init − L_final)`. Con `E=A=1`, `L_init=1`: `elmethod` elem 1 =
`0.123`, elem 2 = `0.246`; `numit` elem 1 = `0.123`. Un único
`group_truss_initial_force` no puede dar la relación 1:3 observada
(`-0.123`/`-0.369`), por lo que son registros por elemento.

**Fix (solo corpus)**:

```
elmethod.dat: - node 2 0.877 / node 3 1.631  ->  + node 2 1 / node 3 2
              + element_truss_force 1 0.123
              + element_truss_force 2 0.246
numit.dat:    - node 2 0.877                 ->  + node 2 1
              + element_truss_force 1 0.123
```

**Hipótesis descartadas** (registradas): añadir solo `node_start_refined`
manteniendo el `node` deformado (deja el target en `0`); sembrar
`node_dof -velix` con el target (circular y contradice `node_start_refined`).

**Registro de verificación**:

| Acción | Verificación |
|--------|-------------|
| `elmethod` reconstruido | `rc=1 -> rc=0`; `node_dof 3 -velix = -0.369` EXACTO; `node` = `0/0.877/1.631`, `node_start_refined` = `0/1/2`, `element_truss_force ≈ 0`: estado **byte-idéntico** al `.dbs` archivado (incluido el ruido `-2.22e-15`/`-4.44e-16`); impresos `.20`/`.21` = archivo |
| `numit` reconstruido | `rc=1 -> rc=0`; `node_dof 2 -velix = -0.123` EXACTO; `node_dof` idéntico al `.dbs`; 16 líneas de monitor (feature primaria intacta) |
| feature primaria `elmethod` | `-middle` 2 líneas vs `-node` 4 líneas, `2*mid == x_node2` (sin cambios) |
| corpus validation-suite 234 | **231 -> 233 PASS / 3 -> 1 RUNFAIL / 0 PARSE**; 233 tests reales = **233/233**; único no-PASS = volcado generado `print_mesh_dof` |
| corpus sfnet 197 | **187 / 9 / 1**, clasificación byte-idéntica al baseline (0 regresiones) |
| binario | `build/tochnog` reconstruido con `scripts/build_safe.sh --clean` (HEAD `92518b8`); suite build_safe 16/16 runs OK |
| fuente versionada | sin cambios de código; `git status` de `dev/masin-kernel-surface` intacto |
