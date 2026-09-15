# Corpus reconstruction from archived `.dbs` dumps

Recovery of the 194 validation-suite tests whose `.dat` inputs were lost on
2026-08-29 (`rm -f *.dat`). Reconstructed from the archived `.dbs` database
dumps with the current binary (`build/tochnog`, HEAD `e983be8`, branch
`dev/masin-kernel-surface`).

> **Update (2026-09-14, third pass — doc-grounded).** The 8 UNRECOVERABLE
> and the 38 PARTIAL were revisited against `ProjectDocs/manual-user` /
> `manual-developer`. Five of the 8 are now recovered (three by a
> reconstruction fix, two by the current documented syntax); three PARTIAL
> were reclassified as reconstruction issues and recovered; four genuine
> binary bugs were isolated (one already reported here, three new). See
> **Third pass** at the end of this document. Evidence:
> `/tmp/opencode/docrecover-20260914-233512/`.

Evidence directories:

- `/tmp/opencode/reconstruct-20260914-184915/` — first pass
  (scripts `work/reconstruct.py`, `work/recon2.py`; per-test runs under `runs/`).
- `/tmp/opencode/recover-partial-20260914-194220/` — second pass, recovery of
  the PARTIAL set (scripts under `work/`; per-test runs under `runs/final/`).

## Method

`tochnog` dumps the whole external database at the end of a run (`print.cc` /
`print_db.cc`, `data_number=-DATABASE`). The dump mixes the parsed input
records (canonicalised, comment-free, reordered) with records the run itself
produced, plus a few lost scalars.

1. **Parser grounding.** The input grammar was read from `input.cc` (the
   init-section keyword list and the data-part `db_number()` lookup). Input
   records are stored verbatim; the dump line format (`keyword [index]
   values`) already matches what the parser expects, so kept lines are
   emitted as-is.
2. **Derived drop-list.** 103 record names are dropped. The core is the 45
   names that appear in the 41 known `.dat`/`.dbs` pairs in *no* `.dat`
   (state/results: `node_dof*`, `element_dof*`, `node_rhside*`,
   `node_start_refined`, `element_matrix*`, `dof_*`, `time_*`,
   `global_solver_*`, `post_*`, `tochnog_version`, `Size`, `Total`, ...),
   extended with the run-produced records seen in the wider dump set
   (`element_beam_*`, `element_interface_*`, `element_mass/volume`,
   `element_geometry`, `crack_stressintensityfactor`, `inverse_*`, ...).
   Output provenance was confirmed in the source (e.g. `elem.cc` PUTs
   `ELEMENT_VOLUME`/`ELEMENT_MASS`, `crack.cc` PUTs
   `CRACK_STRESSINTENSITYFACTOR`).
3. **Initial-state recovery.** `node_dof_start_refined` is the dump's copy of
   the start-of-run DOF (it holds the input `node_dof`). It is re-emitted as
   `node_dof`, otherwise hypoplasticity tests abort with "Illegal porosity".
   For full-mesh-generation tests (`control_mesh_macro`,
   `control_mesh_generate_beam`, `control_mesh_generate_truss`) the generated
   `node`/`element`/`element_group`/`node_start_refined` are dropped instead.
4. **Lost scalars — derived, not invented.** The init-section
   `number_of_integration_points` is stored bare in the dump (its value is
   never written). The value is recovered from the archived `element_dof`
   record length (see below) and re-inserted into the reconstructed input.

## Derivation of the lost `number_of_integration_points`

The lost value is recoverable from archived evidence alone.

In `elem.cc` the element state vector is allocated with a fixed stride

```
mnolnuknwn = npointmax * nuknwn
```

and `db( ELEMENT_DOF, ..., mnolnuknwn, ... )` dumps exactly that many values.
Here `npointmax` is the input `number_of_integration_points` and `nuknwn` is
the number of element unknowns per integration point (a function of the
`derivatives` init set and `ndim` only, independent of `npointmax`). The dump
also contains the `dof_label` record, whose length is `nuknwn - 1`
(the label array has one unused trailing slot). Therefore:

```
nuknwn = len(dof_label) + 1
number_of_integration_points = len(element_dof) / nuknwn
```

Both inputs are read from the archived `.dbs`; no fresh run is required and no
value is invented. The result is accepted only when the division is exact and
the element types are consistent with the resulting count.

### Cross-check with a fresh probe

As an independent check the reconstructed `.dat` (which lacks the record) is
re-run with an explicit `number_of_integration_points K` (K=27). Then
`len(element_dof) = K * nuknwn`, so `nuknwn = len(element_dof)/K` measured on
the current binary. The two nuknwn estimates agree on **72/72** tests, and the
derived counts agree on **72/72**.

### Validation on the known `.dat`/`.dbs` pairs

11 of the 41 known pairs carry the record; for these the true value is known.
The archived-only derivation reproduces it exactly in **11/11** cases:

| test | true | dof_label | nuknwn | element_dof | derived | match |
| --- | --- | --- | --- | --- | --- | --- |
| gforce7 | 9 | 71 | 72 | 648 | 9 | yes |
| gforce7q4 | 4 | 71 | 72 | 288 | 4 | yes |
| msf_beam2d | 9 | 15 | 16 | 144 | 9 | yes |
| msf_cant3d_hex27 | 27 | 104 | 105 | 2835 | 27 | yes |
| msf_sheet3d | 27 | 89 | 90 | 2430 | 27 | yes |
| msf_tunnel3d | 27 | 89 | 90 | 2430 | 27 | yes |
| mstrain_cap_elast | 8 | 149 | 150 | 1200 | 8 | yes |
| mstrain_compression_elast | 4 | 111 | 112 | 448 | 4 | yes |
| tpef_force10 | 8 | 104 | 105 | 840 | 8 | yes |
| tquad6_iface | 9 | 15 | 16 | 144 | 9 | yes |
| tsup_3d | 8 | 104 | 105 | 840 | 8 | yes |

The remaining 30 known pairs have **no** `number_of_integration_points` record
in either file; applying the formula there returns the archived binary's
*default* count (6, 27, ...), not an input value, so those are correctly left
untouched (they were never flagged PARTIAL for a lost value).

### Applied values (72 tests)

`name(derived value)`:

```
cd_geom(4) cd_method(4) iface_3d_order(8) iface_conv3d(8) mat_rel(4) mat_rel_reset(4)
materi_direct(4) materi_direct_auto(8) materi_direct_mc(4) materi_direct_visco(4)
materi_direct_wall(4) mcap1(8) mcap1_comb(8) mcap1_elast(8) mcap2(4) mcap_legacy(4)
mcrunch(4) mcrunch_low(4) mdirect_comp(4) mdirect_gate(4) mdp_shear(4) mesh_act_grav2(4)
mhardsoil_elast(4) mhardsoil_elast2(4) mhardsoil_gp0(4) mhardsoil_gp0_off(4)
mhardsoil_plast(4) mhardsoil_plast_elast(4) mhardsoil_unload(4) mhardsoil_unload_flat(4)
mk0(4) mk0_off(4) mmc_tension(4) mmchs_soft(4) mpower(4) msf_beam2d_pure(9)
msf_hex27_avg(27) msf_parse_3d(8) msf_quad9(9) msf_quad9_noavg(9) msf_sheet3d_hex8(8)
mshf(4) mshf_nof(4) msph(4) msph_flat(4) mstrain_cap(8) mstrain_compression(4)
mstrain_druckprag(4) mstrain_druckprag_elast(4) mvoid(4) mvoid_low(4) myoung6(4)
myoung6_apply(4) myoung6_e2(4) myoung6_e3(4) qsri3d_beam(8) qsri3d_modes(8)
qsri3d_modes_rigid(8) qsri3d_patch_s(8) qsri3d_patch_s_off(8) qsri3d_patch_t(8)
qsri3d_patch_t_off(8) reset_value_linear(4)
cmat_gate(4) gforce10(8) iface_3d(8) iface_3d_slip(8) iface_3d_stress(8)
iface_gen_quad9(9) iface_tangref(8) mesh_act_grav(4) qsri3d_beam_off(8)
```

(The last 9 are the solver-divergence subgroup; see below.)

## Results

| class | first pass | after this pass |
| --- | --- | --- |
| RECOVERED | 85 | **148** (85 + 63) |
| PARTIAL | 101 | **38** (9 + 29) |
| UNRECOVERABLE | 8 | **8** |
| total | 194 | 194 |

- **RECOVERED: 148** — faithful input reconstruction, runs `rc=0`, no lost
  value. The 63 new ones (the former `PARTIAL_run_ok_lost`) now also carry the
  derived `number_of_integration_points`; their fresh `element_dof` length
  matches the archived length in **63/63** cases.
- **PARTIAL: 38**
  - **9** former `PARTIAL_runfail_lost`: the lost value is derived and
    inserted, but the run still fails (Bi-CG/direct-LU divergence, see below).
  - **29** former `PARTIAL_runfail_other`: no lost value; all echo-able input
    records are present (verified by keyword diff against the archived `.dbs`),
    but the current binary does not reproduce the archived target.
- **UNRECOVERABLE: 8** — cannot produce a runnable `.dat` (unchanged).

## RECOVERED (148)

The 85 from the first pass, unchanged:

```
aeg_node aeg_seq bmom1 bmom_2d bmom_noint bmom_switch
bmom_truss cda1 cda_activate cdist_clamp cdist_corr cdist_normal
condif_convec condif_convec_el condif_heat_edge condif_heat_vol condif_heat_vol2 condif_rad
contact cpn1 cpn_angular cpn_geom cpn_sort cpn_zero
dbmeth dbmeth20 dignore dof1 dofid_no dofrhside
dpline1 dpline_eps dpline_geom dpline_group dpline_method dpline_move
dpline_n dpline_time dpoint1 dpoint_time dsmall dsmooth1
dsmooth_n fedge_alias fedge_restrict fproj_tunnel freq_timeint freq_timestep
fvol_elem genbeam1 groundflow_consolidate_off groundflow_flux_edge groundflow_nonsaturated_off groundflow_phreatic_multiple
groundflow_pressure_atm groundflow_pressure_atm_def groundflow_seepage groundflow_total_pressure_limit_dry groundflow_vangenuchten hreltime
iface_mc_1step meshdoff mfactor mlx1 mrepeat_save msf_axisym
msf_errors_2dwarn msf_nor msf_parse msf_print msf_shear mstrain_diprisco_elast
node_pressure node_slide partialname slide_axi smooth1 strain_settle
strain_settle_diag vtk_coord1 vtk_dof1 vtk_dofcalc1 vtk_empty1 vtk_nodmeth1
vtk_other1
```

The 63 new ones (derived `number_of_integration_points`):

```
cd_geom cd_method iface_3d_order iface_conv3d mat_rel mat_rel_reset
materi_direct materi_direct_auto materi_direct_mc materi_direct_visco materi_direct_wall mcap1
mcap1_comb mcap1_elast mcap2 mcap_legacy mcrunch mcrunch_low
mdirect_comp mdirect_gate mdp_shear mesh_act_grav2 mhardsoil_elast mhardsoil_elast2
mhardsoil_gp0 mhardsoil_gp0_off mhardsoil_plast mhardsoil_plast_elast mhardsoil_unload mhardsoil_unload_flat
mk0 mk0_off mmc_tension mmchs_soft mpower msf_beam2d_pure
msf_hex27_avg msf_parse_3d msf_quad9 msf_quad9_noavg msf_sheet3d_hex8 mshf
mshf_nof msph msph_flat mstrain_cap mstrain_compression mstrain_druckprag
mstrain_druckprag_elast mvoid mvoid_low myoung6 myoung6_apply myoung6_e2
myoung6_e3 qsri3d_beam qsri3d_modes qsri3d_modes_rigid qsri3d_patch_s qsri3d_patch_s_off
qsri3d_patch_t qsri3d_patch_t_off reset_value_linear
```

## PARTIAL (38)

### Solver divergence — lost value derived, run still fails (9)

The derived `number_of_integration_points` is inserted and the fresh
`element_dof` length matches the archive, but the current binary diverges
(Bi-CG breakdown, direct-LU retry, `|b-Ax|/|b|` not reduced). The archived run
converged (its `node_dof` satisfies the target). Cause: solver drift at HEAD
(`6eb399d` "correct Bi-CG direct-LU retry band width" and later physics
commits), not a reconstruction defect — the failures persist with the default
quadrature and with the derived value.

```
cmat_gate gforce10 iface_3d iface_3d_slip iface_3d_stress iface_gen_quad9
iface_tangref mesh_act_grav qsri3d_beam_off
```

### Solver divergence — no lost value (10)

Same Bi-CG/direct-LU breakdown signature; no `number_of_integration_points`
record exists for these, so there is nothing to derive.

```
contact_block contact_ctrl_apply contact_heatgen iface_gen iface_gen_geom
iface_gen_geom_off iface_gen_method iface_gen_method_gen iface_mc iface_mc_mem
```

### Converged but target missed — binary physics/feature drift (19)

The run completes, all echo-able input records are present (keyword diff
against the archived `.dbs` is empty apart from `end_*`), but the computed
target differs. Each is attributable to a post-2026-08-29 commit, not to a
missing record:

- `cda_arith`, `cda_copy` — `control_data_arithmetic` / `control_data_copy`
  now mutate the stored group value (archived keeps it unchanged).
- `bt_factor`, `iface_expansion`, `iface_mc_dil`, `iface_mc_dil_1step`,
  `iface_mc_gap`, `iface_mc_num`, `iface_mc_slip`, `iface_mc_slip_1step`,
  `iface_mc_tension`, `creset_iface` — interface/constitutive state differs
  (post-quadrature and Masin-kernel commits).
- `elmethod`, `groundflow_total_pressure_limit`, `numit`,
  `node_force_inertia`, `creset_geom` — computed DOF/rhs is 0 where the archive
  had a value (feature-path regression, e.g. `options_mesh -follow_material`,
  groundflow pressure limit).
- `groundflow_total_pressure_tension` — groundflow velocity sign fix
  (`0838457`).
- `gen1` — full mesh generation (`control_mesh_generate_truss`) produces no
  elements; the archived dump has `element`/`node`, the fresh run has none.

```
bt_factor cda_arith cda_copy creset_geom creset_iface elmethod gen1
groundflow_total_pressure_limit groundflow_total_pressure_tension iface_expansion
iface_mc_dil iface_mc_dil_1step iface_mc_gap iface_mc_num iface_mc_slip
iface_mc_slip_1step iface_mc_tension node_force_inertia numit
```

## UNRECOVERABLE (8)

Unchanged from the first pass; none of these has a `number_of_integration_points`
record to derive.

- **dbmeth21** — `.dbs` is a `-size_tot` dump (no input records); its input is dbmeth20 with `control_print_database_method 21 -sizetot`
- **dbmeth22** — `.dbs` is a `-size_tot_large` dump (no input records); sibling of dbmeth20/21
- **groundflow_interface** — current binary rejects `group_interface_materi_elasti_stiffness` record 10 (record layout changed)
- **iface_condif** — current binary rejects `group_interface_materi_elasti_stiffness` record 10 (record layout changed)
- **mc_pressure_min** — current binary requires `group_materi_plasti_camclay` length 4, dump has 3 (record layout changed)
- **mc_pressure_min_off** — current binary requires `group_materi_plasti_camclay` length 4, dump has 3 (record layout changed)
- **mdiprisco_hist** — current binary renamed diprisco history labels to `dipriscohis*`; input uses `-hisv10` (binary rename)
- **mstrain_diprisco** — current binary renamed diprisco history labels to `dipriscohis*`; input uses `-hisv10` (binary rename)

## Spot-check: run vs archived output

### Derivation evidence (element_dof length)

For every one of the 72 derived tests the fresh `element_dof` length equals the
archived length (72/72). Examples:

| test | element type | derived IP | archived element_dof | fresh element_dof |
| --- | --- | --- | --- | --- |
| mvoid | quad4 | 4 | 44 | 44 |
| msf_quad9 | quad9 | 9 | 144 | 144 |
| mhardsoil_elast | quad4 | 4 | 480 | 480 |
| msf_hex27_avg | hex27 | 27 | 2430 | 2430 |
| iface_3d | hex8 | 8 | 720 | 720 |
| cd_method | quad4 | 4 | 256 | 256 |

### Result reproduction (targets/post values)

Using a significance threshold `|x-y| > 1e-6 + 1e-4*max(|x|,|y|)` over
`node_dof*`, `post_point_dof`, `post_node_result`, `post_*_dof`:

- **27/63** new RECOVERED tests reproduce the archived result values exactly;
  the other **36** differ only by the documented binary drift.
- Baseline: of the 85 first-pass RECOVERED, **66** match exactly and **19**
  drift. So drift is expected even for byte-faithful inputs.

Concrete examples:

| test | record | archived | fresh | note |
| --- | --- | --- | --- | --- |
| gforce7 (known) | node_dof 2, comp 2 | -0.00157436754992 | -0.00157436754992 | exact |
| mvoid | node_dof 3, comp 5 | -0.576923076763 | -0.576923076763 | exact (matches target -0.5769±0.02) |
| msf_quad9 | node_dof 1, comp 6 | 0.243606537815 | 0.211018640117 | drift |
| cd_method | node_dof 3, comp 56 | 0.0957666496068 | -0.114352719345 | drift; target `post_point_dof sigxy 1±0.3` still met |
| iface_3d | node_dof 10, comp 6 | 0.00623029147559 | -4.6721071383e+14 | diverged (target 0±1) |

The `element_dof` length difference reported in the first pass as a **4.5x**
default-quadrature drift is now resolved for the 72 derived tests: with the
archived value restored the length matches. It remains for tests that never
set the record (their archived length encodes the *old* default), which is
expected and not recoverable from the input.

## Caveats

- The archived `.dbs` were produced on 2026-08-29; the current binary contains
  many later physics commits (`06e175d` per-operation quadrature, Masin kernel,
  node_dof parse fix `ddc9a5b`, groundflow velocity sign `0838457`, Bi-CG retry
  `6eb399d`, ...). Even a perfect reconstruction therefore need not reproduce
  the archived numbers, and the 36 drifting new RECOVERED tests are in the same
  situation as 19 of the 85 first-pass RECOVERED tests.
- Records the run overwrites with its final value cannot be recovered as input
  (`node_support_edge_normal_force`, ...). No test uses them; `node_dof` is
  handled through `node_dof_start_refined`.
- Inputs that are never echoed (`control_solver`, `include`,
  `support_edge_normal_element_node`, Abaqus/Gmsh/Feflow import records) are
  lost. `tslv_override` (known pair) is the only observed case.
- `dbmeth21`/`dbmeth22` are `-size_tot` dumps and contain no input.
- The derivation requires the `dof_label` record; all 72 targets and all 41
  known pairs have it. If a future dump lacks it, fall back to the K-probe
  cross-check (which needs a runnable reconstruction).

---

# Third pass (2026-09-14) — doc-grounded recovery and classification

Method: for each remaining test the matching `ProjectDocs/manual-user` /
`manual-developer` feature page was read; the reconstructed input was
compared with the documented **current** syntax and with the archived
`.dbs`. Acceptance uses the concrete values the docs give. All runs use
`build/tochnog` (HEAD `e983be8`), `LC_ALL=C`,
`LD_LIBRARY_PATH=<repo>/external-downloads/numlib-runtime`. Evidence:
`/tmp/opencode/docrecover-20260914-233512/` (`runs2/`, `evidence/`).

## 1. UNRECOVERABLE-8 — result

| test | doc used | reconstructed input | run |
| --- | --- | --- | --- |
| dbmeth21 | `control_print_database_method.md` | not a test: output of `dbmeth.dat` (`control_print_database_method 21 -sizetot`, base `dbmeth` → `dbmeth21.dbs`) | `dbmeth.dat` rc=0; `dbmeth21.dbs` 58 `Size of` + `Total size is` + system matrix |
| dbmeth22 | `control_print_database_method.md` | same `dbmeth.dat`, index 22 `-size_tot_large` | `dbmeth22.dbs` exactly 1 `Size of` (system matrix) |
| groundflow_interface | `group_interface_groundflow_*.md`, `group_interface.md` | input is faithful (no missing record) | **BLOCKED by bug** (see §4.1): `interface_element()` aborts |
| iface_condif | `group_interface_condif_conductivity.md` | input is faithful | **BLOCKED by bug** §4.1; adding a dummy `group_interface_materi_elasti_stiffness 10 1 1` makes it rc=0 |
| mc_pressure_min | `group_materi_elasti_camclay_pressure_min.md`, `materi_plasti_camclay_history.md` | 3-value camclay `group_materi_plasti_camclay 0 0.8 0.02 0.2` (drop the legacy 4th `N=3`) | **RECOVERED** rc=0; sigxx `-0.0093333`, sigyy `-0.2093333` (doc targets, exact) |
| mc_pressure_min_off | same | same 3-value fix | **RECOVERED** rc=0 |
| mdiprisco_hist | `materi_plasti_diprisco_history.md` | `-hisv10` → `-dipriscohis10` (documented Professional naming) | **RECOVERED** rc=0; dipriscohis10 `-152.79537` (target `-152.795±1`) |
| mstrain_diprisco | `materi_strain_plasti_diprisco.md` | same rename | **RECOVERED** rc=0; eppdipyy `-0.0105714`, eppdipxx `0.00350735`, dipriscohis10 `-152.79537` (targets met) |

Recovered `.dat` written to `validation-suite/test-2014/`:
`mc_pressure_min.dat`, `mc_pressure_min_off.dat`, `mdiprisco_hist.dat`,
`mstrain_diprisco.dat` (the 4-value camclay / `hisv` forms were the stale
archive-era syntax). `dbmeth21`/`dbmeth22` need no own `.dat`.

## 2. PARTIAL-38 — classification

Buckets: **R** reconstruction (input defect, fixed here), **B** bug
(doc says X, binary gives Y), **D** expected drift (archived target is
stale and the doc supports the new behaviour).

### 2.1 Reconstruction issues (fixed → now rc=0)

| test | defect | doc | fix / result |
| --- | --- | --- | --- |
| cda_arith | input `group_materi_elasti_young 0 2500` is the **post-run** value; the dump stores the mutated record | `manual-developer/control_data_arithmetic.md` ("young 1000 … final 2500") | initial 1000 → rc=0, target 2500 |
| cda_copy | initial young 500/1000 are the post-run values | `manual-developer/control_data_copy.md` ("young 1000 / 500 … final 500 / 1000") | initial 1000/500 → rc=0, targets 500/1000 |
| node_force_inertia | `node_inertia` is an **output** record but was seeded with the post-run value; `node_mass` was dropped | `node_force.md` ("declare node_inertia 0 0 … fills to -19.83") | `node_inertia 3 0 0` + `node_mass 3 0 2` → rc=0, -19.7515 (target -19.83±0.2) |
| numit | `-follow_material` reference geometry lost: the dump `node` is the **deformed** geometry (`node_start_refined` is the input reference) | `control_print_number_iterations.md` | primary feature OK (exactly 16 monitor lines); the secondary `-node_dof 2 -velix` target (`-0.123 = 0.877-1`) not reproducible — **unresolved** (see §3) |
| elmethod | same `-follow_material` reference loss | `control_print_element_method.md` | primary A/B OK (`-middle` 2 lines / `-node` 4 lines, `2*mid == x_node2`); secondary `-node_dof 3 -velix` (`-0.369`) not reproducible — **unresolved** (§3) |

### 2.2 Genuine bugs (doc contradicts the binary)

| test | doc | observed | source |
| --- | --- | --- | --- |
| creset_geom | `control_reset_family.md`: hisv0 `0.5 → 0` inside the brick (node 1 = 0, node 5 = 0.5) | node 1..4 = **-1**; the written value follows `3*value - 1` (verified: value 0.123 → -0.631, 1 → 2, 2 → 5). The 2026-08-28 binary gives the correct 0/0.5 | `data.cc` `data()` control_reset_dof value block (~1059-1228); suspect commit `64e43d2` |
| groundflow_total_pressure_limit | `groundflow_total_pressure_limit.md`: free nodes capped to 0.5 (else ~1.67) | free-node `pres = 0` — the pressure never develops, **with or without** the limit record | groundflow coupling / `0838457` sign path (follow-up) |
| creset_iface | `control_reset_family.md`: `element_interface_strain_normal ≈ 0` after `control_reset_interface_strain` | `-0.00050005` (target `0 ± 1e-6`) | `control_reset_interface_strain` timing/value (follow-up) |
| gen1 | `control_mesh_generate_truss.md`: generates trusses on the geometry | 0 elements / 0 nodes generated (archived: 2/2). The doc example needs an existing isoparametric element; the geometry-only path generates nothing | `generate.cc` mesh generation (follow-up) |

Plus the two bugs found while recovering the 8 (§4).

### 2.3 Expected drift (archived target stale, doc supports the new value)

| test(s) | doc / commit | archived | fresh |
| --- | --- | --- | --- |
| bt_factor | `integration_quadrature_and_stress_recovery.md` §"When a stored reference value changes"; commit `06e175d` (Lobatto → Gauss) | sigyy 20 | 34.641 (= 20·√3) |
| iface_expansion | `group_interface_materi_expansion_normal.md` + `group_interface.md` 2026-09-04 thermal fix | sigxx 1.984 | 0.998004 |
| iface_mc_dil, iface_mc_dil_1step | `group_interface.md` 2026-09-03 incremental-dilatancy fix | -2.8065 / -2.82675 | -1.39435 / -1.39338 |
| iface_mc_num, iface_mc_slip, iface_mc_slip_1step | `group_interface.md`: `ELEMENT_INTERFACE_FORCE_TANG` now stores the **elastic trial** (`interface.cc:711/942`), not the clamped force | 50 / 0 / 0 | 1e6 / 5000 / 5000 |
| iface_mc_gap | `group_interface.md` 2026-09-04 gap-condition inversion | velix 0 | -100.125 |
| iface_mc_tension | `group_interface.md` tension-limit fix | velix 3.88 | 100.102 |
| iface_mc, iface_mc_mem | `group_interface.md` (full-force assembly + solver honesty) | vely 0 | 3.59e9 |
| iface_3d, iface_3d_slip, iface_3d_stress, iface_tangref | `group_interface.md`, `control_print_interface_stress.md`; solver `solve_iterative_bicg.md` | converged | solver divergence |
| iface_gen, iface_gen_geom, iface_gen_geom_off, iface_gen_method, iface_gen_method_gen, iface_gen_quad9 | `control_mesh_generate_interface.md`; solver honesty | converged | solver divergence (1e13-1e14) |
| contact_block, contact_ctrl_apply, contact_heatgen | `control_contact_apply.md`, `contact_heat_generation.md`; solver honesty | converged | solver divergence |
| cmat_gate, gforce10, mesh_act_grav, qsri3d_beam_off | quadrature `06e175d` + solver honesty | converged | solver divergence / target drift |
| groundflow_total_pressure_tension | `groundflow_velocity_sign.md` (commit `0838457`) | -4.147 | -4.2896 |

The solver-divergence group is drift rather than a bug: `solve_iterative_bicg.md`
§1 documents that the old loop reported success on three non-convergence exits;
the current honest criteria report the true residual. The archived targets were
met only by that dishonest success.

## 3. Unresolved

- `numit`, `elmethod`: the `-follow_material` reference geometry
  (`node_start_refined`) is lost by the reconstruction (the dump `node` is
  the deformed geometry). Re-inserting `node_start_refined` as input did not
  reproduce the archived `-0.123 / -0.369` (neither the current nor the
  2026-08-28 binary); the primary features (16 monitor lines; 2/4 element
  lines) pass. Needs the original input form.
- `gen1`: geometry-only `control_mesh_generate_truss` generates nothing.

## 4. Bug list for follow-up (source-grounded, not fixed here)

1. **`interface.cc` `interface_element()`** — `db_len(GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS, element_group, VERSION_NORMAL)` aborts when the record is absent (non-elastic interfaces: `group_interface_condif_conductivity`, `group_interface_groundflow_*`). Blocks `iface_condif`, `groundflow_interface` (both doc-validated). Backtrace: `parallel_element_loop → elem → interface_element → db_max_index → db_error`. Fix: guard the `db_len` with `db_active_index(...)` (the line above already zeroes `ddum3` and the intent is a 3-value default).
2. **`database.cc:4456` `data_length[GROUP_MATERI_PLASTI_CAMCLAY] = 3`** — the developer manual (`materi_plasti_camclay_history.md`: "variable-length 3..4") and the user manual ("the legacy … four values … still accepted") document the 4-value form, but the parse PUT rejects length 4 (`Length too small … should become at least 4`). `plasti.cc` still reads `plasti_data[3]` for `ldum>=4`, so the legacy path is unreachable. Test is recoverable with 3 values, but the doc/binary mismatch is real. **FIXED in source** (fourth pass, §1: `data_length = 4`).
3. **`data.cc` `data()` control_reset_dof value block** — the reset value is written as `3*value - 1` instead of `value` (creset_geom, regression after 2026-08-28; suspect `64e43d2`).
4. **`control_reset_interface_strain`** — leaves `element_interface_strain_normal = -5.0e-4` instead of ≈0 (creset_iface).
5. **groundflow pressure / `groundflow_total_pressure_limit`** — the free-node pressure does not develop (0 instead of 0.5/~1.67) with a faithful input. **Reclassified (fourth pass, §2): reconstruction defect, not a binary bug** — the input was missing the `groundflow_consolidation_apply -yes` + `inertia_apply -yes` records (defaults changed in `0560dcf`).
6. **`control_mesh_generate_truss`** — geometry-only generation produces no elements (gen1). **Reclassified (fourth pass, §3): reconstruction defect, not a binary bug** — the base `node`/`element 1 -bar2` (an INPUT mesh, not generated) was dropped, so `NODE_NODE` was empty.

## 5. Updated counts

| class | before 3rd pass | after 3rd pass |
| --- | --- | --- |
| RECOVERED | 148 | **155** (+4 from the 8: mc_pressure_min, mc_pressure_min_off, mdiprisco_hist, mstrain_diprisco; +3 from the 38: cda_arith, cda_copy, node_force_inertia) |
| PARTIAL | 38 | **35** (2 unresolved reconstruction + 4 bug + 29 drift) |
| UNRECOVERABLE | 8 | **4** (groundflow_interface, iface_condif blocked by bug #1; dbmeth21/dbmeth22 resolved as outputs of `dbmeth`, no own `.dat`) |

Classification of the 38: 5 reconstruction (3 fixed + 2 unresolved), 4 bug, 29 drift.

> Note: `dbmeth21`/`dbmeth22` were never separate tests; they are the
> `-sizetot` / `-size_tot_large` dumps of the `dbmeth` test. The
> `UNRECOVERABLE (8)` section above predates this finding.

---

# Fourth pass (2026-09-15) — materi/mesh/groundflow follow-up

Three items of the §4 bug list were investigated on a scratch worktree from
`5c80b60` (branch `dev/materi-mesh-groundflow-fixes`). Evidence:
`/tmp/opencode/matmesh-20260915-074922/evidence/`.

## 1. `group_materi_plasti_camclay` — REAL bug, fixed in source

Cause: `database.cc:4456` set `data_length[GROUP_MATERI_PLASTI_CAMCLAY] = 3`
(maximum stride 3) with `fixed_length = 0`. The parser PUTs the record
length it read, so a 4-value record aborted at parse time (`Length too small
of : group_materi_plasti_camclay … should become at least: 4`) and the
legacy GNU branch in `plasti.cc:156` (`if (ldum>=4) N = plasti_data[3]`) was
unreachable.

Fix: `data_length = 4` (maximum stride; `fixed_length` stays 0, so the
record is variable-length 3..4). A 3-value record keeps the `DBL_MAX`
sentinel in slot 3 and `db_len()` returns 3; a 4-value record returns 4.

Validation (fixed binary):

- `mc_pressure_min`, `mc_pressure_min_off` (3-value): rc=0, sigxx
  `-0.0093333`, sigyy `-0.2093333` (documented targets, exact).
- `mc_pressure_min` with a 4th value `N=3`: rc=0, same sigxx/sigyy.
- The legacy closure is now reachable: with the p0 history zeroed, the
  4-value form derives p0 from N and runs rc=0; the 3-value form correctly
  errors "needs initial p0".

## 2. `groundflow_total_pressure_limit` — reconstruction defect, not a bug

Cause: the reconstructed test relied on the pre-`0560dcf` defaults
(`groundflow_consolidation_apply -yes`, `inertia_apply -yes`). Since
`0560dcf` aligned both with the Professional (manual 6.557: "Default switch
is -no"), a faithful input without the records gives pres=0: there is no
material-divergence source and no capacity term, so the limit record has
nothing to clamp.

Fix (test input): add `groundflow_consolidation_apply -yes` and
`inertia_apply -yes`.

Validation (fixed binary):

- `groundflow_total_pressure_limit` (limit 0.5): rc=0, free nodes pres=0.5,
  Dirichlet nodes 0.
- Without the limit record: free nodes pres=1.6667 (the ~1.67 of the
  manual), target 0.5 correctly missed.
- `groundflow_total_pressure_limit_dry` (limit 0): rc=0, pres=0 (the dry
  check skips the consolidation term even with both switches active).

The developer manual `groundflow_total_pressure_limit.md` was corrected
accordingly.

## 3. `control_mesh_generate_truss` (gen1) — reconstruction defect, not a bug

Cause: `generate_beam_truss()` (`generate.cc:491`) generates trusses between
nodes connected by an existing isoparametric element (it reads `NODE_NODE`
and `NODE_START_REFINED`). The reconstruction dropped the whole `node` /
`element` / `element_group` set for this test, including the base
`element 1 -bar2 1 2` (an INPUT element, not generated), so `NODE_NODE` was
empty and nothing could be generated. The manual's own example shows the
required base mesh.

Fix (test input): restore `node 1 0`, `node 2 1`, `element 1 -bar2 1 2`.

Validation (fixed binary): rc=0; `element 2 -truss 1 2` generated,
`node_dof 2 disx = 0.5` (target 0.5).

## Corpus A/B (bug 1 only)

The source change is the `data_length` bump. Full sfnet 197 corpus with the
fixed binary vs the stored baseline: **187 PASS / 10 RUNFAIL / 0 PARSE**,
classification byte-identical to the baseline (no per-test change). Suite
build_safe **16/16**.
