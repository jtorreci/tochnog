# Keyword gap vs Tochnog Professional (2026-09-26)

## Why this report exists

The GNU→Professional convergence cycle closed with **behavioural** convergence:
every test in the sfnet corpus, the reconstructed validation-suite (233/233)
and the Professional validation files runs and matches the Professional
binary. That claim was never about **input-language coverage**. A keyword that
no test uses can be missing and the corpus still passes.

Example: `geometry_list` (manual 6.531, plus section 9.5 "Use -geometry_list
as geometry entity") was not recognised by the GNU parser (implemented
2026-09-26, see update below).

## Method (reproducible)

1. Extract every keyword heading from chapters 4, 5 and 6 of
   `external-downloads/professional/UserManual-professional.txt`.
   The text uses TeX OT1 ligature bytes (0x1B ff, 0x1C fi, 0x1D fl, 0x1E ffi,
   0x1F ffl); they are expanded first, otherwise names such as
   `print_filter` or `node_stiffness` are truncated.
2. Feed each name to the real `db_number()` of the current object files
   (`scripts/keyword_probe.cc`, linked against every `*.o` except `tn.o`).
   This honours all aliases resolved inside `db_number` (for example
   `force_edge*` → `force_element_edge*`).
3. Names rejected by `db_number` but handled as literals in `input.cc`
   (initialisation part: `derivatives`, `number_of_space_dimensions`,
   `materi_acceleration`, `echo`, `end_initia`, `end_data`) count as present.

Build probe:

```sh
g++ -c scripts/keyword_probe.cc -o /tmp/probe.o
g++ -o /tmp/probe /tmp/probe.o $(ls *.o | grep -v '^tn\.o$') \
  <same LINK_FLAGS_AFTER as scripts/build_safe.sh, without -o>
/tmp/probe < pro_keywords.txt
```

Limitation: "accepted by the parser" does not imply "implemented". Some
keywords are registered parse-only (for example `control_data_save`). This
report measures the lower bound of the gap.

## Result

| Metric | Value |
|--------|-------|
| Professional keywords (chapters 4–6) | 1154 |
| Accepted by GNU parser | 832 (72.1 %) |
| Missing | 322 |
| Missing but already discarded (GiD) | 63 |
| Missing and not discarded | **259** |
| Coverage excluding the GiD discard | 76.3 % |

The "Métrica global" header of `SEGUIMIENTO-CONVERGENCIA.md` (309/1094, 28 %)
is stale and should be replaced with these figures.

Full list: `ProjectDocs/keyword-gap-2026-09-26.txt`.

Update 2026-09-26: `geometry_list` implemented; the figures above and the
full list already reflect it (826 → 827 accepted, 265 → 264 pending).

Update 2026-09-27: `materi_damage_apply`, `materi_failure_apply`,
`plasti_apply` and `control_plasti_apply` implemented (verified against the
Professional binary 02-08-2026; see `manual-user/<name>.md`). Probe re-run on
the clean build: 827 → 831 accepted, 264 → 260 pending; the four names were
removed from the full list.

Update 2026-09-30: `element_normal` implemented (generated record used by
`*_direct_normal_automatic`, rule measured on the Professional 02-08-2026;
see `manual-user/element_normal.md`). Probe re-run on the clean build:
831 → 832 accepted, 260 → 259 pending; the name was removed from the full
list.

Update 2026-11-01 (F5 start): the keyword list was re-extracted from the
manual shipped with the Professional binary 02-08-2026 (`doc/user.pdf`,
dated August 2, 2026, `pdftotext -layout`, every `4.x` / `5.x` / `6.x`
heading of chapters 4-6: 1156 names) and probed on the clean build of
`7da9e21`. Accepted 824 + the 6 literal names of the initialisation part;
missing 320 = 63 GiD (discarded) + **257** pending. Differences with the
2024 manual: the `*_node_method` records of `print_vtk`, `print_gmsh`,
`print_gid` and `control_print_*` are renamed `*_node_type` (plus the new
`control_print_dof_line_node_type`), `global_post_point_node_type` becomes
`post_point_node_type`, and `mesh_delete_geometry_moving` is gone. Full
list: `ProjectDocs/keyword-gap-2026-11-01.txt` (supersedes the 2026-09-26
list for F5).

Update 2026-11-02 (end of F5 wave W1, `c135eeb`): the 2026-11-01 list
probed again on the clean build: 34 of the 320 missing names are now
accepted (the W1 switches and control records, F5 T1 U1-U6); accepted
858 + 6 literals, missing **286** = 63 GiD (discarded) + **223** pending.
`control_slide_damping_apply` stays missing (moved to T5 with
`slide_damping`).

Update 2026-11-07 (end of F5 waves W2 output records and W3 materials,
`f427646`): the 2026-11-01 list probed again on the clean build
(`scripts/keyword_probe.cc` linked with the build objects): 131 of the 320
names are accepted (97 more than after W1); accepted 955 + 6 literals,
missing **189** = 63 GiD (discarded) + **126** pending. The material family
is complete: the two names still reported, `materi_history_variable` and
`materi_strain_energy`, are initialisation names normalised by `input.cc`
(they are not `db_number` records).

Update 2026-11-13 (end of F5 wave W5 geometry, mesh tools, mpc, dependency, `16f3721`): the
2026-11-01 list probed again on the build of the same sources (`scripts/keyword_probe.cc` linked
with the build objects): 199 of the 320 names are accepted (68 more than after W2/W3); accepted
1023 + 6 literals, missing **121** = 63 GiD (discarded) + **58** pending. The 58 are the F5 T5
other-physics names (baseline correction, force_edge_diagram, force_gravity_geometry,
`force_volume_element_group_0`, `groundflow_phreatic_*`, `group_beam_*`, `group_contact_spring_*`,
`group_dof_initial*`, `group_time_fill`, `group_truss_elasti_elongation_force_diagram`,
`slide_damping`, `control_slide_damping_apply`, `slide_user`, `strain_volume_*`), the 31
`safety_*` names of T6, and the two initialisation names `materi_history_variable` /
`materi_strain_energy` normalised by `input.cc`. None of the 58 is a T4 (W5) name.
2026-11-14 (F5 T5, owner decision): `force_volume_element_group_0` is discarded, it is a PDF header
artifact of manual 6.504 and not a record name (`force_volume_element_group` is accepted); 57 pending.
2026-11-15 (F5 T5 U3b): `groundflow_phreatic_project` accepted.
2026-11-17 (F5 T5 U2): `slide_damping`, `control_slide_damping_apply`, `slide_user` accepted.
2026-11-17 (F5 T5 U4): `group_dof_initial`, `group_dof_initial_specific_number`, `group_dof_initial_specific_value` accepted. The T4 tracker counted 69 W5 names, the probe shows 68 newly
accepted; the one-name difference in the 2026-11-08 count is not explained.

2026-11-18 (F5 T5 U5): `force_edge_diagram` (6.455) and `force_gravity_geometry` (6.499) accepted.
2026-11-19 (F5 T5 U6): `group_contact_spring_direction_automatic_planes` (6.606) and
`group_truss_elasti_elongation_force_diagram` (6.773) accepted.
2026-11-19 (F5 T5 U7): `group_time_fill` (6.770) accepted.
2026-11-22 (F5 T5 U8): `strain_volume_*` (6.1062-6.1066) computed, `post_strain_volume_absolute` / `_relative` filled.
2026-11-22 (F5 T5 U9): `bounda_baseline_correction`, `bounda_baseline_correction_parameters` accepted.
2026-11-22 (F5 T5 U10): `group_beam_force_moment_plasti` accepted (plastic hinges).
2026-11-22 (F5 T6 U1): `control_safety_slip`, `safety_slip_multi_linear*` (4) accepted.
2026-11-22 (F5 T6 U2): `safety_slip_combined_linear*` (4) accepted.
2026-11-22 (F5 T6 U3): `safety_slip_circle_grid_*`, `safety_slip_circle_line_*` (12) accepted.
2026-11-22 (F5 T6 U4): `safety_slip_set`, `safety_slip_set_result` accepted; control_safety_slip N passes.
2026-11-22 (F5 T6 U5): `control_print_gid_safety_slip_critical` accepted (GiD plot of the slip surfaces).
2026-11-22 (F5 T6 U6): `safety_slip_ellipsoide*` (5), `safety_slip_grd*` (4, plus the output `safety_slip_grd_result`) accepted (3D slip surfaces, GiD triangle mesh).

Update 2026-11-24 (end of F5, build `231f1da`): full re-inventory. The names were extracted again
from `doc/user.pdf` of the Professional 02-08-2026 (`pdftotext -layout`), this time from the table
of contents and the body headings together, keeping the longest name per section number: the body
layout cuts long names (`..._cam_clay_for_bonded_` in 6.700) and the table of contents glues the page
number to them (`..._history104`). The 2026-11-01 list was taken from the body alone and lost the
three bonded-soil Cam clay names that way, so they were never in the F5 scope. List:
`ProjectDocs/pro-keywords-2026-11-24.txt` (section, name; 1147 distinct names). Probe result:

| Metric | Value |
|--------|-------|
| Professional keyword names (chapters 4-6) | 1147 |
| Accepted by `db_number` | 1073 |
| Initialisation literals (`derivatives`, `echo`, `end_data`, `end_initia`, `materi_acceleration`, `number_of_space_dimensions`) | 6 |
| Initialisation names normalised by `input.cc` (`materi_history_variable`, `materi_strain_energy`) | 2 |
| Present in the GNU | **1081 (94.2 %; 99.7 % excluding the discarded names)** |
| Discarded (62 GiD names, `force_volume_element_group_0`) | 63 |
| Missing | **3** |

The 3 missing names are the generalised non-associate Cam clay for bonded soils:
`group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils` (6.700),
`materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils_history` (4.21) and
`materi_strain_plasti_generalised_non_associate_cam_clay_for_bonded_soils` (4.38). The manual gives
no parameters (6.700 reads "index ...") and one theory line ("a modification of the 'Milan' model
of Prof. Roberto Nova"). The Professional binary registers the three names. One-element probe on
the Pro 02-08-2026 (the Cam clay test `materi_strain_plasti_model3` with this group record, plain
elasticity and both initialisation names): fewer than 12 values give "Error detected for data item";
12 to 16 values, and 17 or 20 increasing values, give "Program error, assertion failed"; 17 or more
equal values give "Error found in : Xd data ... The value ... is not larger than ...". No probe ran
to a result, so the parameters and the behaviour cannot be measured on the Professional.

`SEGUIMIENTO-CONVERGENCIA.md`: the 98 checklist entries still marked PENDIENTE were classified with
the same probe (31 accepted, 46 GiD, 11 OCR fragments of the 2024 manual, 5 literals, 1 renamed,
1 removed from the 2026 manual, 3 bonded-soil names still missing). The family counts in the
section headings were not recomputed.

2026-11-24 (bonded-soil Cam clay): the 3 missing names are accepted and implemented from the Nova
bonded-soil model ([user page](manual-user/group_materi_plasti_generalised_non_associate_cam_clay_for_bonded_soils.md));
1084 of 1147 names present, 0 missing apart from the 63 discarded.

2026-11-24 (GiD family reopened, odd/tasks/gid.md): the GiD output and geometry names are accepted and
implemented; probe on the build of `8a53ace`: 1137 accepted + 6 literals + 2 normalised = 1145 of 1147;
discarded `control_print_gid_batch` and `force_volume_element_group_0`.

## Gap by family

| Family | Missing | Notes |
|--------|---------|-------|
| GiD (`control_print_gid_*`, `print_gid_*`, `mesh_gid_*`, `control_mesh_gid_batch`) | 63 | Discarded. Note: `mesh_gid_assign_conditions_*` is the Professional route to fill `geometry_list`. |
| `control_mesh_*` | 34 | extrude (5), truss_distribute_mpc (7), remove/remove_really (5), delete_geometry options (4), refine_locally options (3), convert quad9→quad6 / tria6→tria3, split_element_from/to, renumber offsets, duplicate_element_group, map |
| `safety_slip*` + `control_safety_slip` | 32 | Whole slip-surface safety-factor feature (circle grid/line, ellipsoide, multi_linear, combined_linear, grd, set). |
| Material groups (`group_materi_*`, `materi_plasti_*`) | 29 | umat (3), stress_null (3), hypo niemunis_visco / wolffersdorff_niemunis / minimum_void_ratio, generalised cam clay for bonded soils, hyper reduced_polynomial / simo_taylor, elasti borja_tamagnini / c / young_user, history_variable_user, visco_power, residual_stiffness, mpc, diprisco_density, `materi_plasti_maximum_iterations` (only the `group_` form exists) |
| Geometry | 22 | `geometry_moving*` (6) + `control_geometry_moving` + `mesh_delete_geometry_moving`, `geometry_hexahedral`, `geometry_tetrahedral`, `geometry_boundary`, `geometry_exclude`, `geometry_method`, `geometry_*_part`, `*_eps_iso`, `geometry_element_geometry*`, `geometry_mpc` |
| Post (`post_*`, `control_post*`) | 21 | post_calcul options (7), post_data (3), post_count, post_node/point options |
| Element records (`element_*`) | 19 | mostly print/output records: element_intpnt_* (5), element_interface_intpnt_* (4), truss strain, beam direction/force, geometry_present, boundary (element_normal done 2026-09-30) |
| `*_apply` switches | 14 | interface_gap_apply, linear_calculation_apply, element_group_apply, truss_rope_apply, group_plasti_apply, ... and their `control_*` forms (plasti_apply, control_plasti_apply, materi_damage_apply, materi_failure_apply done 2026-09-27) |
| Print (`print_*`) | 14 | print_vtk_* (4), print_gmsh_* (3), print_frd_freecad/prepomax, print_precision, print_data_name, print_element_geometry_present* |
| Initialisation dofs (`materi_*`) | 11 | materi_strain_plasti_* (5), materi_strain_total_*kappa (4), materi_strain_energy, materi_history_variable |
| MPC | 8 | mpc_element_group_* (5), mpc_geometry_method/tolerance, control_mpc_element_group |
| `control_timestep_*` | 5 | until_data/minimum/maximum, adjust_minimum_iterations, iterations_automatic_minimum_maximum_wished |
| groundflow | 4 | phreatic_bounda / level_n / only / project (+ groundflow_apply above) |
| dependency | 4 | dependency_geometry/method/number/type |
| Other | 42 | bounda_baseline_correction (2), group_dof_initial* (3), group_beam_* (3), group_contact_spring_* (3), strain_volume_* (3), processors_maximum/partition, slide_damping/user, mesh_boundary/correct, control_zip, control_system_call, control_input, message, mrange, mstring, ... |

## Suggested next cycle

Convergence should now be tracked with two numbers, not one:
behavioural (corpus pass rate) and language (probe coverage). The probe should
run in CI so a regression in keyword registration is caught.

Priority proposal (value for geotechnical users first, cost second):

1. ~~`geometry_list` and `-geometry_list` as geometry entity~~ — DONE
   2026-09-26 (see `manual-user/geometry_list.md`).
2. `*_apply` switches and `control_timestep_until_*` — mostly cheap switches.
3. Material groups actually used by Professional examples (umat, stress_null,
   hypo niemunis variants, bonded-soil cam clay).
4. `geometry_moving` family and `control_mesh_*` topology tools.
5. `safety_slip` family — large, self-contained feature; plan it separately.
6. Output records (`element_*`, `print_*`, `post_*`) — low risk, many items.

Before implementing each keyword, follow the AGENTS.md rule: find a usage
example (manual, Professional repo, sfnet) and record the verification.
