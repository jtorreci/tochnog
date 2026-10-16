# group_integration_method -reduced (developer notes)

F4 T6f C80, 2026-10-16. User page:
[group_integration_method_reduced](../manual-user/group_integration_method_reduced.md).

## Where it is implemented

| file | function | role |
|---|---|---|
| `tochnog.h`, `tochnog-mod.h` | enum | `REDUCED`, `GROUP_INTEGRATION_METHOD_REDUCED_FACTOR` appended after `UPDATED_JAUMANN` (both headers in sync); `stress_recovery_reduced_iso`, `integration_reduced_factor()` declared |
| `database.cc` | `db_initialize` | name `"reduced"`; record `group_integration_method_reduced_factor`: DOUBLE_PRECISION, length 1, `version_all` 1, required `GROUP_TYPE` (like `group_integration_method`) |
| `polynom.cc` | `integration_reduced_factor()` | reads the factor, default 10 (manual 6.621), `db_error` for f <= 0 (or NaN) |
| `polynom.cc` | `pol()` tensor branch | the rule table below; sets `stress_recovery_gauss`, `stress_recovery_reduced_iso`; on-axis Gauss switch |
| `miscel.cc` | `recovery_fit_shape()`, `recovery_fit_factor()` | the kind 1 fit places the quad4/hex8 points at +-1/f; the cache key includes `riso` |
| `calcul_force.cc` | `msf_element_rule()` | replica of the rule and of the on-axis switch (section forces) |
| `delete.cc` | `mc_element_rule()` | replica of the rule (mesh cut forces; no axisymmetric groups) |

## Rule (Professional 02-08-2026, 135 probes, residual <= 2e-13)

In `pol()` the per-direction loop has a new branch before the generic
fallback:

- `integration_method==-REDUCED && npol==2 && npoint_per_dir==2 &&
  (QUAD4 || HEX8)`: iso = +-1/f, weight_iso = 1/2 each (the GNU 1D weights
  sum to 1; the Professional has 1 each on [-1,1]). `stress_recovery_gauss
  = 1`, `rule_reduced = 1`.
- Every other `-reduced` case falls to the pre-existing generic branches:
  fewer points than nodes -> Gauss, otherwise Lobatto. This gives bar2
  default 1 point = Gauss, bar2 `-maximal` = Lobatto, quad4/hex8
  `-minimal` = 1 point, bar3/quad9/hex27 default = Lobatto -1,0,1, and
  `-minimal` = 2-point Gauss. These are the measured Professional rules.
  tria3/tet4/prism6 never read the method. bar4/quad16/hex64 take the
  same fallback (Lobatto, `-minimal` Gauss); not measured.
- After the loop: for `-reduced` quad9/hex27 with the Lobatto rule,
  `stress_recovery_gauss = 1`. `recovery_fit_kind()` then returns kind 1
  (physical linear fit at the Lobatto points), not kind 3 (point value of
  `-lobatto`). bar2/bar3 keep the `-lobatto` / Gauss flags (C67).
- `stress_recovery_reduced_iso = 1/f` only when `rule_reduced`.

Nodal recovery: the Professional uses kind 1 (physical linear fit plus the
arithmetic mean) for every quad/hex `-reduced` rule, including f = 1 with
the points at the nodes. The probes: `-lobatto` and `-reduced` f = 1 give
the same velocities but different nodal stresses.

## Axisymmetric on-axis switch (C71 extension)

Professional probe (axisymmetric quad4 column touching r = 0, f = 1 and
f = 2): the element with a node on the axis gets the 2x2 Gauss points
(iso +-0.57735) for both factors. The off-axis elements keep +-1/f. So
the Professional switches on the NODE radius, not on the point radius.
`pol()` checks `|new_coord[inol*ndim]| < EPS_COORD` for every node when
`rule_reduced` and the group is axisymmetric, then sets
`integration_method = -GAUSS` and goes back to `build_rule`. The quad9 /
hex27 `-reduced` Lobatto rule is caught by the existing C71 point check
(its points are at the nodes). `msf_element_rule()` repeats the node
check. Test `reduced_quad4_axisym1` (f = 2).

## Default factor

`integration_reduced_factor()` returns 10 when the record is absent
(manual 6.621). The Professional 02-08-2026 uses 1: without the record its
results equal the f = 1 run to all digits. This is a deliberate class A
divergence (owner 2026-10-16; DIVERGENCES C80, SEGUIMIENTO Descarte). Test
`reduced_default_factor1`: its targets are the Professional node_dof of
the explicit f = 10 run.

## Verification

- Tests (Professional targets 4e-6): `reduced_quad4_f2_1`,
  `reduced_quad4_f10_1`, `reduced_hex8_f2_1`, `reduced_quad9_minimal1`,
  `reduced_quad9_maximal1`, `reduced_quad4_axisym1`, and
  `reduced_default_factor1` (GNU default). On HEAD `0ac600d` all give
  rc 1 ("I do not know : -reduced"). GREEN: every node_dof equals the
  Professional within 5e-12 (relative, column floor 1e-6 of the global
  maximum).
- The earlier probe set (`gen.py`: bar2, bar3, quad4, quad9, hex8, hex27,
  tria3, tet4, prism6 x f = 2, 4 x default/min/max) matches within 5e-12,
  except two groups. bar3: 1.3-2.3 %, the stress recovery, C67.
  prism6 `-minimal`: the known C76.

- Clean build; `build_safe.sh` TN_CHECK_ONLY 279 OK / 0 FALLO, runs
  226/226; `audit_db_buffers.py` exit 0. Corpus A/B on omen (966 runs,
  baseline the C77 build): no run changed (masked `.dbs` identical); no
  corpus input uses `-reduced`; the 7 new tests rc 0.

## Hardcoded data and pending refactorings

- Default factor 10 and the threshold `EPS_COORD` (1e-9) of the axis
  switch are hardcoded.
- The rule is duplicated in `pol()`, `msf_element_rule()` and
  `mc_element_rule()`. A shared per-direction rule helper would remove
  three copies of the method branching (also for `-lobatto` and SRI).
- The integer method values (`-GAUSS`, `-LOBATTO`, `-REDUCED`) are the
  generic name enums, not a dedicated method enum.
- Not done: the bar3 `-reduced` linear-fit recovery (C67), the
  Professional fallback at very large f (element average), and the
  selective / B-bar extension (BACKLOG C80.1).
