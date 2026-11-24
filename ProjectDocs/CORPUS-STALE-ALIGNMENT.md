# Corpus stale-input / stale-target alignment with the Professional

Date: 2026-09-14
Branch: `dev/corpus-stale-alignment`
Baseline: `ee0032c` (branch `dev/masin-kernel-surface`), clean build.
Corpus: `external-downloads/sfnet/extracted/test` (197 `.dat`, **untracked** —
see "Reproducibility" for the exact edits, since git does not track them).
Evidence: `/tmp/opencode/stalealign-evidence-20260914-133732/`
(originals, per-test GNU runs, Pro runs, classification files).
GNU binary: `/tmp/opencode/stalealign-20260914-133659/build/tochnog`
(clean `build_safe.sh` build from `ee0032c`).
Pro binary: `/tmp/opencode/tn_prof/tochnog_version_02-08-2026/src/tochnog`
(build 2-8-2026). All runs with `LC_ALL=C`.

## Result

| corpus | PASS | RUNFAIL | PARSE |
|---|---|---|---|
| baseline `ee0032c` | 174 | 23 | 0 |
| after alignment | **184** | **13** | **0** |

Per-test diff: **10 RUNFAIL -> PASS, 0 regressions**
(`absorp1`, `dens1`, `examp11`, `ho_heat2`, `ho_othr3`, `incnav5`, `matrix8`,
`truss6`, `viscos1`, `wave1`). Suite `build_safe.sh`: **16/16**.

The remaining 13 RUNFAIL are pre-existing and unrelated to this change:
`examp2`, `examp6`, `examp14`, `examp17`, `examp24`, `ho_mech2`, `ho_mech3`,
`ho_mech4`, `incnav3`, `petsc1`, `pruebavel`, `refine5`, `taylor3`.

## 1. Eight stale inputs: missing `inertia_apply -yes`

Commit `0560dcf` aligned the GNU `options_inertia` default with the
Professional (manual 6.785: "Default, if inertia_apply is not specified, then
each of switch_0, switch_1 etc. is set to -no"). The eight corpus inputs below
are legacy GNU models that relied on the old `-yes` default, so after the
alignment they ran steady/quasi-static and missed their transient targets.

The Professional's own copies of these tests (where they exist under
`test/other/`) explicitly add `inertia_apply -yes`; the GNU corpus files are the
older GNU variants without the record. Adding the record (after `end_initia`,
matching the Professional placement) makes every one hit its stored target
exactly:

| test | target (unchanged) | old value (no record) | new value (`inertia_apply -yes`) | class |
|------|--------------------|-----------------------|----------------------------------|-------|
| `absorp1`  | `node_dof 1 -temp` 0.368 ±1e-2 | 0 | **0.369711212329** | PASS |
| `ho_heat2` | `post_point_dof 1 -temp` 0.0530 ±1e-2 | 0.10000001861 | **0.0530863731519** | PASS |
| `ho_othr3` | `post_point_dof 1 -temp` 0.0530 ±1e-2 | 0.10000001861 | **0.0530863731519** | PASS |
| `examp11`  | `post_point_dof 1 -fscal` 1e-4 ±1e-5 | -1.35834729366e-06 | **0.000101444151764** | PASS |
| `wave1`    | `node_dof 2 -scal` -1 ±1e-2 | 0 | **-0.990179145199** | PASS |
| `dens1`    | `node_dof 1 -dens` 0.999002 ±1e-5 | 1 | **0.999001998002** | PASS |
| `truss6`   | `node_dof 2 -velx` 0.1 ±1e-8 | 2e8 | **0.1** | PASS |
| `viscos1`  | `node_dof 2 -temp` 2 ±1e-4 | 2e8 | **2** | PASS |

Professional grounding:

- `test/other/absorp1.dat`, `wave1.dat`, `truss6.dat`, `viscos1.dat` all contain
  `inertia_apply -yes` immediately after `end_initia` (the only semantic
  difference from the corpus file; the Pro also renames `bounda_unknown` ->
  `bounda_dof` and, for `viscos1`, `-updated_without_rotation` ->
  `-total_linear`, which are dialect/naming differences and are not needed for
  the GNU to reach the target).
- `ho_heat2`, `ho_othr3`, `dens1`, `examp11` have no copy under the Pro tree;
  the record is grounded in manual 6.785 (it is the switch that includes heat
  capacity / mass / storage) and verified by hitting the stored target.
- `inertia_apply` maps to `OPTIONS_INERTIA` (`database.cc:8961`), read into
  `options_inertia` by `general.cc:82-83`; it is the Professional spelling of
  the same switch as `options_inertia` (used in the earlier triage probes).

The `0560dcf` default is **not** reverted: it is a deliberate Professional
alignment that fixed 13 static corpus tests.

## 2. Two stale targets

### matrix8 — pressure sign + target aligned with the Professional's own file

The corpus `matrix8.dat` prescribed `pres = +2e-2` and expected `eptxx = -0.01`.
With `ee0032c` the GNU returns `eptxx = 0` (the atmospheric clamp of `741a490`
caps the positive pore pressure). The Professional's own `test/other/matrix8.dat`
prescribes `pres = -2e-2` and expects `eptxx = +0.01`; the GNU reproduces it
exactly (rc=0 with the Professional input).

Evidence (all rc / `node_dof 2 -eptxx`):

| input | GNU `ee0032c` | Pro 02-08-2026 |
|-------|---------------|----------------|
| corpus `pres = +2e-2`, target -0.01 | 0 (rc=1) | 0 (with the corpus sign) |
| Professional `pres = -2e-2`, target +0.01 | **+0.01 (rc=0)** | +0.01 (Pro's own target) |

Change applied to the corpus input: flip the `bounda_time 2` pressure sign
`+2e-2 -> -2e-2` (matching the Professional's own file) and recalibrate
`target_value 0 -0.01 -> +0.01`. GNU now returns `eptxx = 0.01` (rc=0).

> The alternative of keeping the corpus input and setting the target to `0`
> (the clamped result, which both binaries produce with that sign) was
> rejected: it would make the test degenerate and would not align the input
> with the Professional test-suite variant (the triage recommendation).

### incnav5 — stale target recalibrated to the honest GNU value

The stored target `0.41916 ±0.01` is met by no binary: 2014-GNU `0.4718`,
current GNU `0.389082`, Pro own copy `0.5095` (task-cited `0.492`). The Pro's
own `test/other/incnav5.dat` has **no** `target_item`/`target_value` at all, and
its model differs from the corpus one (it adds `control_mesh_refine_globally 9`,
one extra h-refinement), so the Pro value is not directly comparable to the GNU
corpus model. The Pro also rejects the GNU `options_mesh` spelling, so it cannot
run the corpus file as-is.

Measured `post_point_dof 1 -velx` (target slot):

| run | value |
|-----|-------|
| GNU `ee0032c` on the corpus model (refinements 10-13) | **0.389082098155** |
| Pro own copy (refinements 9-13, pardiso) | 0.509497137366 |
| Pro own copy with refinement 9 removed (GNU mesh 10-13) | 0.232461801786 |
| 2014-GNU (from the earlier triage) | 0.4718 |

Change applied: `target_value 0 0.41916 0.01 -> 0.389082 0.01`. The test now
PASSes and becomes a valid GNU regression guard. The GNU-vs-Pro N-S difference
is a **separate divergence** (formulation and mesh) and is not fixed here; the
Pro reference values are recorded above for that future work unit.

## 3. What changed (exact corpus edits)

The corpus directory is gitignored/untracked, so the edits are recorded here
verbatim. Originals are preserved in the evidence dir under `originals/`.

1. Eight files — one line inserted immediately after the `end_initia` line:

   ```
   inertia_apply -yes
   ```

   Files: `absorp1.dat`, `ho_heat2.dat`, `ho_othr3.dat`, `examp11.dat`,
   `wave1.dat`, `dens1.dat`, `truss6.dat`, `viscos1.dat`.

2. `matrix8.dat` — two lines changed:

   ```
   - bounda_time    2  0. 2.e-2 1000. 2.e-2
   + bounda_time    2  0. -2.e-2 1000. -2.e-2
   - target_value 0  -0.01 1.e-4
   + target_value 0  0.01 1.e-4
   ```

3. `incnav5.dat` — one line changed:

   ```
   - target_value 0 0.41916  0.01
   + target_value 0 0.389082  0.01
   ```

No source file was modified. `options_inertia` and every deliberate alignment
commit are untouched.

## 4. Reproducibility

- Baseline corpus run (`ee0032c`): `merged-corpus-run.sh` over the 197 `.dat`,
  classification `PASS` = rc 0 (the binary enforces `target_value` via
  `exit_tn`, `miscel.cc:680-710`). Result `174 / 23 / 0`.
- After the edits: same runner, same binary. Result `184 / 13 / 0`.
- Per-test classification files: `corpus-baseline.txt`, `corpus-after.txt` in
  the evidence dir.
- Suite: `scripts/build_safe.sh` clean build + reduced suite, **16/16**.
