# control_bounda_relax

## Files and functions

- `database.cc` — registrations: `control_bounda_relax` (INTEGER, 1 value, CONTROL) and `control_bounda_relax_geometry` (INTEGER, 2 values, CONTROL, `data_required = CONTROL_BOUNDA_RELAX`).
- `tochnog.h` / `tochnog-mod.h` — enums `CONTROL_BOUNDA_RELAX`, `CONTROL_BOUNDA_RELAX_GEOMETRY` appended before `LAST_DUMMY` (both headers, same order).

## Implementation details

Parse-only registration: both records are accepted and stored in the
database; no code consumes them yet. The storing machinery of the
Professional (per-node record of the RHS/reaction of the prescribed
dofs on the relax geometry, consumed later by the relaxing
`bounda_force`) is PENDING.

Why parse-only is enough for the current corpus tests (measured on the
Professional binary 25-10-2023, A/B with and without the relax records
on `relax1.dat`/`relax2.dat`): the relaxation observable in both tests
is driven entirely by the `bounda_force` + `bounda_time` time function
(`bounda.cc` sets `node_rhside = factor*load`, and the load window
`[time0, time1]` of the record makes the prescription expire at the end
of its window). With and without `control_bounda_relax*` the Pro
produces identical results, and so does this GNU version.

## Status

- `relax1` rc=0 (target `node_dof 2 -disx = 0.5`; the model uses
  `materi_displacement` + `group_materi_memory -total_linear`).
- `relax2` RUNFAIL with a documented, ORTHOGONAL root cause (NOT the
  relax machinery): the model is velocity-only (no `materi_displacement`,
  no `group_materi_memory`). This GNU version defaults `options_mesh` to
  `-follow_material` for velocity-only analyses (`top.cc`), and
  `locate()` (`locate.cc`) moves the reference mesh by `v*dt` every
  iteration: the constitutive integration then runs over the deformed
  geometry and produces large-strain kinematics. Measured: the reaction
  at node 1 follows `sigma = E*ln(1+u)` (0.7188 at time 1, 0.7677 at
  time 3) instead of the linear `E*u` of the Professional (1.0/1.0,
  whose velocity-only reference mesh stays fixed; `options_mesh` is a
  legacy GNU record the Professional 25-10-2023 rejects). A/B: adding
  `options_mesh -fixed_in_space` makes the GNU reproduce the
  Professional EXACTLY (reaction 1.0 at t=1, 0.0 at t=2, 1.0 at t=3,
  rc=0). Aligning the velocity-only default with the Professional is a
  core-kinematics work unit with corpus-wide blast radius (PENDING).

## Hardcoded parameters / pending refactorings

- Implement the storing semantics per the Professional manual 6.108
  (store the RHS of prescribed dofs on the relax geometry into a
  per-node record) and the consumption (relaxing `bounda_force`
  scaled by the stored value) — no corpus test currently exercises the
  stored value, so a dedicated verification example is needed.
- Resolve the velocity-only follow-material default divergence
  (documented above).
