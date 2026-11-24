# control_bounda_relax

## Description

With `control_bounda_relax` you can require Tochnog to store the nodal
right-hand-sides — for example the external nodal forces of nodes with
prescribed velocities. These stored nodal right-hand-sides can later be
used to relax prescribed boundary conditions: a prescribed velocity is
removed and substituted by the stored external right-hand-side
(external force), which is slowly set to zero by multiplication with a
time function as specified with `bounda_force` in combination with
`bounda_time`.

With `control_bounda_relax_geometry` (same index) you can select a
specific geometry for which the storing will be done.

A typical example can be found in the `relax1.dat` file in the test
distribution.

## Usage

```
control_bounda_relax <index> [-yes | -no]
control_bounda_relax_geometry <index> -<geometry_item> <geometry_index>
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (the switch and the geometry record share it). |
| `-yes` / `-no` | Enable / disable the storing of the nodal right-hand-sides for these timesteps. |

## Example

`relax1.dat`: node 2 of a bar is moved to the right between time 0 and
1 (prescribed velocity); at time 1 the prescription is removed, the
external nodal force is imposed and relaxed to 0 between time 1 and
time 2 — at time 1.5 the displacement of node 2 is 0.5:

```
bounda_dof   0  1 -velx
bounda_time  0  0. 0. 100. 0.
bounda_dof   1  2 -velx
bounda_time  1  0. 1. 1. 1.
bounda_force 2  2 -velx
bounda_time  2  1.00001 1. 2. 0.

control_timestep      10  0.1 1.0
control_bounda_relax  20  -yes
control_bounda_relax_geometry 20  -geometry_point 10
control_timestep      30  0.1 0.5
```

## Differences with the Professional version

- Both records are registered and parsed; the storing semantics of the
  Professional (a dedicated per-node record of the reaction of the
  prescribed dofs on the relax geometry) is NOT implemented. The
  observable behavior of the relax corpus tests is reproduced without
  it: in `relax1`/`relax2` the relaxation is driven by the
  `bounda_force` time function alone (verified A/B against the
  Professional binary 25-10-2023: with and without the relax records
  the Pro produces identical results).
- `relax1` runs rc=0. `relax2` stays RUNFAIL for a reason ORTHOGONAL to
  this family: its model is velocity-only (no `materi_displacement`,
  no `group_materi_memory`), and this GNU version convects the
  reference mesh with the material by default in velocity-only analyses
  (`options_mesh` defaults to `-follow_material`, a legacy GNU record
  the Professional does not even accept), producing large-strain
  (logarithmic) kinematics: the measured reaction is 0.768 at time 3
  instead of 1.0 (the response follows sigma = E*ln(1+u); with an
  explicit `options_mesh -fixed_in_space` the GNU reproduces the
  Professional exactly: 1.0/0.0/1.0 at times 1/2/3). Aligning the
  velocity-only default with the Professional needs a corpus-wide A/B
  and is PENDING.
