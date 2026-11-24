# control_materi_updated_apply

## Description

`control_materi_updated_apply index switch` (manual Professional 6.152)
changes the material memory for the timesteps of the `control_timestep`
record with the same index:

- `-no`: any `-updated` (and `-updated_without_rotation`) material
  memory is set to `-updated_linear`.
- `-yes`: any NON-specified material memory (group without a
  `group_materi_memory` record) is set to `-updated`. Explicit memory
  records are not changed.

Without the record nothing changes: an unspecified memory is
`-updated_linear` (see [group_materi_memory_default](group_materi_memory_default.md)).

It is useful to run a first phase with small-deformation kinematics and a
later phase with updated (large rotation) kinematics, or the other way
round, in one input file.

## Usage

```
control_timestep 0  0.01 1.
control_materi_updated_apply 0  -yes
```

Parameters: `index` = the control index (same as the `control_timestep`
it applies to); `switch` = `-yes` or `-no`.

## Example

Single quad4, simple shear gamma = 1, young 1, poisson 0, no memory
record, plus `control_materi_updated_apply 0 -yes`: sigxx 0.231098,
sigxy 0.420739, sigyy -0.228598 (the `-updated` values). Input:
`validation-suite/test-2026/materi_memory_default_updated_apply1.dat`.

## Differences with the Professional version

- The Professional 02-08-2026 binary ignores `-yes` for an unspecified
  memory (it gives sigxx 0, sigxy 0.5, sigyy 0 in the example above, for
  displacement, `-follow_material` and `-fixed_in_space` kinematics).
  The GNU follows the manual. Deliberate divergence (decision
  2026-10-04).
- `-no` behaves the same in both versions.
