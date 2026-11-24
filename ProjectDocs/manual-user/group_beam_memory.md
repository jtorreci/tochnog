# group_beam_memory

## Description

`group_beam_memory <group> -option` chooses how a beam follows the
deformation (manual Professional 6.591):

- `-total_linear` (default): geometrically linear beam. The local beam
  axes stay those of the initial geometry. This gives the classical
  strength-of-materials results: a cantilever of length 1 with tip
  displacement 0.01 gets the force `12 E I v / L^3` exactly.
- `-updated`: the local beam axes follow the deformed beam (the chord
  from the first to the second node in the current geometry), and the
  forces of the previous step are rotated with the beam. Useful for
  large rotations of a beam; for small deformations the results differ
  from `-total_linear` by the chord rotation (the cantilever above gets
  0.2399880 instead of 0.24, the chord rotates by 0.01).
- `-updated_without_rotation`: GNU name of the same behaviour as
  `-total_linear`.

The default was `-updated` in the GNU version until 2026-11-02 (F5 T2
V4d, owner decision, following the Professional). Inputs that relied on
the old default must now give `group_beam_memory <group> -updated`.

## Example

`group_beam_memory_default1.dat` (validation-suite/test-2026): 3D
cantilever, E 2, Izz 4, tip `vely` 0.01 prescribed, no
`group_beam_memory`:

```
element_beam_force_moment 1 0 -0.24 0 0 0 -0.24 0 0.24 0 0 0 0
```

`group_beam_memory_updated1.dat`: the same with
`group_beam_memory 0 -updated`: -0.2399880.

## Differences with the Professional version

- The Professional 02-08-2026 stops with `materi_displacement should not
  be initialised.` for `-updated` / `-updated_linear` on these inputs;
  the GNU computes the updated beam (BEAM-MEMORY-UPDATED,
  [DIVERGENCES-PROFESSIONAL](../DIVERGENCES-PROFESSIONAL.md)).
