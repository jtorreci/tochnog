# group_materi_memory -updated_jaumann

## Description

`group_materi_memory -updated_jaumann` (manual Professional 2.2.1 and
6.686) is an updated Lagrange memory with a first-order (Jaumann-like)
rotation of the stresses:

- the incremental strain is the LINEAR engineering strain of the
  incremental deformation tensor F on the current configuration (as
  `-updated_linear`);
- the stresses of the previous time point plus the extra stresses are
  rotated with R = 0.5 (F - F^T) + I (the first-order rotation, not the
  polar rotation that `-updated` uses).

It is useful for large rotations in an updated formulation where the
first-order Jaumann update of the Professional must be reproduced.

| Memory | strain increment | stress rotation |
|--------|------------------|-----------------|
| `-updated_linear` | linear | none |
| `-updated_jaumann` | linear | R = I + 0.5 (F - F^T) |
| `-updated` | from the polar decomposition | polar R |

## Usage

```
group_materi_memory  index  -updated_jaumann
```

`index` is the element group. The value is accepted with
`materi_velocity` (with or without `materi_velocity_integrated`) and also
together with `materi_displacement`, as in the Professional.

## Example

Single quad4 in simple shear gamma = 0.5, young 1, poisson 0:

```
echo -no
number_of_space_dimensions 2
materi_velocity
materi_velocity_integrated
materi_stress
end_initia
node 1 0 0
node 2 1 0
node 3 0 1
node 4 1 1
element 1 -quad4 1 2 3 4
group_type 0 -materi
group_materi_elasti_young 0 1.0
group_materi_elasti_poisson 0 0.0
group_materi_memory 0 -updated_jaumann
mesh -fixed_in_space -fixed_in_space
bounda_dof 0 -ra 1 2 -ra -velx
bounda_dof 1 -ra 3 4 -ra -velx
bounda_dof 2 -ra 1 2 3 4 -ra -vely
bounda_time 1 0. 0.5 1. 0.5
control_timestep 0 1.0 1.
control_timestep_iterations 0 20
end_data
```

Result (one step): sigxx 0.125, sigxy 0.234375, sigyy -0.125; in 100
steps 0.0618335 / 0.2396333 / -0.0618335. Both equal to the Professional
02-08-2026. Self-checking inputs:
`validation-suite/test-2026/materi_memory_updated_jaumann_shear1.dat` and
`materi_memory_updated_jaumann_shear2.dat`.

## Notes

- Implemented in F3 (2026-10-07, BACKLOG C20). Before, the GNU rejected
  the value ("Error in data part").
- Tendons treat it like `-updated_linear` (tendon direction not rotated).

## Differences with the Professional version

- **explicit `-follow_material`: in-step geometry** (known difference).
  On a `mesh -follow_material` the GNU builds the increment on the
  END-of-step geometry, the Professional on the START-of-step geometry,
  so one large step differs (2-quad4 column sheared by 1 in one step:
  GNU sigxx 0.0207, Pro -0.0216; with `-fixed_in_space` both -0.0216).
  The difference vanishes with small steps. Status: Open, LOW. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
