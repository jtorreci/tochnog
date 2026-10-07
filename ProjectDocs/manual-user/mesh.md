# mesh (short alias of options_mesh)

## Description

`mesh` (manual Professional 6.799+) is the short alias of
`options_mesh`: it sets the mesh motion for the calculation:

```
mesh -fixed_in_space -fixed_in_space
```

| Value | Meaning |
|-------|---------|
| `-fixed_in_space` | the mesh does NOT follow the material (Eulerian-ish) |
| `-follow_material` | the mesh follows the material (Lagrangian) |

The GNU had a dead `MESH` placeholder record (no type/class); routing
the keyword through the exact-name loop returned `MESH` with
`no_index=0` and the parser tried to read `-fixed_in_space` as an index.
The keyword is now routed to `OPTIONS_MESH` (which has `no_index=1` and
an INTEGER data of length ndim).

## Default

Without a `mesh` / `options_mesh` record the mesh is `-fixed_in_space`
(manual Professional 6.806): the nodes do not move, also with
`materi_velocity` only. Give `mesh -follow_material ...` explicitly for an
updated-Lagrange (moving) mesh. With `materi_displacement` the mesh is
always `-fixed_in_space`, whatever the record says. Until F3
(2026-10-07) the GNU default was `-follow_material`, so inputs written
for the old GNU that relied on the moving mesh now need the explicit
record.

Example (velocity-only bar stretched by 50 percent, no `mesh` record):
`sigxx` 0.5, as the Professional (the old GNU default gave 0.407136).
Self-checking inputs: `validation-suite/test-2026/mesh_default_fixed_in_space1.dat`
(stretch) and `mesh_default_fixed_in_space2.dat` (rigid rotation: zero
stress).

## Usage

```
mesh -fixed_in_space -fixed_in_space
```

## Notes

- slope_classical_numerical uses this keyword and now parses/runs.

## Differences with the Professional version

- **default `mesh` kinematics** (C8): RESOLVED in F3 (2026-10-07). The
  default is now `-fixed_in_space`, as the Professional.
- **explicit `-follow_material`: in-step geometry** (known difference,
  documented). GNU evaluates the velocity gradient of a `-follow_material`
  step on the END-of-step geometry (nodes already advanced by the
  iteration), the Professional on the START-of-step geometry. The results
  converge with the number of steps: 2-quad4 column sheared by 1 with
  `-updated_linear`, node 3: one step GNU 0.47222 / Pro 0.5; 100 steps
  GNU 0.52048 / Pro 0.52093. Status: Open, LOW (F3, 2026-10-07). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
