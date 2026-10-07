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

## Usage

```
mesh -fixed_in_space -fixed_in_space
```

## Notes

- slope_classical_numerical uses this keyword and now parses/runs.

## Differences with the Professional version

- **default `mesh` kinematics** (known difference, not decided; BACKLOG C8).
  GNU: `-follow_material` (nodes move with `materi_velocity` and no
  `materi_displacement`). Professional 02-08-2026: `-fixed_in_space` (manual
  6.806); 1D bar +50% `-updated`: Pro 0.5, GNU 0.40714. Status: Open,
  MEDIUM, effort large (corpus relies on the moving mesh); tied to C18/C19.
  Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
