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
