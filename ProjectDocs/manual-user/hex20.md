# hex20

## Description

`-hex20` is the conventional 20-node serendipity hexahedron (8 corners +
12 mid-edge nodes), the format produced by most mesh generators and FEA
codes (Abaqus C3D20, GiD, etc.) and also accepted by Tochnog
Professional. The GNU has no real hex20 element routine: it
**auto-converts** every `-hex20` volume element to the 27-node Lagrange
`-hex27` before the first calculation step, inserting the 6 face-centre
nodes and the body-centre node.

The conversion is invisible in the input: you keep writing
`element i -hex20 n1 ... n20` with the conventional node ordering and the
program rewrites the element internally. Face centres shared between
neighbouring elements are merged (one single node per shared face), so
stacked hex20 meshes stay connected.

## Improvement over plain serendipity treatment

As with `-quad8` (see [quad8](quad8.md)), accepting the conventional
serendipity input and upgrading it internally to the complete Lagrange
element is a deliberate improvement over treating `-hex20` with its plain
20-node formulation:

- The plain serendipity hex20 has an **incomplete** quadratic polynomial:
  it lacks the interior (face/body) degrees of freedom, so it cannot
  represent the full quadratic field and is less accurate under bending
  and non-constant stress gradients.
- The Lagrange `-hex27` used internally keeps the same 20 boundary nodes
  and **adds the 7 interior nodes** (6 face centres + the body centre),
  completing the quadratic interpolation on the same input mesh.
- The numerical results match Tochnog Professional, whose own suite
  states textually that hex20 elements "will be automatically converted
  to ... hex27 volume elements".

## Node ordering

The 20 nodes follow the conventional/Professional order:

1. base plane corners: BL, BR, TL, TR
2. top plane corners: BL, BR, TL, TR
3. base plane mid-edge nodes: BM, LM, RM, TM
4. top plane mid-edge nodes: BM, LM, RM, TM
5. vertical mid-edge nodes at (x0,y0), (x1,y0), (x0,y1), (x1,y1)

(BL = bottom-left, BM = mid-edge of the bottom edge, LM = mid-edge of
the left edge, etc.; the top plane repeats the base ordering at the
higher z.)

## Uso

Place the elements in the data part like any other element:

```
element 1  -hex20 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
```

## Notes

- The conversion runs automatically at every step (idempotent); no
  control record is needed. It also runs before `control_mesh_convert`,
  so a `-quad8` interface face between two hex20 solids is converted
  after the volumes are already hex27 (see [hex18](hex18.md)).
- The stress/strain interpolation of the converted element is the
  Lagrange hex27 one, exactly as the Professional does.
- The element works in 3D only (`number_of_space_dimensions 3`).

## Related

- [quad8](quad8.md) — the 2D counterpart (serendipity quad8 ->
  Lagrange quad9).
- [hex18](hex18.md) — the quadratic 3D interface element that pairs the
  quad9 faces of two hex27 solids.
- `control_mesh_convert` — conversion of interface elements.
