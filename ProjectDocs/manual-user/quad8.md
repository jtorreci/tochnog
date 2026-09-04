# quad8

## Description

`-quad8` is the conventional 8-node serendipity plane element (corners
plus mid-edge nodes), the format produced by most mesh generators and
FEA codes (Abaqus CP8, GiD, etc.) and also accepted by Tochnog
Professional. The GNU has no real quad8 element routine: it
**auto-converts** every `-quad8` volume element to the 9-node Lagrange
`-quad9` before the first calculation step, inserting the centre node.

The conversion is invisible in the input: you keep writing
`element i -quad8 n1 ... n8` with the conventional node ordering
(corners BL, BR, TL, TR followed by the mid-edge nodes BM, LM, RM, TM)
and the program rewrites the element internally.

## Improvement over plain serendipity treatment

Accepting the conventional serendipity input and upgrading it internally
to the complete Lagrange element is a deliberate improvement of this
version over treating `-quad8` with its plain 8-node formulation:

- The plain serendipity quad8 has an **incomplete** quadratic polynomial:
  it lacks the interior (centre) degree of freedom, so it cannot represent
  the full quadratic displacement field (in particular the `xi^2*eta^2`
  bubble term) and is less accurate under bending and non-constant stress
  gradients.
- The Lagrange `-quad9` used internally keeps the same 8 boundary nodes
  and **adds the centre node**, completing the quadratic interpolation.
  The same input mesh therefore receives a strictly richer element at no
  cost to the user.
- The same upgrade applies to the 3D counterpart `-hex20` (serendipity
  20-node hexahedron): it is elevated to the complete 27-node Lagrange
  `-hex27` (see [hex20](hex20.md)).
- The numerical results are verified identical to Tochnog Professional,
  which follows the same auto-conversion approach (its own suite states
  textually that "the quad8 ... will be automatically converted to quad9
  surface elements"). GNU node_dof agrees with the Professional `.dbs` on
  the corpus patch tests.

## Uso

Place the elements in the data part like any other element:

```
element 1  -quad8 1 2 3 4 5 6 7 8
```

The 8 nodes follow the Professional order: corners bottom-left,
bottom-right, top-left, top-right; then the mid-edge nodes bottom,
left, right, top.

## Notes

- The conversion runs automatically at the first step; no control record
  is needed (also works together with `control_mesh_convert` for
  quadratic interfaces and with `control_mesh_extrude`, which then
  extrudes the quad9 to hex27).
- The stress/strain interpolation of the converted element is the
  Lagrange quad9 one, exactly as the Professional does.
- `-quad8` used as an INTERFACE element (group with
  `group_interface -yes`, the 3D `interface_quad8_hex20` family) is NOT
  converted to quad9: those are facial elements handled by the interface
  conversion (lifted to the quadratic 3D interface, see [hex18](hex18.md)).
- `-hex20` (3D serendipity, elevated to `-hex27`) is implemented by the
  same Lagrange upgrade: see [hex20](hex20.md).

## Related

- `control_mesh_convert` — conversion of quadratic interfaces
  (`-bar3` -> `-quad6`, and the 3D `-quad8` face -> `-hex18`).
- `control_mesh_extrude` — quadratic extrusion quad9 -> hex27.
- `hex20` — the 3D serendipity element elevated to hex27.
- `hex18` — the quadratic 3D interface element.
