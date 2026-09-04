# quad8

## Description

`-quad8` is the Professional 8-node serendipity plane element. The GNU
does not have a real quad8: like the Professional, it **auto-converts**
every `-quad8` volume element to the 9-node Lagrange `-quad9` before the
first calculation step, inserting the centre node. This makes the corpus
`.dat` files of the Professional suite (which state textually that "the
quad8 ... will be automatically converted to quad9 surface elements")
run without manual mesh editing.

The conversion is invisible in the input: you keep writing
`element i -quad8 n1 ... n8` with the Professional node ordering
(corners BL, BR, TL, TR followed by the mid-edge nodes BM, LM, RM, TM)
and the program rewrites the element internally.

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
  conversion.
- `-hex20` is not implemented (pending).

## Related

- `control_mesh_convert` — conversion of quadratic interfaces
  (`-bar3` -> `-quad6`).
- `control_mesh_extrude` — quadratic extrusion quad9 -> hex27.
