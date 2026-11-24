# hex18

## Description

`-hex18` is the quadratic 3D interface element: two `-quad9` sides of 9
nodes each (ns1 = 9 facing node pairs), the 3D counterpart of the 2D
`-quad6` interface. It connects the quad9 faces of two quadratic solids
(`-hex27`, or `-hex20` auto-converted to hex27) with a joint whose
strains are the displacement differences between the facing pairs.

The element can appear in the input directly (interface3 of the suite
writes `element 1 -hex18 ...` between two quad9 faces at z=0 and z=1),
but its typical origin is the **auto-conversion** of a `-quad8` face:
when `control_mesh_convert` meets a `-quad8` facial interface element
(`group_interface -yes`, the `interface_quad8_hex20` family), it lifts it
to `-hex18`:

- the 9th side-1 node is the face centre (the average of the 4 corners) -
  the node that the hex20 -> hex27 conversion created on the shared face,
  looked up by coordinates (created if missing);
- the side-1 ordering is rewritten to the GNU `-quad9` tensor order
  `[BL, BM, BR, LM, C, RM, TL, TM, TR]`;
- the 9 nodes of side 2 are new copies of side 1 shifted along the
  interface normal, and the solid on the other side of the interface is
  reconnected to them.

## Node ordering

The 18 nodes are the 9 nodes of side 1 (quad9 tensor order) followed by
the 9 nodes of side 2 (same order):

```
element i  -hex18 s1 s2 s3 s4 s5 s6 s7 s8 s9  t1 t2 t3 t4 t5 t6 t7 t8 t9
```

side 1: BL, BM, BR, LM, CENTRE, RM, TL, TM, TR
side 2: the same positions on the opposite side.

## Integration weights

The interface is integrated per facing pair with the 2D product of the
1D Gauss-Lobatto rule (1, 4, 1) per direction over the quad9 face:
corner pairs weight 1/36, mid-edge pairs 4/36 and the centre pair 16/36
of the element measure (face area). This matches the load distribution
the Professional uses on quadratic faces (interface3 of the suite loads
the top side with -1 on the corners, -4 on the mid-edge nodes and -16 on
the centre, and every integration point reports the same uniform
stress -9).

## Uso

Interface group data is the same as for the linear interfaces
(`group_interface -yes` + `group_interface_materi_elasti_stiffness`):

```
element 1  -hex18 1 2 3 4 5 6 7 8 9  11 12 13 14 15 16 17 18 19
group_interface 0  -yes
group_interface_materi_memory 0  -total_linear
group_interface_materi_elasti_stiffness 0  0.1e11  0.5e10  0.5e10
```

## Notes

- A `-quad8` face is only converted when its element group has
  `group_interface -yes`; a `-quad8` of a material group is a volume
  element and is converted to `-quad9` by the quad8 conversion instead.
- The conversion needs the neighbouring volumes to be quadratic already:
  with `-hex20` solids the automatic hex20 -> hex27 conversion runs first
  (every step, before `control_mesh_convert`), so the shared face centre
  exists when the quad8 face is lifted to hex18.

## Related

- [hex20](hex20.md) — the 3D serendipity solid auto-converted to hex27.
- [quad8](quad8.md) — the 2D quadratic solid conversion.
- `control_mesh_convert` — the record that triggers the quad8 -> hex18
  conversion.
- `group_interface` — marks an element group as an interface.
