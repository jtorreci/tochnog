# control_mesh_convert

## Description

`control_mesh_convert` converts low-dimensional interface elements to
their isoparametric equivalent, creating the nodes of the opposite
interface side. This makes it easy to obtain a mesh with interface
elements between two blocks of material (e.g. between a pile and the
soil): generate a mesh with `-bar2` elements in the interface (with a
preprocessor like GID), add the interface group data, and use
`control_mesh_convert` to generate the interface elements.

Currently implemented (2D): `-bar2` -> `-quad4` (linear) and
`-bar3` -> `-quad6` (quadratic, 3+3 nodes). The quadratic case needs the
adjacent volume elements to be quadratic too (the Professional suite
uses `-quad8` solids, which are auto-converted to `-quad9`; the
interface bar3 sits on the shared edge).

3D: `-tria3` -> `-prism6`, `-quad4` -> `-hex8` (linear interfaces
between 3D solids) and `-quad8` -> `-hex18` (quadratic facial
interface: a `-quad8` face between two `-hex20`/`-hex27` solids, e.g.
the `interface_quad8_hex20` family, is lifted to the 18-node interface
with 2 quad9 sides; see [hex18](hex18.md)). The `-hex20` solids are
auto-converted to `-hex27` first, so the shared face centre node exists
when the quad8 face is converted.

## Uso

Place it in the data part, with the same `icontrol` index as the
`control_timestep` record:

```
element 1  -bar2 101 102
element_group 1  10
...
group_interface 10  -yes
...
control_mesh_convert 0  -yes
control_mesh_convert_element_group 0  0 1
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `0`       | Index of the control record. Must match the `control_timestep` index. |
| `-yes`    | Convert interface elements (`-bar2` -> `-quad4`, `-bar3` -> `-quad6`, `-tria3` -> `-prism6`, `-quad4` -> `-hex8`, 3D `-quad8` -> `-hex18`). |

## Related

- `control_mesh_convert_element_group index g_0 g_1 ...` — element groups
  located on ONE side of the interfaces; the neighbours in these groups
  keep the original nodes, the others are reconnected to the new nodes.
- `group_interface` — marks an element group as an interface.

## Example

A bar2 interface shared by two blocks is converted to a quad4, creating
the second side of the interface and reconnecting the block on the other
side:

```
element 1  -quad4 1 2 3 4
element 2  -bar2 2 4
element 3  -quad4 2 4 5 6
element_group 1  0
element_group 2  10
element_group 3  0
group_interface 10  -yes
group_interface_materi_elasti_stiffness 10  1000.0 0.0 0.0
control_mesh_convert 0  -yes
control_timestep    0  0.1 1.
```
