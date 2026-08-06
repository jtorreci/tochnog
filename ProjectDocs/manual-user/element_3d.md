# element_3d

## Description

Solid 3D elements available in this GNU version:

| Element | Nodes | Geometry |
|---------|-------|----------|
| `-tet4` | 4 | Linear tetrahedron (4 corner nodes). |
| `-tet10` | 10 | Quadratic tetrahedron (4 corners + 6 edge midpoints). |
| `-hex8` | 8 | Linear hexahedron / brick (8 corner nodes). |
| `-prism6` | 6 | Linear prism / wedge: triangular base (nodes 1-3) and triangular top (nodes 4-6). |

The model must run in 3D, so set `number_of_space_dimensions 3`, and give
enough integration points for the element:

| Element | Required `number_of_integration_points` |
|---------|------------------------------------------|
| `-hex8` | 8 (Gauss 2×2×2) |
| `-prism6` | 6 (2 z-levels × 3 area points) |
| `-tet4` | 1 (minimal) or 4 (maximal) |
| `-tet10` | 10 |

## Usage

```
node <inod> x y z
element <ielem> -tet4|-tet10|-hex8|-prism6 <node list>
```

Nodes carry three coordinates; the element record lists its node numbers in the
order shown above (prism: 3 base nodes then 3 top nodes; hex8: 4 bottom nodes
then 4 top nodes; tet10: 4 corners then 6 edge midpoints).

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `inod` | Node number. |
| `x y z` | Node coordinates in the three space dimensions. |
| `ielem` | Element number. |
| element type | `-tet4`, `-tet10`, `-hex8` or `-prism6`. |
| node list | Connectivity: node numbers of the element. |

## Example

A single `-hex8` brick:

```
number_of_space_dimensions 3
end_initia

node 1  0. 0. 0.
node 2  1. 0. 0.
node 3  1. 1. 0.
node 4  0. 1. 0.
node 5  0. 0. 1.
node 6  1. 0. 1.
node 7  1. 1. 1.
node 8  0. 1. 1.
element 1  -hex8 1 2 3 4 5 6 7 8
end_data
```

A single `-prism6` (3 base nodes + 3 top nodes):

```
number_of_space_dimensions 3
end_initia

node 1  0. 0. 0.
node 2  1. 0. 0.
node 3  0. 1. 0.
node 4  0. 0. 1.
node 5  1. 0. 1.
node 6  0. 1. 1.
element 1  -prism6 1 2 3 4 5 6
end_data
```
