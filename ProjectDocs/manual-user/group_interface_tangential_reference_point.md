# group_interface_tangential_reference_point

## Description

Reference point that defines the first tangential direction of a 3D
interface element (manual Professional 6.636). The frame is built as:

1. The normal direction to the interface plane (unchanged).
2. A vector from the element centroid to the reference point.
3. The part of that vector PERPENDICULAR to the normal is the first
   tangential direction t1.
4. `t2 = normal x t1` (second tangential direction).

Without this record the GNU default applies (t1 along the first side-1
edge). If the reference point lies on the normal line through the
centroid (degenerate perpendicular part), the default frame is kept.
2D interfaces are unaffected (single tangent).

Use it to control in which basis the tangential forces
(`element_interface_force_tang` / `_tang2`) and the Mohr-Coulomb
friction components are expressed.

## Usage

```
group_interface_tangential_reference_point <index> <x> <y> <z>
```

## Example

```
group_interface 10  -yes
group_interface_tangential_reference_point 10  1. 0.5 5.
```

3D interface with normal +x (default t1=+y, t2=+z). A reference point
in +z rotates t1 to +z: a shear load in y then lands entirely in
`element_interface_force_tang2` while `_tang` stays ~0 (test
`iface_tangref`; with the default frame the shear lands in `_tang`).
