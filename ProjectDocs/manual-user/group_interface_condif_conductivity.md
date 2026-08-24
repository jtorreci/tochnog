# group_interface_condif_conductivity

## Description

Heat flow through the interface in its thickness direction per unit
temperature difference between the two sides (manual Professional
6.625): `q = k * (T_side1 - T_side2)`, assembled on the temperature
dofs of the facing node pairs with a symmetric tangent (the thermal
analog of `group_interface_groundflow_permeability`).

The conductivity k is NOT the material conductivity: it is the
conductivity of the layer simulated by the interface, incorporating the
thermal thickness of the interface. Units: [power]/[temperature*length]
in 2D, [power]/[temperature*length*length] in 3D.

Requires `condif_temperature`.

## Usage

```
group_interface_condif_conductivity <index> <k>
```

## Example

```
group_type 10  -condif
group_interface 10  -yes
group_interface_condif_conductivity 10  1.0
```

Left block with prescribed T=2, right block thermally isolated
(conductivity 1e-6): the interface heat flux alone drives the right
block to T=2 in steady state (test `iface_condif`; without the record
it stays at 0).
