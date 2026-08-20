# group_interface_groundflow_permeability

## Description

Specifies the permeability per unit length (2D) or unit area (3D) for
interface elements in a ground flow analysis. The interface couples the pore
pressures on its two sides with a through-flux:

```
q = pe * (pres_side1 - pres_side2)
```

The flux is assembled on the pressure dofs of the facing node pairs of the
interface, so water can flow across an interface that separates two flow
domains (e.g. a low-permeability barrier between two aquifers).

## Usage

```
group_interface_groundflow_permeability <element_group> <pe>
```

## Parameters

| Parameter | Meaning                                                       |
|-----------|---------------------------------------------------------------|
| `element_group` | Element group, see `element_group`.                   |
| `pe`       | Interface permeability per unit length (2D) or area (3D).    |

## Example

```
group_type 10  -groundflow
group_interface 10  -yes
group_interface_groundflow_permeability 10  1.0
```
