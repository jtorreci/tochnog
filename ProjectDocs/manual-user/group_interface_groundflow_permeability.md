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

## Per unit length (fixed 2026-10-06, audit K04 `ground9`)

`pe` is the permeability per unit length (2D) or per unit area (3D),
manual Professional 6.628: the total flow through the interface is
`pe * L * (p_side1 - p_side2)`. Before 2026-10-06 the GNU gave each facing
node pair the full `pe`, i.e. `2*pe` for a quad4 interface whatever its
length (Pro test `ground9`: middle pressure 1/3 instead of 0.5). Now each
pair carries its tributary length (Lobatto weight x side length), as the
mechanical interface springs; results equal the Professional 02-08-2026.
Tests: `validation-suite/test-2026/interface_flux_measure_ground1.dat`.
