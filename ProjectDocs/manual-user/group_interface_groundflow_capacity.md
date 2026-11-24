# group_interface_groundflow_capacity

## Description

Specifies the capacity for interface elements in a ground flow analysis. The
interface stores water when its pore pressure changes (storage term on the
pressure dofs of the interface nodes).

## Usage

```
group_interface_groundflow_capacity <element_group> <C>
```

## Parameters

| Parameter | Meaning                                                |
|-----------|--------------------------------------------------------|
| `element_group` | Element group, see `element_group`.            |
| `C`        | Storage capacity of the interface.                     |

## Example

```
group_type 10  -groundflow
group_interface 10  -yes
group_interface_groundflow_capacity 10  1.0
```

## Per unit length, implicit (fixed 2026-10-06, audit K04)

The capacity is per unit length (2D) / area (3D) like the permeability:
each facing node pair stores `C * w_i * L`, split evenly over its two
nodes, and the storage term is implicit (as the continuum
`group_groundflow_capacity`). This reproduces the Professional 02-08-2026
transient of `groundflow_interface` (node 5: 1.0, 1.75, 1.998 at t = 1,
3, 10) and makes it independent of the interface length. Before, the term
was explicit and unweighted (0.32 / 0.81 / 1.65). Test:
`validation-suite/test-2026/interface_capacity_transient1.dat`.
