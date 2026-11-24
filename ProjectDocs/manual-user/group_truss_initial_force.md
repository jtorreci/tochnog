# group_truss_initial_force

## Description

With `group_truss_initial_force` you can give the trusses of a truss
element group an initial internal (axial) force from the moment the
element comes to life (manual Professional 6.776). The incremental
force update of the truss (`F = F_old + E*A/L * dL`) accumulates on top
of this initial value, so a truss with a pre-tension keeps it while it
deforms elastically.

The record does NOT create an initial elastic strain: the output
`element_truss_strain` stays 0 for a stationary truss with an initial
force (measured on the Professional `.dbs` of `truss12.dat`).

## Usage

```
group_type <group_index> -truss
group_truss_initial_force <group_index> <initial_force>
```

## Parameters

| Parameter       | Meaning |
|-----------------|---------|
| `initial_force` | Axial force of the truss at the start of its life (tension positive). |

## Example

See `truss12.dat`: a truss of E=A=1 between the nodes (0,0)-(0,1) with
`group_truss_initial_force = 1` inside a fully fixed tria3 produces the
nodal right-hand-sides `node_rhside: node 1 vely = +1`, `node 3 vely =
-1` (verified against the Professional `.dbs`).
