# mpc_node_number / mpc_node_factor

## Description

Multi point constraints between nodal degrees of freedom (manual
Professional 6.874/6.875). The dof `dof_0` of node `node_0` (the
**slave**) is constrained to a linear combination of dofs of the
**master** nodes:

```
dof_0(node_0) = factor_10 * dof_10(node_1) + factor_11 * dof_11(node_1)
                + factor_20 * dof_20(node_2) + ...
```

The factors come from `mpc_node_factor` in the same order as the master
dofs appear in `mpc_node_number`; factors that are not specified default
to 1.

The slave dof is treated as a KNOWN quantity (like a boundary condition):
its value is recomputed every equilibrium iteration from the current
master values. The slave's equilibrium equation is not solved and no
force redistribution to the masters takes place. Boundary conditions with
`bounda_dof`/`bounda_time` must NOT be specified on slave nodes
(manual 6.875).

## Syntax

```
mpc_node_number <index> <node_0> <dof_0> <node_1> <dof_10> [<dof_11> ...] <node_2> <dof_20> ...
mpc_node_factor <index> <factor_10> [<factor_11> ... <factor_20> ...]
```

- `node_0` is the slave node, `dof_0` the constrained dof (a principal
  dof keyword such as `-velx`, `-vely`, `-temp`, `-pres`).
- The master nodes follow: after each node number, ALL its dofs are
  listed, then the next node number starts a new master group.
- `mpc_node_factor` lists the multiplication factors in exactly the same
  (flattened) order.

## Example

```
node 1  0
node 2  1
node 3  2
element 1  -bar2 1 2
element 2  -bar2 2 3

bounda_dof   10  1 -velx
bounda_time  10  1.

mpc_node_number 10 2 -velx 1 -velx
mpc_node_factor 10 2.
```

Node 2 (slave) gets `velx_2 = 2 * velx_1`; with the boundary condition
`velx_1 = 1` the slave ends at `velx_2 = 2` (manual example, corpus
test `mpc1`).

## Tests

Corpus tests that pass with this record: `mpc1` (node_dof identical to
the Professional .dbs: slave node 2 velx = disx = 2, sigxx = 0.5).
`mohr_coul_direct2/4/7` apply the tie correctly (the tied nodes are
exactly equal, verified) but stay RUNFAIL on the direct plastic material
family (mohr_coul_direct1/3/5/6/8 — no mpc — fail the same way).

## Notes

- Only principal dofs can be constrained (manual 6.875).
- The slave value is recomputed every iteration from the CURRENT master
  values, so the constraint is exact at convergence.
- `mpc_apply` (manual 6.859, switch `-yes`/`-no`, default `-yes`) is
  registered as a known keyword but NOT yet consumed: all mpc records
  are always active.
- `mpc_geometry`/`mpc_element_group` (automatic generation of
  `mpc_node_number`/`mpc_node_factor` records) remain PENDING; see
  `mpc_linear_quadratic.md` for the implemented auto-generation case.
