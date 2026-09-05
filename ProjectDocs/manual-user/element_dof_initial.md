# element_dof_initial

## Description

With `element_dof_initial` you can give an element a "past" dof field:
when the element comes the first time to live it assumes it had in the
past the dofs of this record (manual Professional 6.422). The record
influences the inertia terms of the transient integration — for a
temperature analysis the initial temperature of the element:

- birth step of a conduction element with an initial temperature
  `T0`: the capacity term integrates `C*(T - T0)/dt` instead of
  `C*(T - T_old)/dt`, so a born element "remembers" the initial field
  (phased-analysis / initial-condition modelling).

You specify one value per element dof (the same value applies to all
nodes of the element); with a single value and several dofs the value
is reused for all of them.

## Usage

```
element_dof_initial <element_number> <dof_0> <dof_1> ...
```

## Example

See `phase1.dat`: a 1D conduction bar (`condif_temperature`, rho=1,
capacity=2, conductivity=1, lumped inertia 1) with the left node
prescribed to T=1 and `element_dof_initial 1 0.4` reaches the dynamic
equilibrium `node_dof 2 -temp = 0.7` at t=1: the conduction term
`1*(1-0.7)` balances the inertia `1*(0.7-0.4)` (target, verified
against the Professional `.dbs`).
