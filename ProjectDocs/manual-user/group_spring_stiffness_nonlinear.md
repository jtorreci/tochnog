# group_spring_stiffness_nonlinear

Nonlinear spring stiffness diagram (manual Professional 6.768): pairs
`epsilon0 k0 epsilon1 k1 ...` of total spring strain (= total spring elongation)
vs spring stiffness. The index is the element group.

When the record is present (and the linear `group_spring_stiffness` is absent),
the stiffness of each increment is read from the diagram at the MIDPOINT total
strain of the increment. For a piecewise linear diagram and a linear strain path
this integrates the spring force exactly (spring6 of the corpus checks
`element_spring_force = 0.5` EXACT for `k(eps) = eps` over an elongation of 1).

Example (spring6.dat): stiffness grows 0 -> 10 between strain 0 and 10:
```
group_type 0  -spring
group_spring_memory 0  -total_linear
group_spring_stiffness_nonlinear 0  0. 0. 10. 10.0
```
See also element_spring_strain.
