# group_spring_memory

Spring memory model selector (manual Professional 6.765). Value: `-updated_linear`,
`-total_linear` or `-updated`; the index is the element group (see element_group).

The GNU spring law is incremental on the current configuration, which matches
`-updated_linear`. In the 1D corpus models (spring1/spring6) `-total_linear`
coincides with it, so the record is accepted and stored; a geometrically nonlinear
`-updated` large-rotation spring is not implemented (the two-noded spring of
spring.cc already rotates with its nodes).

Example:
```
group_type 0  -spring
group_spring_memory 0  -total_linear
group_spring_stiffness 0  1.0
```
Unlocks spring1 (rc=0).
