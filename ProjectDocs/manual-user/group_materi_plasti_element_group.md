# group_materi_plasti_element_group

Frictional slip of granular materials on other materials (concrete, steel,
...; manual Professional 6.698). For a granular element group you list the
groups of the "hard" materials (group_0 group_1 ...). The friction angle phi,
the dilatancy angle phiflow and the cohesion c of the granular material are
reduced with a factor (default 2./3.) for granular elements which are a DIRECT
NEIGHBOR of an element of one of the listed groups.

The index is the element group of the granular material (see element_group).

REGISTERED in this batch; the neighbor-based strength reduction is PENDING
(consumption belongs to the plasti/stress sprint). See also
group_materi_plasti_element_group_factor.
