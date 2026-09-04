# element_spring_strain

Output record (one DOUBLE per element): the total spring strain, defined as the
total spring elongation (new length minus the initial length of the spring;
manual Professional 6.768). Written at every spring increment, sibling of
element_spring_force. Readable with control_print / target_item, e.g.:
```
target_item  10  -element_spring_strain 1 0
target_value 10  1. 1.e-8
```
