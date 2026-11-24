# check_element_node

## Descripción

`check_element_node` validates the mesh by checking that no element references
a node more than once. If an element has a duplicate node, the run aborts with
an error identifying the offending element. The check is active by default.
Use `check_element_node -no` to disable it when duplicate nodes are intended
(for example in an intentionally degenerate or periodic mesh).

## Uso

Place it in the data part, as a keyword line:

```
check_element_node -no
```

## Parámetros

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `-yes`    | Enable the duplicate-node check (default). Aborts on duplicates. |
| `-no`     | Disable the check; duplicate nodes are accepted.               |

## Ejemplo

Minimal input that disables the duplicate-node check:

```
control_geometry
   cartesian
control_time
   end 1.0
   dt 0.1
control_print
   history 0
   step 10
materi_elasti_young
   0 2.e7
materi_elasti_poisson
   0 0.3
materi_density
   0 2500.
check_element_node
   -no
```

The check runs after `end_data` closes the data part. If a duplicate node is
found (with the check enabled) the run stops with an error like:

```
Error: element <ielem> has duplicate nodes.
```
