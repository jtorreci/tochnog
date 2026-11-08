# mesh_correct

## Description

Checks and corrects the node order of interface elements (Professional
manual 6.816). With `-yes`, a `-quad4` or `-hex8` of a `group_interface`
group whose two sides are given along the wrong isoparametric direction
is renumbered so that the interface opens along eta (quad4) or zeta
(hex8), the direction the interface element reads its two sides from.

## Usage

```
mesh_correct switch
```

`switch` is `-yes` or `-no` (default `-no`). The record has no index and
is given in the data part.

## What is corrected

The opening direction is the isoparametric direction with the shortest
distance from the first node (zero for coincident nodes). When it is not
the interface direction, the nodes are reordered as the Professional
02-08-2026 does:

| element | opening along | new order |
|---------|---------------|-----------|
| quad4 `n1 n2 n3 n4` | xi | `n1 n3 n2 n4` |
| hex8 `n1 .. n8` | xi | `n2 n4 n6 n8 n1 n3 n5 n7` |
| hex8 `n1 .. n8` | eta | `n1 n2 n5 n6 n3 n4 n7 n8` |

An element with no strictly shorter direction (for example a quad4 whose
second side is reversed, `3 4 6 5`) is left as it is, also by the
Professional.

## Example

```
element 3 -quad4 3 5 4 6
element_group 3 1
group_interface 1 -yes
mesh_correct -yes
```

Nodes 3 and 5 coincide, so the element becomes `-quad4 3 4 5 6`. Tests:
`validation-suite/test-2026/mesh_correct1.dat` .. `4.dat`,
`mesh_correct_off1.dat`.

## Differences with the Professional version

None found on the probes. Quadratic interface elements (`-quad9`,
`-hex27`) are not corrected.
