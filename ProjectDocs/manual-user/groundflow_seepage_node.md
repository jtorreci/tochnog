# groundflow_seepage_node

## Description

Does the same as the `groundflow_seepage_geometry` record, but now you specify
node numbers at which the seepage condition holds.

```
groundflow_seepage_node <index> <node_0> <node_1> ...
```

`node_0` is the first node number, `node_1` the second, etc.

## Parameters

| Parameter | Meaning                                                         |
|-----------|-----------------------------------------------------------------|
| `index`   | Record index.                                                    |
| `node_0 ...` | Node numbers where the seepage condition holds.             |

## Example

```
groundflow_seepage_node 0  1 2
bounda_dof 20  -ra 1 2 -ra -pres
bounda_time 20 0.0
```

Nodes 1 and 2 are a seepage edge: water can only leave the domain there.
