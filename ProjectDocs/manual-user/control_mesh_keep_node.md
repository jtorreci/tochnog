# control_mesh_keep_node

## Description

Deletes ALL the nodes of the mesh EXCEPT the ones listed. Useful to visualize
only a specific set of nodes (and their elements) in isolation.

NOTE: elements that reference deleted nodes will no longer be valid and will
fail. Use this record for visualization only.

## Usage

```
control_mesh_keep_node <index> node_0 [node_1 ...]
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `node_n` | Number of the node to keep. |

## Example

Keep only nodes 1 and 2, deleting all the others:

```
control_mesh_keep_node 0 1 2
```
