# control_mesh_mirror

## Description

Mirrors the mesh about the plane x=0, y=0 or z=0. Every node is duplicated with
the chosen coordinate negated, and every element is duplicated with the new
nodes, so the model ends with the double of nodes and elements. Useful to
generate symmetric meshes from a half model without remeshing.

## Usage

```
control_mesh_mirror <index> -x|-y|-z
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `-x` | Mirror about the plane x=0. |
| `-y` | Mirror about the plane y=0. |
| `-z` | Mirror about the plane z=0. |

## Example

Mirror a 2D half model about x=0:

```
number_of_space_dimensions 2
end_initia

node 1  0.  0.
node 2  1.  0.
node 3  1.  1.
element 1  -bar2 1 2

control_mesh_mirror 0  -x
end_data
```

After this record the model has 6 nodes and 2 elements, symmetric about x=0.
