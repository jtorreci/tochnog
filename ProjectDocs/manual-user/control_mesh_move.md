# control_mesh_move

## Description

Moves ALL nodes of the mesh by a constant plus a linear transformation of the
current coordinates. The number of nodes and elements stays the same; only the
positions are shifted. Useful to reposition an existing mesh without remeshing.

In the x-direction a node is moved over

    c_x + l_xx*x + l_xy*y + l_xz*z

and analogously for the other axes. Only the coefficients for the number of
space dimensions of the model must be given.

## Usage

```
control_mesh_move <index> c_x l_xx l_xy l_xz c_y l_yx l_yy l_yz c_z ...
```

The coefficients are grouped per dimension. In a 1D model only the x block
(`c_x l_xx`) is read; in a 2D model the x and y blocks
(`c_x l_xx l_xy c_y l_yx l_yy`); in a 3D model all three blocks.

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `index` | Control item index (usually `0`). |
| `c_x` | Constant displacement applied to the x-coordinate of every node. |
| `l_xx` | Linear coefficient of the current x-coordinate applied to x. |
| `l_xy` | Linear coefficient of the current y-coordinate applied to x. |
| `l_xz` | Linear coefficient of the current z-coordinate applied to x. |
| `c_y` | Constant displacement applied to the y-coordinate of every node. |
| `l_yx` | Linear coefficient of the current x-coordinate applied to y. |
| `l_yy` | Linear coefficient of the current y-coordinate applied to y. |
| `l_yz` | Linear coefficient of the current z-coordinate applied to y. |
| `c_z` | Constant displacement applied to the z-coordinate of every node. |

For each dimension `d` the new coordinate is

    new_d = old_d + c_d + l_dx*x + l_dy*y + l_dz*z

## Example

Shift a 1D mesh by +1 in x:

```
number_of_space_dimensions 1
end_initia

node 1  0.
node 2  1.

control_mesh_move 0  1.  0.
end_data
```

Shift a 2D mesh by 2 in x and -1 in y:

```
control_mesh_move 0  2.  0.  0.  -1.  0.  0.
```
