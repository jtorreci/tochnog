# post_point_node_type

## Description

`post_point_node_type index node_type` (manual Professional 6.959): which node coordinates locate `post_point` (same index). `-node`: the current coordinates, so in an updated calculation you see the material that is at the point now. `-node_start_refined` (default): the initial coordinates, the material that was there at the start.

## Usage

```
post_point 5 3.05 0.5
post_point_node_type 5 -node
```

## Example

`validation-suite/test-2026/post_point_node_type1.dat` (mesh -follow_material, the right edge moves from 3 to 3.1).

## Differences with the Professional version

The Professional 02-08-2026 gives the `-node_start_refined` result for `-node` too; the GNU follows the manual (DIVERGENCES POST-POINT-SEARCH).
