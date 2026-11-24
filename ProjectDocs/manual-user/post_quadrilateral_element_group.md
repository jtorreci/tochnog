# post_quadrilateral_element_group

## Description

`post_quadrilateral_element_group index element_group` (manual Professional 6.963): the dof values of `post_quadrilateral` (same index) are taken from this element group only; the points of the quadrilateral outside the group are skipped and the average is over the found points.

## Usage

```
post_quadrilateral 1 1 0 3 0 3 1 1 1
post_quadrilateral_element_group 1 0
```

## Example

`validation-suite/test-2026/post_point_search1.dat`: velx 0.0880956.

## Differences with the Professional version

None measured (same value).
