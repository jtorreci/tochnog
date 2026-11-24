# post_point_element_group

## Description

`post_point_element_group index element_group` (manual Professional 6.954): the element holding `post_point` with the same index is searched only in this element group. A point outside the group's elements is not found (no `post_point_dof`).

## Usage

```
post_point 1 2.5 0.5
post_point_element_group 1 0
```

## Example

`validation-suite/test-2026/post_point_search1.dat`.

## Differences with the Professional version

Same search result. When a point is not found the Professional writes a `post_point_dof` of zeros; the GNU writes none (pre-existing).
