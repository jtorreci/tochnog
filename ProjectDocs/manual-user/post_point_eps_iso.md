# post_point_eps_iso

## Description

`post_point_eps_iso index eps` (manual Professional 6.957): the tolerance on the isoparametric coordinates with which `post_point` (same index) is accepted in an element. Default 1.e-3; increase it (say 0.1) for a point on or just over the border of the mesh.

## Usage

```
post_point 3 3.0015 0.5
post_point_eps_iso 3 0.1
```

## Example

`validation-suite/test-2026/post_point_search1.dat`: a point 0.0015 beyond the mesh is found with 0.1 and not with the default.

## Differences with the Professional version

The Professional 02-08-2026 accepts the record but does not find the point with 0.1 either; the GNU follows the manual, a deliberate divergence decided by the owner on 2026-11-05 (DIVERGENCES POST-POINT-SEARCH, BACKLOG C111).
