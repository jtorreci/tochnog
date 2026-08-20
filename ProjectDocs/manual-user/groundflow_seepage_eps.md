# groundflow_seepage_eps

## Description

Tolerance for the groundflow seepage condition. If the inner product of the
groundflow water flow direction with the normal outside the material is
smaller than `eps`, the seepage status is set to closed, and the total pressure
condition is not applied (so that the boundary is really closed for water
flow). If not specified, `eps` is set to `0.1`.

## Usage

```
groundflow_seepage_eps <eps>
```

## Parameters

| Parameter | Meaning                                                         |
|-----------|-----------------------------------------------------------------|
| `eps`     | Tolerance on the outward flow. Default `0.1`.                   |

## Example

```
groundflow_seepage_eps 0.05
```

The seepage edge is closed unless the outward water flow exceeds the tolerance
0.05.
