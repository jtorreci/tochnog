# group_condif_conductivity

## Implementation

- `database.cc`: `GROUP_CONDIF_CONDUCTIVITY`, `DOUBLE_PRECISION`,
  `data_length` 1, class `CONDIF`, required `GROUP_TYPE`.
- `check.cc`: requires the `condif_temperature` unknown.
- `general.cc` (`type == -CONDIF` branch): read with `get_group_data()`
  into the scalar `condif_conductivity` (`GET_IF_EXISTS`).

Legacy GNU record (sfnet 2014); no change in this convergence cycle.

## Hardcoded parameters

- Length 1: the conductivity is a scalar (isotropic).

## Differences with the Professional version

- **One value per direction not accepted** (known difference, not decided,
  BACKLOG C29). The Professional manual 6.599 documents
  `group_condif_conductivity index kx ky kz` and the Professional test
  `condif11` uses `0.1 0.1`; the GNU `data_length` 1 makes it a parse error.
  Proposed fix: accept ndim values (isotropic when one value is given), as
  `group_groundflow_permeability`. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).

## Pending refactorings

- Anisotropic conductivity (C29).
