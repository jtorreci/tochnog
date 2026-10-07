# group_condif_conductivity

## Description

Heat conductivity `k` of the elements of a group in the convection-diffusion
(`condif_temperature`) equation. It is a legacy GNU record (present in the
sfnet 2014 fork); this page exists to document its difference with the
Professional version.

## Usage

```
group_condif_conductivity <index> <k>
```

## Parameters

| Parameter | Meaning |
|---|---|
| `index` | Element group. |
| `k` | Isotropic heat conductivity (one value). |

Requires `condif_temperature` in the `initia` part.

## Example

```
group_type 0 -condif
group_condif_conductivity 0 0.1
group_condif_capacity 0 1.
group_condif_density 0 1.
```

## Differences with the Professional version

- **One value per direction not accepted** (known difference, not decided,
  BACKLOG C29). The Professional (manual 6.599 `kx ky kz`, test `condif11`)
  accepts an anisotropic conductivity with one value per direction; the GNU
  accepts one isotropic value only, and `group_condif_conductivity 0 0.1 0.1`
  is a parse error ("I do not know : 0.1"). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
