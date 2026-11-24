# group_groundflow_nonsaturated_eps_permeability

## Description

Floor for the relative permeability reduction in the non-saturated ground
water flow model (van Genuchten). The non-saturated law lowers the
permeability relative to the linear permeability specified in
`group_groundflow_permeability`:

```
ki = krel(S) * ksat,i
```

With this record you specify the lowest allowed value of `krel`. When `krel`
would drop below `eps`, `eps` is used instead. This avoids an (almost)
impermeable soil in strongly unsaturated zones.

## Usage

```
group_groundflow_nonsaturated_eps_permeability <element_group> <eps>
```

## Parameters

| Parameter | Meaning                                                |
|-----------|--------------------------------------------------------|
| `element_group` | Element group, see `element_group`.            |
| `eps`     | Minimum relative permeability factor.                  |

## Example

```
group_groundflow_nonsaturated_eps_permeability 0  1.e-4
```

The relative permeability never goes below `1.e-4`, even at very high suction.
