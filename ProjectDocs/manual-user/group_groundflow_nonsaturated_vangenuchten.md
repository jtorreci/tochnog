# group_groundflow_nonsaturated_vangenuchten

## Description

Parameters of the van Genuchten non-saturated ground water flow model. The
saturation is a function of the pore-pressure head:

```
phi_p = -pres / (density * |g|)
S(phi_p) = Sresidu + (Ssat - Sresidu) * (1 + (ga*|phi_p|)^gn) ^ ((1-gn)/gn)
```

The total capacity is the saturated capacity plus a non-saturated part, and the
total permeabilities are a relative factor of the saturated ones:

```
c  = csat + n * dS/dphi_p
ki = krel(S) * ksat,i
```

with the effective saturation and the relative permeability (Mualem):

```
Se   = (S - Sresidu) / (Ssat - Sresidu)
krel = Se^gl * (1 - (1 - Se^(gn/(gn-1))) ^ ((gn-1)/gn)) ^ 2
```

Requires the saturated `group_groundflow_capacity` and
`group_groundflow_permeability` as usual, the porosity in
`group_groundflow_porosity`, and `groundflow_saturation` initialised in the
initialisation part.

## Usage

```
group_groundflow_nonsaturated_vangenuchten <element_group>
    <Sresidu> <Ssat> <ga> <gl> <gn>
```

## Parameters

| Parameter   | Meaning                                                     |
|-------------|-------------------------------------------------------------|
| `element_group` | Element group, see `element_group`.                    |
| `Sresidu`   | Residual saturation.                                         |
| `Ssat`      | Saturated saturation (normally 1.0; less with trapped air). |
| `ga`        | van Genuchten constant `ga`.                                 |
| `gl`        | van Genuchten constant `gl` (exponent of the relative permeability). |
| `gn`        | van Genuchten constant `gn`.                                 |

## Example

```
group_groundflow_nonsaturated_vangenuchten 0
                        0.1
                        1.0
                        0.1
                        1.0
                        2.0
```

Since the van Genuchten law is highly non-linear, convergence may be hard;
check `post_node_rhside_ratio` and consider including inertia. Use
`groundflow_nonsaturated_apply -no` or
`control_groundflow_nonsaturated_apply -no` to run saturated only.
