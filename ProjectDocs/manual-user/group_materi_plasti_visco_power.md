# group_materi_plasti_visco_power

## Description

Viscoplasticity data for the power (creep) model (manual Professional
6.748). The record is used in combination with a plasticity model (e.g.
`group_materi_plasti_tension`, `group_materi_plasti_mohr_coul_direct`):

```
group_materi_plasti_visco_power <index> <eta> <p>
```

- `<index>` is the element group (see `group_type`).
- `<eta>` is the fluidity constant of the power model.
- `<p>` is the power.

## Theory

The viscoplastic strain rate of the power model (manual Professional,
section 2.2.11) reads

```
eps_dot_kl^plas = eta * (f)^p * d f_flow / d sigma_kl
```

where `f` is the yield function of the active plasticity model and
`f_flow` its plastic potential. With `f` in the natural stress units of
the problem, the record carries NO reference stress: `f` enters the law
as-is (an implicit reference of 1.0).

The exponential sibling law (`group_materi_plasti_visco_exponential`,
section 2.2.11) reads `eps_dot_kl^plas = gamma * p * exp(alpha*f) *
d f_flow / d sigma_kl` with `p` the (compressive) pressure.

## Legacy GNU layout

The record also still accepts the legacy GNU layout with a reference
stress as third value:

```
group_materi_plasti_visco_power <index> <eta> <p> <f_ref>
```

with `eps_dot_pl = eta * (f/f_ref)^p`. The Professional tests of the
sfnet corpus (`visc_pl1`, `visc_pl2`, `validation_12`) use the two-value
Professional form.

## Example

```
group_materi_plasti_tension 0 1.0
group_materi_plasti_visco_power 0 1.e-1 1.0
```

A 1-D tension bar with fluidity 0.1: the stress overshoots the tensile
strength while loading and relaxes back to the yield surface on the
viscoplastic time scale `1/eta` (Perzyna behaviour; `visc_pl1` of the
corpus).
