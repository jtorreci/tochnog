# group_materi_elasti_young_power

## Description

`group_materi_elasti_young_power` (manual Professional 6.662, theory
section 2.2.2) makes the Young modulus a **power law of the pressure
state**:

    E = E0 + E1 * (p/p1)^alpha,      with the conditions  E >= E2  and  E <= E3

where `p` is the pressure, `p = -(sig11 + sig22 + sig33)/3` (positive in
compression: a compressive state raises E). The record takes **6
parameters**: `E0 E1 E2 E3 p1 alpha`.

> **Semantic change (Sprint 10, lote 6)**: this is the Professional
> convention. The old GNU form of the record (3 parameters,
> `young = young0 * |p/p0|^alpha`) is **no longer accepted** — inputs
> written for the GNU 3-parameter form must be converted to the
> 6-parameter Professional form.

The record combines with `group_materi_elasti_poisson` (the poisson
coefficient is taken from that record) and with
`group_materi_elasti_poisson_power`, `group_materi_elasti_shear_factor`,
`group_materi_elasti_stress_pressure_history_factor` and
`control_materi_elasti_k0`. When `group_materi_elasti_young_power` is
present the power law **replaces** the constant
`group_materi_elasti_young` (it does not add to it — see developer
manual for the GNU accumulation gotcha).

## Uso

```
group_type 0  -materi
group_materi_elasti_young 0  1000.0   ( fallback, overridden by the power law )
group_materi_elasti_poisson 0  0.3
group_materi_elasti_young_power 0
                        1000.   ( E0 )
                        500.    ( E1 )
                        800.    ( E2: E >= E2 )
                        3000.   ( E3: E <= E3 )
                        1.      ( p1, must be > 0 )
                        1.      ( alpha )
```

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_elasti_young_power` | `E0 E1 E2 E3 p1 alpha` | Power-law Young modulus: `E = E0 + E1*(p/p1)^alpha` with the caps `E >= E2` and `E <= E3`. `p1` must be positive; `p = -sig_mean` (positive in compression). |

## Física

The pressure-dependent modulus follows the Professional theory
(section 2.2.2): in a confined compression (oedometer) the modulus
grows with the mean compressive stress, which makes the lateral stress
growth superlinear with the vertical strain. The law is evaluated with
the CURRENT stress state each time the elastic stiffness is computed.

## Validation

- `myoung6.dat`: oedometer with `E0 = 1000`, `E1 = 500`, `E2 = 800`,
  `E3 = 3000`, `p1 = 1`, `alpha = 1`, `nu = 0.3`, `eps_zz = -0.001`.
  The self-consistent analytic fixed point (`E = 1000 + 500*p`,
  `p = E*eps/(3(1-2*nu))` -> `E = 1714.29`, `p = 1.4286`) gives
  `sigma_zz = -2.3077`, `sigma_xx = -0.9890`. The incremental code
  lands at `E_eff ~ 1237` (72% of the fixed point): measured
  `sigma_xx = -0.713886`, `sigma_yy = -1.66573` (window targets that
  exclude the linear base `-0.5769/-1.3462`). See the developer manual
  for the deviation GOTCHA.
- `myoung6_e2.dat`: `E2 = 3000` forces `E = 3000` exactly:
  `sigma_xx = -1.7308`, `sigma_yy = -4.0385` (analytic, cap active).
- `myoung6_e3.dat`: `E3 = 800` forces `E = 800` exactly:
  `sigma_xx = -0.4615`, `sigma_yy = -1.0769` (analytic, cap active).
- `myoung6_apply.dat`: with `materi_elasti_young_power_apply -no` the
  constant young `E0 = 1000` is applied at all times:
  `sigma_xx = -0.5769`, `sigma_yy = -1.3462` (the linear base).

## Notas

- Requires `materi_stress`, `materi_velocity` and `materi_strain_total`
  in the initialization part (same check as the GNU record).
- The base `group_materi_elasti_young` record is read as a fallback and
  is replaced when the power record exists.
