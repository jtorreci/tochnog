# group_materi_maxwell_chain

## Description

Linear viscoelastic material built from `n` parallel Maxwell chains. Each chain
`m` has a stiffness `E_m` and a relaxation time `t_m`. The number of chains `n`
must equal `materi_maxwell_stress` set in the initialization phase.

At each time step, every chain updates its viscous (Maxwell) stress
`new_msig[m]` from the previous rotated value `rotated_old_msig[m]`, the
strain increment `inc_epe` and the relaxation factor `1-exp(-dtime/t_m)`:

```
inc_msig = C(E_m,nu) . inc_epe * t_m/dtime - rotated_old_msig[m]
new_msig[m] = rotated_old_msig[m] + inc_msig * (1 - exp(-dtime/t_m))
```

In a `TOTAL` formulation the chain stresses are added to the total stress; in
an `INCREMENTAL` formulation only the increment is added. The viscoelastic
tangent is assembled from the same stiffness `C(E_m,nu)` of each chain.

Useful for creep/relaxation of polymers, bituminous and sealing materials, or
creeping geosynthetics, where a few parallel Maxwell chains approximate the
relaxation spectrum.

Requires `materi_stress`, `materi_maxwell_stress` and a velocity/displacement
formulation with `group_materi_memory` (e.g. `-updated_without_rotation`).

## Usage

```
materi_maxwell_stress <n>

group_materi_maxwell_chain <element_group>  E_0 t_0 E_1 t_1 ... E_{n-1} t_{n-1}
```

## Parameters

| Parameter | Meaning |
|-----------|---------|
| `E_m` | Stiffness (Young-like modulus) of chain `m`. |
| `t_m` | Relaxation time of chain `m`. |
| `materi_maxwell_stress` | Number of Maxwell chains `n` (initialization). It equals `group_materi_maxwell_chain` data length / 2. |

`group_materi_maxwell_chain` must be used together with
`group_materi_elasti_poisson` (the Poisson ratio is shared by all chains).

## Example

```
materi_velocity
materi_strain_total
materi_maxwell_stress 1
materi_stress
end_initia
...
group_materi_memory 0 -updated_without_rotation
group_materi_maxwell_chain 0  1.  1.
```

This is the `viscel1.dat` regression test: a single Maxwell chain with
`E_0=1, t_0=1`.

## Notes

- `group_materi_maxwell_chain_nonlinear` is REGISTERED but its model routine
  `visco_elastiticity_nonlinear()` in `visconon.cc` is an EMPTY stub — it does
  nothing. Do not use it; only the linear `group_materi_maxwell_chain` is
  functional.
