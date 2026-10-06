# direct_plastic_strain — plastic strain and yield value of the direct cut-offs

## What it does

The direct plasticity laws (`group_materi_plasti_tension_direct`,
`group_materi_plasti_mohr_coul_direct`, `group_materi_plasti_compression_direct`,
full principal mode and plane mode with `_normal` / `_normal_automatic`,
also with `_visco` and `_wall`) cut the stress directly. Since 2026-10-06
(audit K13) the cut is also written to the state output, as the
Professional does:

- `materi_strain_plasti` (`eppxx` ...): the plastic strain that keeps
  sigma = C (eps - eps_p), i.e. the increment `C^-1 (sigma_trial -
  sigma)` with C the elastic stiffness of the group. `materi_strain_elasti`
  (`epe..`) is reduced by the same amount.
- `materi_plasti_f` (`f`): the largest yield value of the direct laws on
  the final stress: tension `sig_max - sigy` (or `sig_n - sigy` on the
  tension plane), Mohr-Coulomb `0.5(s1-s3) + 0.5(s1+s3) sin(phi) - c
  cos(phi)` (or `|tau| - (c - sig_n tan(phi))` on the MC plane),
  compression `sigy - sig_min`; `0` when a law cut the stress (also a
  partial `_visco` cut).

Before, the GNU left `eppxx` = 0 and `f` = -1e10 (the "no yield law"
value), so the output violated sigma = E (eps - eps_p).

## Why it is useful

Plastic-strain contour plots, plastic-zone detection (`f` = 0) and
post-processing that rebuild the elastic strain now work for the direct
laws as for the incremental ones. The stress itself is unchanged.

## Input

Nothing new: declare the dofs in the initia part.

```
materi_strain_plasti
materi_plasti_f
```

## Example

Single quad4, plane strain, E 1, nu 0.25, uniaxial eps_xx 2e-3, full
tension cut-off sigy 1e-3:

```
group_materi_plasti_tension_direct 0 1.e-3
```

gives sigxx 1e-3, eppxx 1.4e-3, eppyy = eppzz -0.35e-3, epexx 0.6e-3,
f 0 (`validation-suite/test-2026/direct_plastic_strain_tension1.dat`; MC
and elastic-state cases: `direct_plastic_strain_mc1`,
`direct_yield_value1`). `visc_pl4` (Professional test, 1D tension
visco): eppxx 0.5 = eps - sigma/E.

## Differences with the Professional version

- **Shear components (deliberate).** The Professional writes the shear
  plastic strain of a direct cut as `dsig_xy / G`, twice the tensor
  value, although its own total and incremental-law plastic strains are
  tensor components (von Mises shear: both binaries store the tensor
  value). The GNU stores the tensor value `dsig_xy / 2G`, so sigma = C
  (eps - eps_p) holds component by component. Example
  (`direct_plastic_strain_mc1`): GNU eppxy 2.04114892e-3, Professional
  4.08229785e-3; every normal component and `f` are equal.
- **Visco stress path.** With `_visco` the GNU cut differs from the
  Professional (BACKLOG C33: the relaxation factor uses the previous
  time step, and the principal tension cap ignores `_visco`); the stored
  plastic strain follows each binary's own stress.
