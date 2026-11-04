# materi_strain_plasti_mohr_coul

## Description

`materi_strain_plasti_mohr_coul` (manual Professional 4.41) adds the plastic
strain `eps_kl_plas` **specifically for the Mohr-Coulomb model** to the
`node_dof` records. It uses the 6-component tensor layout of
[`materi_strain_plasti`](materi_strain_plasti.md) with its own dof
(labels `eppmoxx eppmoxy eppmoxz eppmoyy eppmoyz eppmozz`).

The dof is filled only in element groups that use the Mohr-Coulomb model
(`group_materi_plasti_mohr_coul` or `group_materi_plasti_mohr_coul_direct`). In those groups it equals `materi_strain_plasti`; in all
other groups it stays 0.

## Usage

In the initialization part:

```
materi_strain_plasti_mohr_coul
```

## Parameters

None.

## Example

`validation-suite/test-2026/materi_strain_plasti_model2.dat`: a Mohr-Coulomb model with
`materi_strain_plasti` and the four per-model dofs. The mohr_coul dof equals the
`epp` columns at every node, and the other three stay 0.

## Differences with the Professional version

The Professional 02-08-2026 accepts `materi_strain_plasti_vonmises`,
`_mohr_coul` and `_tension` but adds no dof, and it rejects
`materi_strain_plasti_camclay` ("Error in initialisation part"). GNU follows
the manual (owner decision 2026-11-03, DIVERGENCES STRAIN-PLASTI-MODEL-DOFS).
When a group combines several models (for example Mohr-Coulomb with a
tension cut-off), every active model's dof receives the total plastic
strain of the point. GNU does not split the strain per yield surface.
