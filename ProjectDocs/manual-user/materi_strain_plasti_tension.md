# materi_strain_plasti_tension

## Description

`materi_strain_plasti_tension` (manual Professional 4.42) adds the plastic
strain `eps_kl_plas` **specifically for the tension model** to the
`node_dof` records. It uses the 6-component tensor layout of
[`materi_strain_plasti`](materi_strain_plasti.md) with its own dof
(labels `epptexx epptexy epptexz eppteyy eppteyz epptezz`).

The dof is filled only in element groups that use the tension model
(`group_materi_plasti_tension` or `group_materi_plasti_tension_direct`). In those groups it equals `materi_strain_plasti`; in all
other groups it stays 0.

## Usage

In the initialization part:

```
materi_strain_plasti_tension
```

## Parameters

None.

## Example

`validation-suite/test-2026/materi_strain_plasti_model4.dat`: a tension model with
`materi_strain_plasti` and the four per-model dofs. The tension dof equals the
`epp` columns at every node, and the other three stay 0.

## Differences with the Professional version

The Professional 02-08-2026 adds the same dofs for `materi_strain_plasti_vonmises`,
`_mohr_coul` and `_tension`, with the same labels and the same values (0
under another model), and it rejects `materi_strain_plasti_camclay`
("Error in initialisation part"). GNU keeps the camclay dof as the manual
says (DIVERGENCES STRAIN-PLASTI-MODEL-DOFS; the V7 reading "adds no dof"
was wrong, BACKLOG C108).
When a group combines several models (for example Mohr-Coulomb with a
tension cut-off), every active model's dof receives the total plastic
strain of the point. GNU does not split the strain per yield surface.
