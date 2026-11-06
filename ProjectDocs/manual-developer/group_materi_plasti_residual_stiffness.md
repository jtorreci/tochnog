# group_materi_plasti_residual_stiffness (developer)

## Where

- `database.cc`: record `GROUP_MATERI_PLASTI_RESIDUAL_STIFFNESS`
  (DOUBLE, length 1, MATERI, `GROUP_TYPE`); enum appended in `tochnog.h` /
  `tochnog-mod.h`.
- `stress.cc` `set_stress()`: the flag `plastic_ip` is set when the first
  yield test starts plastic iterations (`plasti_iter==1`, `f>=EPS_F`) and
  when the K13 direct cut stores a plastic strain. Before the material
  matrix is filled, the factor is read with `get_group_data`
  (`GET_IF_EXISTS`); 0 or less stops with
  `Error: group_materi_plasti_residual_stiffness must be larger than 0.`

## Implementation

With C = Cmem + Cuser + Chyper + Chypo and S = C_mmnn, the volumetric
part is V_ijkl = delta_ij delta_kl S/9 (K delta delta for isotropic C).
A plastic point gets C' = V + factor (C - V) in `ddsdde`; the
compressibility term is volumetric and not scaled. The stresses and the
right-hand side are untouched, so only the iteration path changes.

The idea matches the GNU residual stiffness of `slide_plasti_*`,
`group_interface_materi_*` and `support_edge_normal_plasti_*` (matrix only).

## Verification

Probes on the `plasti_apply1` square (Mohr-Coulomb, total_linear):
50 iterations: disx 0.0182467717 with and without factor 0.1; one
iteration: sigxx 0.0691216 (no record, = Professional) and -0.0101589
(factor 0.1; Professional 0.0691216). Tests
`group_materi_plasti_residual_stiffness1`-`3`, `_error1`.

## Pending

- "Plastic element" is read per integration point in the current update;
  a hypoplastic point is never plastic.
- For an anisotropic C the coupling terms are scaled with the deviatoric
  part.
