# group_materi_plasti_hypo_intergranularstrain (6 parameters)

## Implementation

- **Registration** (database.cc):
  `GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN` data_length changed
  from 5 to 6 (manual Professional 6.715: R m_R m_T beta_r chi theta).
- **Kernel buffer** (hypoplas.cc): `LENGTH_INTERGRANULARSTRAIN` 5 -> 6
  (the local array is read with GET_AND_CHECK, so 6 values are stored).
- **Kernel interface**: the hypo_ kernel signature carries five
  intergranular parameters (epi_R, epi_mr, epi_mt, epi_betar, epi_chi);
  the sixth value (theta) has no dedicated slot: the kernel evaluates
  the `rho^theta f_d N S_hat` term with `chi`, which the manual
  recommends to equal theta for monotonic loading.
- **Alias** (database.cc, db_number):
  `group_materi_plasti_hypo_strain_intergranular` (Professional
  spelling) -> GROUP_MATERI_PLASTI_HYPO_INTERGRANULARSTRAIN.

## Verification

- hypo2/hypo3 of the corpus (6 values) pass their targets.
- hypo10 (ISA extension) stays PARSE: it needs
  `materi_strain_isa_c`/`materi_strain_isa_eacc` and
  `group_materi_plasti_hypo_strain_isa`, not implemented.

## Pending

- theta as an independent exponent (a separate kernel term) is not
  implemented; documented partial.
