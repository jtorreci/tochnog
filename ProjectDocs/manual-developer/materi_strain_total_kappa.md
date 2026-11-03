# materi_strain_total_kappa

## Implementation (F5 T2 V7)

- **Input:** `input.cc` parses the word, sets the flag `materi_strain_total_kappa` and the index
  `eptk_indx`. It adds one SCALAR dof of type `-MATERI_STRAIN_TOTAL_KAPPA`. A check
  after the initia part requires `materi_strain_total`
  (check_unknown ... CHECK_USAGE_AND_ERROR).
- **Declarations:** the enum and globals are in `tochnog.h` /
  `tochnog-mod.h` / `initia.cc`.
- **Labels:** `database.cc` gives the label `eptk` and the name.
- **Fill:** `materi.cc`, "Total strain kappas" block before the C68
  integration-point states. The measures of the symmetric `new_ept` are
  computed: the size, the deviatoric size, and the principal values from
  `matrix_eigenvalues`. Then `new = max(old, measure)` and the block adds
  `volume * w_rec * (new - old) / dtime`.
- **Integration-point state:** the dof is an integration-point state.
  `elem.cc` masks it in both places, so it is read back from and written
  to ELEMENT_INTPNT_DOF, and the history maximum is per point. Its nodal
  value uses `strain_recovery_weight` (`general.cc` recovery list).
- **Hardcoded:** the size measure `sqrt(0.5 a:a)`, the same as kappa in
  `stress.cc`; the compression kappa is stored as a positive number.

## Verification

Analytic: `materi_strain_total_kappa1` (homogeneous strain, peak then
unloading) and `materi_strain_total_kappa_error1` (missing
materi_strain_total).

## Pending

- No consumer yet: the manual's softening through `dependency_diagram` on
  these dofs works only if `dependency_item` accepts them as a dof
  (not verified).
