# materi_plasti_hypo_history

## Implementation

- **Enum**: `MATERI_PLASTI_HYPO_HISTORY` in tochnog.h and tochnog-mod.h
  (in sync; next to `MATERI_PLASTI_CAP1_HISTORY` /
  `MATERI_PLASTI_DIPRISCO_HISTORY`).
- **Globals** (initia.cc): flag `materi_plasti_hypo_history`; the dof
  index is the shared `hisv_indx` (same dof as the generic history
  variables).
- **Parser** (input.cc, branch shared with `materi_history_variables`
  and `materi_plasti_diprisco_history`): when the keyword is
  `materi_plasti_hypo_history`, `materi_history_variables` is forced to
  8 and the dof type is `-MATERI_PLASTI_HYPO_HISTORY` with
  `8*nder` slots:
  ```c
  materi_history_variables = 8;
  hisv_indx = unknown_indx;
  n = materi_history_variables;
  array_set( &dof_type[hisv_indx], -MATERI_PLASTI_HYPO_HISTORY, n*nder );
  array_set( &dof_scal_vec_mat[unknown_indx], -SCALAR, n*nder );
  ```
- **Registration** (database.cc): name[] of the initia + basename
  `hyhis` with counter 0..7 in the DOF_LABEL table (the `hyhis` block
  restarts the counter at `hisv_indx`, pattern of the `hisv` block).
- **Unknown handling** (general.cc): the `hyhis` unknowns are integrated
  like the generic history variables (`unknown_belongs_to_type = 1`,
  `inertia = 1.`, `conv_part = 1.` for MATERI type).
- **Self-checks** (check.cc): the hypoplasticity group records
  (`GROUP_MATERI_PLASTI_HYPO_*`) now accept `materi_plasti_hypo_history`
  OR `materi_history_variables`
  (`check_unknown_atleastone`).
- **Kernel mapping**: the Masin branch (hypoplas.cc) reads the void
  ratio from `new_hisv[e_slot]` and the sensitivity from
  `new_hisv[s_slot]` with `e_slot/s_slot = 0/4` when
  `materi_plasti_hypo_history` is set (Professional layout, manual
  4.23) and `6/7` for the legacy layout. The intergranular strain delta
  (statev[0..5] of the Masin kernels) comes from the epi dof
  (`materi_strain_intergranular`) in the hyhis layout, not from the
  history slots.

## Verification

- hypo2/hypo4 of the corpus: rc=0 (targets of the Professional).
- hypo7/8/9/12/13 run and converge close to the Professional targets
  (3-6%); the residual difference is the accuracy of the masin.c kernel
  port (see the pending calibration note in SEGUIMIENTO).
- The reset `-hyhis0`/`-hyhis4` in the corpus tests resolves to the
  expected slots (verified against the Professional .dbs: e evolves
  from the reset value, structure s stays at the reset value).

## Pending

- The Niemunis visco models (`group_materi_plasti_hypo_niemunis_visco`,
  `group_materi_plasti_hypo_wolffersdorff_niemunis`) are not
  implemented (no kernel term in the GNU hypo.c); hypo6/hypo11 stay
  RUNFAIL.
- The ISA extension (`materi_strain_isa_c`, `materi_strain_isa_eacc`,
  `group_materi_plasti_hypo_strain_isa`) is not implemented; hypo10
  stays PARSE.
