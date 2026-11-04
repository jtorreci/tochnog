# post_calcul_label

## Implementation (F5 T2 V8)

- **Record:** `database.cc`, INTEGER, variable length, no index, class
  POST. It is never stored: the database holds numbers and enum names only.
- **Input:** `input.cc` (`input()`, data part). For this record every
  value is read and ignored (`ival[iv] = 0`, no enum lookup), and the store
  loop is skipped (`ready=(idat==POST_CALCUL_LABEL)`). A Professional
  `.dbs` with the record can therefore be read.
- **Output:** `print_db.cc`, `print_database()`, before `end_data` of a
  full database print (`-yes`, `-everything`, last database): when
  `POST_CALCUL_SCAL_VEC_MAT` exists, it writes `post_calcul_label` with
  `post_calcul_names[0..ncalcul-1]` (`calcul.cc`, `calculate()`).

## Pending

- The GNU generic names (`aept`) differ from the Professional `.dbs`
  (`a_ept`); renaming them would change the vtk / gmsh / node print names
  and `target_item` labels (DIVERGENCES POST-CALCUL-LABEL-NAMES).
- A single `control_print -post_calcul_label` to the screen is not
  supported (the record is only in the `.dbs`).

## Verification

`post_calcul_modify1`: the `.dbs` line is checked in
`scripts/build_safe.sh`; the input copy of the Professional record runs.
