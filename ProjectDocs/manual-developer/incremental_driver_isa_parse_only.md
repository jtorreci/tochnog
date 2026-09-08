# incremental_driver ISA records — parse-only boundary

## Approved scope

GNU recognizes and stores only these Professional initialisation records:

- `materi_strain_isa_c`: six symmetric matrix components.
- `materi_strain_isa_eacc`: one scalar component.

This is syntax/DOF metadata support, not ISA or Niemunis constitutive support.
No material implementation, dispatch, state evolution, or incremental-driver
behavior reads these records.

## Professional layout measured

The Professional `incremental_driver_syntax.dat` declares the records after
`materi_strain_intergranular`. Its emitted `.dbs` shows the following order and
metadata:

| Record | Components | `dof_scal_vec_mat` | Professional labels |
| --- | ---: | --- | --- |
| `materi_strain_isa_c` | 6 | `-matrix` | `-episa_cxx`, `-episa_cxy`, `-episa_cxz`, `-episa_cyy`, `-episa_cyz`, `-episa_czz` |
| `materi_strain_isa_eacc` | 1 | `-scalar` | `-episa_eacc` |

Evidence inspected:

- Input: `/tmp/opencode/tn_prof/tochnog_version_02-08-2026/test/incremental_driver/incremental_driver_syntax.dat`.
- Professional parser output: `/tmp/incdrv_validation_clean/professional/incremental_driver_syntax/incremental_driver_syntax.dbs`.

## Unsupported execution

An input that declares these records can parse, and the stored DOFs remain
unused by GNU's incremental driver. Any result that requires ISA/Niemunis
constitutive evolution remains unsupported; no substitute model is inferred or
provided by this registration.
