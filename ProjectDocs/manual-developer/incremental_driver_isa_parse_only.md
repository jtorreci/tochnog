# incremental_driver ISA records — parse-only boundary

## Approved scope

GNU recognizes and stores only these Professional ISA metadata records:

- `materi_strain_isa_c`: six symmetric matrix components.
- `materi_strain_isa_eacc`: one scalar component.
- `group_materi_plasti_hypo_strain_isa`: indexed MATERI data with two doubles.

> **Limitation — not ISA/Niemunis support.** These records do not provide a
> material implementation, dispatch, state evolution, or computed hypo result.
> They must never masquerade as simulation output.

This is syntax/metadata support. When an active
`group_materi_plasti_hypo_strain_isa` record reaches `incremental_driver`, GNU
uses the narrow validation-only exit documented below.

## Professional layout measured

The Professional `incremental_driver_syntax.dat` declares the records after
`materi_strain_intergranular`. Its emitted `.dbs` shows the following order and
metadata:

| Record | Components | `dof_scal_vec_mat` | Professional labels |
| --- | ---: | --- | --- |
| `materi_strain_isa_c` | 6 | `-matrix` | `-episa_cxx`, `-episa_cxy`, `-episa_cxz`, `-episa_cyy`, `-episa_cyz`, `-episa_czz` |
| `materi_strain_isa_eacc` | 1 | `-scalar` | `-episa_eacc` |
| `group_materi_plasti_hypo_strain_isa` | 2 doubles | indexed `MATERI` | fixture: `0 20.0 0.017` |

The Professional `.dbs` confirms the two-value layout. Their names and
constitutive interpretation are not inferred by this parse-only registration.

Evidence inspected:

- Input: `/tmp/opencode/tn_prof/tochnog_version_02-08-2026/test/incremental_driver/incremental_driver_syntax.dat`.
- Professional parser output: `/tmp/incdrv_validation_clean/professional/incremental_driver_syntax/incremental_driver_syntax.dbs`.

## Validation-only incremental_driver exit

An active `group_materi_plasti_hypo_strain_isa` metadata record is the only
trigger; there is no input flag. It makes the whole GNU `incremental_driver`
invocation validation-only. GNU validates the driver experiment records but does
not dispatch a material, execute ISA/Niemunis or hypo physics, perform a Newton
solve, or emit result rows. Instead it creates `incremental_driver_result.txt`
with the visible marker:

```
syntax-only / ISA not executed
```

The file also states that it contains no computed results and the program returns
rc=0. The trigger is global so a file cannot contain both an ISA metadata record
and a partial numerical result. Normal Mohr-Coulomb-only incremental-driver
execution is unchanged.

Any result requiring ISA/Niemunis constitutive evolution remains unsupported; no
substitute model is inferred or provided by this registration.
