# group_materi_plasti_cap2

## Description

`group_materi_plasti_cap2` is the **Professional name** (manual
Professional 6.692) of the cap plasticity model of the GNU fork, registered
under the legacy keyword `group_materi_plasti_cap`. Both names resolve to
the SAME material law and data layout — `c phi alpha R` plus the
`epsilonp_v pb` table. The alias is resolved inside `db_number()`
(database.cc), so the physics in plasti.cc is shared and never duplicated.

The cap model is an ellipsoidal cap plasticity (Drucker-Prager-type yield
plus a hardening cap) used for soils and granular materials: cohesion `c`,
friction `phi`, the cap-shape parameter `alpha` and the cap radius `R`,
with a hardening table that maps the plastic volumetric strain
`epsilonp_v` to the hardening pressure `pb`.

## Uso

```
group_type 0  -materi
group_materi_plasti_cap2 0
                        1.e6    ( c )
                        0.5     ( phi, rad )
                        0.5     ( alpha )
                        0.1     ( R )
                        0.      ( epsilonp_v, pb table: >= 2 pairs )
                        100.
                        1.
                        500.
```

The table needs **at least two** `epsilonp_v pb` pairs. The record length is
variable (DATA_ITEM_SIZE).

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_plasti_cap2` | `c phi alpha R e1 p1 e2 p2 ...` | Cohesion, friction angle (rad), cap shape `alpha`, cap radius `R`, and the `epsilonp_v` vs `pb` hardening table. |

On the isotropic compression path the yield function reduces to
`f = p - pb`, so the element is elastic while the pressure `p` stays below
the current `pb` of the table.

## Validation

- `mcap2.dat` / `mcap_legacy.dat` (validation-suite/test-2014): confined
  compression (oedometer) with `c = 1e6` (huge cohesion -> `f = p - pb < 0`,
  elastic regime), `E = 1000`, `nu = 0.3`, `eps_zz = -0.001`. Both the
  Professional keyword (`mcap2`) and the legacy keyword (`mcap_legacy`)
  give **sigma_xx = -0.5769 EXACT** — the same value: `cap2 == cap`.

## Notas

- Requires `materi_stress` and `materi_strain_plasti` in the initialization
  part (the cap yield uses the plastic strain `epsilonp_v`).
- `group_materi_plasti_cap1` (a different cap family) remains PENDING.
