# group_materi_plasti_cap2 (alias Professional del cap legacy)

## Implementación

- **Keyword**: `group_materi_plasti_cap2` resuelta como **alias en
  `db_number()`** (database.cc) al enum legacy `GROUP_MATERI_PLASTI_CAP`:

  ```c
  else if ( !strcmp( str, "group_materi_plasti_cap2" ) )
    return GROUP_MATERI_PLASTI_CAP;
  ```

  Mismo mecanismo que `group_materi_plasti_druck_prag` -> `DRUCKPRAG`
  (Sprint 10 lote 2): el `name[]` del enum se queda con el nombre GNU
  (`group_materi_plasti_cap`) y el nombre Professional se resuelve en
  db_number. **plasti.cc NO se toca** — no hay duplicación de física.
- **Datos** (database.cc, `GROUP_MATERI_PLASTI_CAP`): DOUBLE_PRECISION,
  length DATA_ITEM_SIZE (variable), data_class MATERI, data_required
  GROUP_TYPE. Layout: `c phi alpha R` + tabla `epsilonp_v pb` (>= 2 pares).
- **Física**: bloque `group_materi_plasti_cap` en `plasti_rule()`
  (plasti.cc:186-231): invariantes `p = -sigm`, `t = sqrt(3)*sig_eq`;
  `pa = (pb - R*c)/(1 + R*tan(phi))`; yield
  `f = sqrt((p-pa)^2 + (R*t/tmp)^2) - R*(c + pa*tan(phi))` con
  `tmp = 1 + alpha - alpha/cos(phi)`. En el camino isotrópico (t=0) se
  reduce a `f = p - pb` — condición elástica `p < pb`.
- **check.cc**: `GROUP_MATERI_PLASTI_CAP` requiere `materi_stress` +
  `materi_strain_plasti` (sin cambios — usa enums, no strings).

## Validación

- `mcap2.dat` / `mcap_legacy.dat` (validation-suite/test-2014):
  oedometro confinado, `E = 1000`, `nu = 0.3`, `eps_zz = -0.001`,
  cap con `c = 1e6` (f = p - pb = 0.83 - 100 < 0 -> elástico).
  `sigma_xx = -0.5769` EXACTO con AMBAS keywords (cap2 == cap).
- Verificado que ningún test legacy (test-2014, sfnet) usa el nombre
  `group_materi_plasti_cap`, por lo que el alias no rompe nada.

## Gotchas

- El valor de `pb` de la tabla se obtiene con `table_xy()` sobre
  `epp_vol` (la suma de las deformaciones plásticas diagonales en valor
  absoluto); con `c` enorme el cap queda lejos y la respuesta es elástica.
- `group_materi_plasti_cap1` sigue PENDIENTE (familia distinta).
