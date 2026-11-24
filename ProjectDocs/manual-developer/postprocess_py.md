# postprocess.py

## Implementación

- **Ubicación**: `tools/postprocess.py` (Python 3, `argparse`, `pandas`,
  `numpy` opcional, `matplotlib` opcional).
- **Carga**: `load_sqlite()` lee todas las tablas del `.sqlite` a
  DataFrames y `meta` a un dict. `pivot_primary()` convierte el formato
  largo `primary_data(node,dof,t,value)` a ancho `(node,t)` x columnas
  de dof. `wide_with_derived()` fusiona además la tabla `derived` por la
  clave `(node,t)`.
- **Subcomandos**: `info`, `stats`, `plot`, `line`, `user` (dispatch por
  diccionario en `main()`).
- **Línea de perfil**: se proyecta cada nodo sobre el segmento `p0->p1`,
  se filtran los nodos con distancia perpendicular `<= tol` y parámetro
  `s` en `[0,1]`; la salida ordena por `s`.

## Dependencias y requisitos de build

- Python 3.12 con `pandas` (verificado). `numpy`/`matplotlib` solo para
  `plot` y `line --plot`.
- Requiere que el `.sqlite` se haya generado con `SQLITE_USE=1`
  (compilación con soporte SQLite; ver `control_print_tabular` developer).

## Notas

- Las columnas `sigX_00..sigX_22` de un dof matriz: tochnog registra cada
  componente de stress como un dof matriz que alias al mismo tensor, por
  lo que `sigxx_11 == sigyy_00` etc. son el mismo valor físico.
- La tabla `user_data` usa formato largo `(node, variable, t, value)`;
  `cmd_user` escribe con `INSERT OR REPLACE`.

## Pendiente

- La interpolación sobre la línea usa solo los nodos del almacén (sin
  refinamiento); para perfiles suaves haría falta interpolación espacial
  (p.ej. scipy `griddata`).
- `stats` no filtra por rango de tiempo (agrega sobre todos los `t`).
- El `line` ignora `z`; para geometrías 3D habría que generalizarlo.
