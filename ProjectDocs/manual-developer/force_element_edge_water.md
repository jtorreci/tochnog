# force_element_edge_water (force_edge_water) — layout dual

## El conflicto de layouts (documentado en SEGUIMIENTO, sección force_edge)

El nombre `force_element_edge_water` (Professional `force_edge_water`,
manual 6.489) tiene DOS layouts incompatibles:

1. **Legacy GNU**: `force_element_edge_water <index> <rho> <g> <dirx>
   [<diry> <dirz>]` — la densidad, la gravedad y la dirección de la
   presión se dan explícitas (el force5.dat original del GNU-2014 usa
   `1. 4.9 -1. 0.`).
2. **Professional**: `force_edge_water <index> -yes` — un SWITCH puro: la
   presión hidrostática se calcula automáticamente como `rho*g*delta_z`
   con `rho = groundflow_density`, `g = la componente vertical de
   force_gravity` (con signo; negativa hacia abajo) y `delta_z` la
   profundidad bajo `groundflow_phreatic_level`; actúa normal al borde
   del elemento en dirección INTERIOR (= -normal saliente).

Ambos layouts conviven en el corpus (force5 y excavate1 usan el switch;
ningún test del corpus usa el layout legacy, pero el GNU-2014 original de
force5 sí).

## Implementación (2026-09-07)

### Parse del switch (input.cc)

`FORCE_ELEMENT_EDGE_WATER` es un record DOUBLE_PRECISION y el parser
genérico rechaza el token `-yes`. En la rama DOUBLE del lector de
valores de `input()` se acepta para ESTE record un token `-yes`/`-no`
(con el guion) y se almacena como el enum negativo (como en los records
INTEGER): `dval = -(double)db_number("yes")`. Cualquier otra palabra
sigue siendo un error de parse.

El record se registra ahora con `fixed_length = 0` (longitud variable):
1 valor para el switch, `ndim+2` para el layout legacy (el máx. de
valores se mantiene en `data_length`).

### Consumidor (area.cc)

En la rama de carga de `FORCE_ELEMENT_EDGE_WATER` (area.cc), el layout se
discrimina por `ldum==1 && (long)values[0]==-YES`:

- **Auto (-yes)**: `rho` del record `groundflow_density` (default 1.0),
  `grav[] = force_gravity_calculate()`, y el vector de carga =
  `normal` saliente del lado (ya acumulado y normalizado por el
  mecanismo de area()) multiplicado por la presión CON SIGNO
  `rho * grav[ndim-1] * delta_z` (grav negativa → carga hacia el
  interior del elemento; convención consistente con la presión estática
  de groundfl.cc). Solo cuando `delta_z > 0` (nodo bajo el nivel
  freático) y `ndim > 1`.
- **Legacy**: sin cambios (`values[0]*values[1]*delta_z` a lo largo de
  `values[2..ndim]` normalizado).

### Bug estructural corregido: restricciones _element/_element_group/
_element_side muertas

El bloque de restricciones de elemento de las familias
`force_element_edge*` (companions `_element`, `_element_group`,
`_element_side`) estaba anidado DENTRO de la rama condif
(`if type==CONDIF_RADIATION || CONDIF_CONVECTION || ...`) del bucle de
nodos de area(): para los tipos force la condición exterior es falsa y
el bloque NUNCA se ejecutaba (código muerto). Consecuencia: las
restricciones de elemento/grupo/lado de la familia force_edge no se
aplicaban en runtime — el registro `force_edge_water_element_group` de
force5 no filtraba y las cargas de agua de los dos elementos que
comparten la interfaz se CANCELABAN (cada elemento tira de su nodo
común hacia su interior: −F de la pared y +F del suelo → neta 0 →
sigxx = 0 en lugar de −9.8). El bloque se movió al inicio del bucle por
nodo, ejecutándose para todos los masters force (`force_edge_is_master`),
antes del dispatch por tipo. Las restricciones `_node`/`_element_node`
se siguen aplicando dentro de cada rama de carga (como antes).

### Keyword nuevo: groundflow_phreatic_level_static

`groundflow_phreatic_level_static` (manual 6.57x, la versión single del
`_multiple_static`): registro INTEGER no_index registrado en database.cc/
tochnog.h para que excavate1 (mecánica sin dof de presión) parsee. La
semántica de flujo (fijar la presión total = estática en los nodos del
nivel) pertenece a la familia groundflow (groundfl.cc) y queda pendiente;
A/B contra el Professional (2026-09-07): quitando el record Y el
`group_type -groundflow` de excavate1 el Pro sigue dando rc=0 (el target
solo depende de la carga de agua + la gravedad + el borrado de
elementos), así que el registro parse-only es suficiente para el test.

## Validación

- force5 (corpus, layout Professional -yes + `_element_group 1`): rc=0 —
  post_point sigxx = −9.8071 (target −9.8 ± 0.1), idéntico al
  Professional (−9.807155). A/B: sin el fix de las restricciones el
  resultado era sigxx = 0 (cargas canceladas).
- force6 (corpus): rc=0 (mejorado por el fix de restricciones — el A/B
  con el binario pristino da rc=1).
- Layout legacy (A/B local, force5 con `1. 9.8 -1. 0.` y sin
  restricción): la carga alcanza el equilibrio (el factor 2 por la carga
  de ambos elementos reproduce la semántica original del GNU-2014 con
  g=4.9).
- excavate1: el parse de `force_edge_water ... -yes` + `_time` + la
  keyword `_static` funcionan; el test queda RUNFAIL por blockers de la
  familia groundflow (ver diagnóstico en el reporte del sprint).

## Pendiente

- `force_edge_water -no` se parsea pero no tiene semántica de apagado
  (nunca se usa en el corpus).
- La semántica de flujo de `groundflow_phreatic_level_static` (familia
  groundflow).
