# group_materi_plasti_bounda / group_materi_plasti_bounda_factor

Reducción de la fricción de materiales granulares en paredes (manual
Professional 6.231/6.232). Alias Professional de los records GNU
`group_materi_plasti_boundary`/`group_materi_plasti_boundary_factor`
(resolución en db_number, database.cc).

## Semántica (6.231/6.232)

- `group_materi_plasti_bounda <index> <i0> <i1> ...`: los valores son los
  ÍNDICES de los records `bounda_dof` que definen las paredes: "an
  element is on a wall when at least one of the velocities
  (displacements) of the element is prescribed (via bounda_dof)".
- `group_materi_plasti_bounda_factor <index> <factor> [<factor_i1> ...]`:
  factor de reducción (default 2./3. cuando el record está pero el factor
  no; con un solo valor se aplica a todas las paredes). Un factor 0 =
  pared lisa (sin fricción).
- Se reduce phi, phi_flow y c (mohr_coul / mohr_coul_direct / druck_prag
  / hardsoil), M (camclay) o los incrementos desviadores
  (hypo_*).

## Implementación

### Detección de pared (group.cc)

`group_materi_plasti_boundary_evaluate(nodes, nnol, element_group,
plasti_on_boundary)` (llamada por elem.cc una vez por elemento antes del
bucle de puntos de integración) mantiene la semántica legacy (los valores
del record como grupos de elemento: el elemento está en la pared si su
propio `element_group` está en la lista) Y añade la semántica del
Professional (6.231): cuando un valor coincide con un record `bounda_dof`
ACTIVO, el elemento está en la pared si uno de sus nodos está bounded
(`NODE_BOUNDED`) en las partes de velocidad o desplazamiento
(`vel_indx..vel_indx+ndim`, `dis_indx..dis_indx+ndim`).

### Consumo del factor

- Leyes incrementales (plasti.cc, `plasti_rule`): `phi *= factor;
  phi_flow *= factor;` cuando `plasti_on_boundary` (la cohesión en los
  bloques clásicos de plasti.cc sigue su propio manejo — revisar por ley).
- Leyes directas (stress.cc): `materi_direct_cutoff` (modo plano, con
  `_normal`) multiplica `phi` y `c` por el factor del grupo;
  `materi_direct_full_mc` (modo espectral) usa los records `_wall` cuando
  `plasti_on_boundary` (los parámetros alternativos explícitos de la
  familia _wall); el factor de plasti_bounda no se aplica en el modo
  espectral (la semántica _wall es la del Professional para ese modo).

## Validación (2026-09-07)

A/B discriminante (compresión + cortante, modo plano con normal (0,1),
`-total_linear`): sin el record la fricción aguanta (sigxy = 2500 =
elástico: max_fric = 4228 > 2500); con
`group_materi_plasti_bounda 0 10` (10 = índice del bounda_dof que fija
los nodos inferiores) y factor 0, sigxy = 0 (pared lisa). El mecanismo
completo (detección por bounda_dof + factor en runtime) está activo.

## Estado de los tests del corpus

- mohr_coul_direct8: parsea, la reducción de pared se aplica, pero queda
  RUNFAIL por la cinemática `-updated` del default GNU con el incremento
  de cortante 100 % (rotación polar; la versión `-total_linear` pasa).
  Detalle en group_materi_plasti_mohr_coul_direct.md.
- La reducción en las leyes incrementales clásicas (mohr_coul, druck_prag,
  hardsoil, camclay, hypo) usa el mismo flag `plasti_on_boundary`
  calculado por group.cc.
