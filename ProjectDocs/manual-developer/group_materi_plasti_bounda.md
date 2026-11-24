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

### Factor consumption (C114, 2026-11-06)

- Incremental laws (`plasti.cc`, `plasti_rule`): `phi`, `phi_flow` and `c`
  times the factor at the five Mohr-Coulomb-type sites (mohr_coul, its
  hardening/softening variant, druckprag, matsuokanakai, legacy mohrcoul).
  The factor is `group_materi_plasti_boundary_factor`, or
  `group_materi_plasti_mpc_factor` for a wall found through mpc records
  (`plasti_on_boundary` 2); default `BOUNDARY_REDUCTION_FACTOR` 2/3.
- Direct laws (`stress.cc`): `materi_direct_full_mc` and
  `materi_direct_cutoff` use the `_wall` record when present (no factor on
  top), else phi, c and phi_flow times `direct_wall_factor()` (same factor
  selection and default 2/3).
- Before C114 the cohesion was not reduced in `plasti.cc`, the spectral
  direct mode ignored the factor and the plane mode had default 1 and did
  not reduce phi_flow. Professional 02-08-2026 probes (t3/c2) fixed the
  rule: an explicit group with phi and c times 2/3 reproduces the wall
  result exactly (Mohr-Coulomb, hardening/softening, Drucker-Prager, both
  direct modes).

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
- Wall code exists only for the Mohr-Coulomb-type laws of `plasti.cc` and
  the direct laws of `stress.cc` (C114). hardsoil, camclay (M), hypo
  (deviatoric increments) and the incremental tension (sigy) listed by the
  manual have no wall reduction in the GNU (BACKLOG C115).
- The detection is broader than in the Professional (BACKLOG C115): any
  velocity-bounded node marks the element, whatever bounda_dof index
  bounded it, and a listed value also matches element groups.
- `group_materi_plasti_bounda_factor` holds one value (data_length 1); the
  per-index factors of the manual are not supported (BACKLOG C115).
