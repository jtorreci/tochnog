# group_materi_elasti_k0 (+ hook de control_materi_elasti_k0)

## Implementación

- **Enum nuevo**: `GROUP_MATERI_ELASTI_K0` en tochnog.h y tochnog-mod.h
  (en sync, orden alfabético entre `ELASTI_COMPRESSIBILITY` y
  `ELASTI_LADE`).
- **Registro** (database.cc): DOUBLE_PRECISION, length 1, data_class
  MATERI, data_required GROUP_TYPE: `K0`. El control
  `CONTROL_MATERI_ELASTI_K0` (INTEGER, length 1, CONTROL) YA existía
  como parcial de Sprint 9 (sin hook); se conecta ahora.
- **Hook en `set_stress()`** (stress.cc), ANTES del bloque young:
  ```c
  k0_active = 0;
  if ( db_active_index( CONTROL_MATERI_ELASTI_K0, 0, VERSION_NORMAL ) ) {
    idum[0] = 0;
    db( ICONTROL, 0, idum, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
    db( CONTROL_MATERI_ELASTI_K0, idum[0], &k0_control_swit, ddum, ldum,
      VERSION_NORMAL, GET_IF_EXISTS );
    if ( k0_control_swit==-YES &&
         get_group_data( GROUP_MATERI_ELASTI_K0, gr, element,
           new_unknowns, &k0_elasti, ldum, GET_IF_EXISTS ) ) {
      if ( k0_elasti>0.95 ) k0_elasti = 0.95;
      k0_active = 1;
    }
  }
  ```
  Y en los bloques young y young_power, tras leer el poisson:
  ```c
  if ( get_group_data( GROUP_MATERI_ELASTI_POISSON, gr, element,
      new_unknowns, &poisson, ldum, GET_IF_EXISTS ) && k0_active )
    poisson = k0_elasti/(1.+k0_elasti);
  ```
- **Guarda estricta**: el override aplica SOLO con
  young/young_power Y poisson presentes (el return de get_group_data del
  poisson se captura) Y control -yes Y record k0 existente — la
  combinación exacta del manual 6.650.
- **check.cc**: `GROUP_MATERI_ELASTI_K0` requiere `materi_stress` +
  `materi_velocity`.

## Física

Manual 6.650: con K0 + control -yes, `nu = K0/(1+K0)` (derivado de
`K0 = nu/(1-nu)`), K0 > 0.95 truncado a 0.95. En un oedometro
lateralmente confinado el ratio sigma_lateral/sigma_axial =
`nu/(1-nu) = K0` — "K0 stresses".

## Validación

- `mk0.dat` (control -yes, K0=0.5): oedometro `E=1000`, `eps_zz=-0.001`:
  nu = 1/3 -> sigma_xx = -0.75, sigma_yy = -1.5 EXACTOS
  (sigma_xx/sigma_yy = 0.5 = K0).
- `mk0_off.dat` (control -no): nu base 0.3 -> sigma_xx = -0.5769,
  sigma_yy = -1.3462 — A/B discrimina.

## Gotchas

- El control se lee con el patrón de `control_materi_gate_off`
  (general.cc): ICONTROL + índice del control (idu m[0] inicializado a 0
  porque GET_IF_EXISTS no escribe si el record falta).
- El control es INTEGER (swit -yes/-no): el buffer del valor debe ser
  long int, NO double (error de compilación corregido en el desarrollo).
- La combinación con `group_materi_elasti_hardsoil` del manual NO está
  implementada (hardsoil no existe en el GNU) — pendiente documentado.
- young_polynomial NO recibe el override (el manual solo lista
  young/young_power) — decisión deliberada.
