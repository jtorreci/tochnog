# materi_elasti_young_power_apply

## Implementación

- **Enum**: `CONTROL_MATERI_ELASTI_YOUNG_POWER_APPLY` en tochnog.h y
  tochnog-mod.h (en sync; entre `CONTROL_MATERI_ELASTI_K0` y
  `CONTROL_MATERI_FAILURE_APPLY`).
- **Registro** (database.cc): INTEGER, length 1, data_class CONTROL.
- **Alias** (db_number, database.cc): el nombre del manual 6.801
  `materi_elasti_young_power_apply` (sin prefijo control_) resuelve a
  `CONTROL_MATERI_ELASTI_YOUNG_POWER_APPLY` (patrón de los demás gates
  control_materi).
- **Consumo** (stress.cc, dentro del bloque
  `GROUP_MATERI_ELASTI_YOUNG_POWER`): helper estándar
  `control_materi_gate_off( CONTROL_MATERI_ELASTI_YOUNG_POWER_APPLY )`
  (general.cc, ICONTROL-indexado, default activo):
  ```c
  if ( control_materi_gate_off( CONTROL_MATERI_ELASTI_YOUNG_POWER_APPLY ) )
    young = young0;               // E = E0 constante
  else {
    young = young0 + young1 * scalar_power( p/p1, alpha );
    if ( young<young2 ) young = young2;
    if ( young>young3 ) young = young3;
  }
  ```

## Física

Manual 6.801: "If switch is set to -no, any nonlinearity in young
dependent on a power law will be ignored; simply the constant young as
encountered in the group_materi_elasti_young_power records will be
applied at all times." La constante del RECORD es `E0` (primer
parámetro) — el registro young_power es autosuficiente (por eso el
bloque limpia C/Cmem: la ley ES el módulo, no se suma al young).

## Validación

- `myoung6_apply.dat`: oedometro E0=1000 E1=500 E2=800 E3=3000 p1=1
  alpha=1 con el gate -no -> E = E0 = 1000 constante: sigxx -0.5769 /
  sigyy -1.3462 EXACTO (= la base lineal; A/B contra myoung6 cuyo
  punto fijo da -0.713886/-1.66573). El gate discrimina.

## Gotchas

- El gate aplica SOLO al young_power: el `young_polynomial` (dependiente
  de deformación) y el `young_strainstress` no se ven afectados (así lo
  dice el manual: "young dependent on a power law").
- El check de la ley (p1>0) se ejecuta SIEMPRE (incluso con el gate
  -no): un p1<=0 sigue siendo un error de input.
