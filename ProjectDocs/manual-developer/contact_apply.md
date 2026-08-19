# contact_apply / contact_plasti_friction / contact_target_*

## Implementación

- **Algoritmo**: `parallel_contact()` in `contact.cc` (preexistente). Los
  contactores son todos los nodos; los targets son geometrías
  (`CONTACT_GEOMETRY`) o elementos.
- **Keywords** (data_class CONTACT, `no_index` para los globales) en
  `database.cc`:
  - `contact_apply` (INTEGER 1, no_index) — gate: si `-no`, `parallel_contact`
    retorna antes del cálculo (manual 6.99).
  - `contact_plasti_friction` (DOUBLE 2, no_index) — `phi c`; `friction_limit
    = max(c + Fn*tan(phi), 0)` reemplaza el `mu*Fn` de `contact_friction`
    (manual 6.104). En el bloque de slip, el límite de fricción se calcula con
    la ley plástica cuando el record está presente.
  - `contact_target_element_group` (INTEGER variable, no_index) — filtro de
    targets por grupo (manual 6.105): en el loop de targets por elemento, si
    el grupo del elemento no está en la lista, se salta.
  - `contact_target_geometry`/`_switch` (INTEGER 2/1) — alias de
    `CONTACT_GEOMETRY`/`_SWITCH`: al inicio de `parallel_contact`, si los
    records target-geometry están activos y los base no, se copian.
- **Enums**: `CONTACT_APPLY`, `CONTACT_PLASTI_FRICTION`,
  `CONTACT_TARGET_ELEMENT_GROUP`, `CONTACT_TARGET_GEOMETRY`,
  `CONTACT_TARGET_GEOMETRY_SWITCH` en `tochnog.h`/`tochnog-mod.h`
  (alfabético).

## Validación

- `contact`: quad4 con target `geometry_line`, `contact_penalty_velocity`,
  `contact_plasti_friction`, `contact_apply -yes` — el modelo con la familia
  completa se parsea y corre de forma estable (RC=0, determinista). El gate
  `contact_apply -no` se verifica por inspección (retorna antes).
- Nota: el manual advierte que el algoritmo de contacto es experimental; la
  detección de penetración de la geometría depende del caso (en este modelo
  la penetración no se corrige — limitación del algoritmo de geometría).
