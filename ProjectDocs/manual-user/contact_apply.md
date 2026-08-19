# contact_apply / contact_plasti_friction / contact_target_*

## Description

The `contact_*` records configure the **contact algorithm**: contacting
nodes (the contacter) are kept on the outward-normal side of a target
surface (a geometry or target elements). The contact penalty puts a spring
between the contacter and the target when penetration occurs.

The GNU already implements `contact_geometry` (target surface),
`contact_penalty_velocity/pressure/temperature`, `contact_friction`,
`contact_stick`, `contact_relaxation` and `contact_heatgeneration`. The
following records complete the family (manual 6.99-6.107):

- **`contact_apply index switch`** — enables (`-yes`) or disables (`-no`)
  the contact algorithm for all timesteps. Default: enabled when any
  contact data is present.
- **`contact_plasti_friction phi c`** — plastic (Mohr-Coulomb) friction on
  the contact surface: `max_fric = max(c + Fn*tan(phi), 0)`. Takes
  precedence over the simple `contact_friction mu*Fn`. The friction opposes
  the tangential sliding and its energy is partly converted to heat
  (`contact_heatgeneration`).
- **`contact_target_element_group g_0 g_1 ...`** — restricts the target
  elements to the listed element groups.
- **`contact_target_geometry` / `contact_target_geometry_switch`** — alias
  of `contact_geometry`/`contact_geometry_switch` (the manual names).

## Uso

Place it in the data part:

```
start_define
  contact_face geometry_line 1
end_define
contact_face 0. 1.  1. 1.  0.01
contact_geometry 0  -geometry_line 1
contact_penalty_velocity 100000.
contact_plasti_friction 0.785398 0.5
contact_apply -yes
```

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `contact_apply` | `-yes`/`-no` | Enable/disable the contact algorithm. |
| `contact_plasti_friction` | `phi c` | Mohr-Coulomb friction (`max(c + Fn*tan(phi), 0)`). |
| `contact_target_element_group` | groups | Target element groups. |
| `contact_target_geometry` | geometry | Target geometry (alias of `contact_geometry`). |
| `contact_penalty_velocity` | penalty | Normal penalty for penetration. |
| `contact_friction` | mu | Simple Coulomb friction coefficient. |

## Related

- `contact_geometry` / `contact_geometry_switch` — the target surface.
- `contact_penalty_pressure/temperature` — pressure/temperature penalties.
- `contact_stick`, `contact_relaxation`, `contact_heatgeneration`.

## Estado de implementación

- **Implementado**: `contact_apply` (gate), `contact_plasti_friction`
  (Mohr-Coulomb), `contact_target_element_group` (filtro de grupos),
  `contact_target_geometry`/`_switch` (alias). El algoritmo base
  (`parallel_contact` en contact.cc) ya existía. Validado con `contact`
  (el modelo con la familia se parsea y corre de forma estable). Nota: el
  manual advierte que el algoritmo de contacto es experimental; la
  detección de penetración de la geometría depende del caso.
