# group_interface

## Description

`group_interface` marks an element group as an **interface element**
group. Interface elements model joints or discontinuities between blocks
of material (e.g. between a pile and the soil). Their strains are the
displacement differences between the two opposite sides of the element,
not field gradients.

This is the first phase of the interface family (Carril A). Currently the
elastic interface law and the Fase 3 constitutive features are
implemented:

- `group_interface_materi_elasti_stiffness kn kt,first kt,second`:
  `stress_normal = kn * strain_normal`,
  `stress_shear = kt * 2 * strain_shear`.
- `group_interface_gap gap`: physical gap of the interface. Sign
  convention: **compression = positive normal strain**. The interface is
  CLOSED (full stiffness) when the accumulated normal strain
  `strain_normal > gap`, OPEN (residual stiffness) when
  `strain_normal <= gap`. A **negative** gap is a real gap: the interface
  stays open until compression exceeds |gap|. Without the record the
  interface is always closed.
- `group_interface_materi_residual_stiffness factor`: stiffness fraction
  used when the interface is open (default 0.01).
- `group_interface_materi_plasti_tension_direct tension_limit`: tensile
  limit on the TOTAL normal force `|kn*strain_normal|`; the interface opens
  in traction when the limit is exceeded (and it was still closed).
- `group_interface_materi_plasti_mohr_coul_direct phi c phi_flow`:
  cumulative Mohr-Coulomb friction (manual Professional 6.632, angles in
  RADIANS). The **presence of the record activates the law**: the trial is
  the ACCUMULATED tangential displacement times the stiffness
  (`kt*gamma_total`, the history keeps accumulating across plastic steps)
  and the resulting tangential force is clamped to the yield limit
  `max_fric = |c + Fn*tan(phi)|`, where `Fn` is the accumulated normal
  force (POSITIVE under compression). While sliding the tangential
  stiffness drops to 0 and the interface transmits the full limit force
  (spring pattern, the same as any other element: the assembled right
  hand side carries the FULL accumulated forces, so multi-step runs
  reach the static equilibrium in one step and report the total
  reactions - verified against the Professional `.dbs` of
  `mohr_coul_direct3`: node_rhside = weight*kn*eps_acc). With
  `phi=0, c=0` the limit is 0 → free sliding; without the record the
  interface is purely elastic. `phi_flow` is the dilatancy angle
  (non-associated flow): the INCREMENTAL plastic slip of each step opens
  the interface in both sliding directions
  (`strain_normal += -dgamma_inc*tan(phi_flow)`; the increment, not the
  total accumulated slip, so the dilatancy is not re-counted every step).
- `control_reset_dof -sigxx/-sigyy/-sigzz value` (normal stress reset):
  also initialises the accumulated normal strain of the interface
  elements whose normal aligns with the reset axis
  (`sigma_n := value`, i.e. `ELEMENT_INTERFACE_STRAIN_NORMAL :=
  value/kn`), matching the Professional behaviour (verified on
  `mohr_coul_direct4`: the `-sigyy -1` reset pre-stresses the horizontal
  interface to `sigma_n = -1` and the yield limit becomes
  `c + |Fn|*tan(phi) = 1.20271` instead of the unconfined cohesion c).
- `group_interface_materi_memory memory_type`: memory model of the
  constitutive law. `-updated_linear` (default) recomputes the interface
  normal/tangent from the current (deformed) configuration each step;
  `-total_linear` uses the time-0 reference geometry
  (`NODE_START_REFINED`) and keeps the original interface orientation.

In 2D the interface element is a quadrilateral with 4 nodes: nodes 0-1
form side 1 and nodes 2-3 form side 2.

## Uso

Place it in the data part, in the element group definition:

```
group_type 10  -materi
group_interface 10  -yes
group_interface_materi_elasti_stiffness 10  1000.0  0.0  0.0
```

## Parámetros

| Parameter | Meaning |
|-----------|---------|
| `10`      | Element group index. |
| `-yes`    | Activate interface behaviour for this group. |

## Parámetros (Fase 3)

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_interface_gap` | `gap` | Physical gap (negative = real gap, closes under compression). Closed when strain > gap, open (residual) when strain <= gap. Default without record: always closed. |
| `group_interface_materi_residual_stiffness` | `factor` | Stiffness fraction of an open interface (default 0.01). |
| `group_interface_materi_plasti_tension_direct` | `tension_limit` | Opens in traction when the total normal force `\|kn*strain_normal\|` exceeds the limit (and the interface was closed). |
| `group_interface_materi_plasti_mohr_coul_direct` | `phi c phi_flow` | phi = friction angle (rad), c = cohesion, phi_flow = dilatancy angle (rad). The record's presence activates the cumulative Mohr-Coulomb law on the TOTAL accumulated tangential force (`trial = kt*gamma_total` clamped to `\|c + Fn*tan(phi)\|`, Fn = accumulated normal force, positive under compression); phi=0,c=0 gives free sliding. |
| `group_interface_materi_memory` | `memory_type` | `-updated_linear` (default) or `-total_linear`. Memory model of the interface law; `-total_linear` fixes the normal/tangent to the time-0 geometry. |
| `control_reset_interface` / `control_reset_interface_strain` | `index geometry` | Reset the accumulated normal strain (and, for `_interface`, the tangential force histories) of the interface elements located in the geometry. Fires in the control step of its own index. |

## Convergence record (2026-09-03, interface Mohr-Coulomb direct)

Closed `mohr_coul_direct3` (shear test of an interface with
`group_interface_materi_plasti_mohr_coul_direct 0. 0. 0.785398`, target
node_rhside ±100). Three corrections, all verified against the
Professional binary `.dbs`:

1. **Full-force assembly (spring pattern)**: the interface now assembles
   the FULL accumulated forces into the element right-hand side
   (`stress_normal = kn*strain_normal_acc`, `stress_shear` = the clamped
   accumulated trial), not the last-step increment. The old incremental
   rhs reported the last-step increment as the reaction
   (`mohr_coul_direct3`: −50 instead of −100 per node) and made
   multi-step runs creep one increment per step under a constant load
   (`interface9`: 10 equal steps vs the Professional's single-step static
   equilibrium).
2. **Incremental dilatancy**: the plastic return of the accumulated-trial
   formulation gives the TOTAL plastic slip; the normal-strain update now
   uses only the plastic slip INCREMENT of the step, so the dilatancy
   ratchet is exact (`mohr_coul_direct3` sigma_n = −200 =
   kn*(−1e-6 − 1e-6*tan(pi/4)) for a total slip of 1e-6, not −250).
3. **Normal force sign in the yield limit**: `Fn` is positive under
   compression (`max_fric = |c + Fn*tan(phi)|`; the Professional plateau
   of `mohr_coul_direct4` is c + |Fn|*tan(phi) = 1.20271 with c=1,
   phi=0.2, sigma_n = −1), and a `control_reset_dof -sigxx/-sigyy/-sigzz`
   reset pre-stresses the aligned interface normal stress
   (`sigma_n := value`, `ELEMENT_INTERFACE_STRAIN_NORMAL := value/kn`).
4. **`control_reset_interface` / `_interface_strain` fixed**: the reset
   was nested inside the `control_reset_dof` gate (it never fired when
   only the interface reset record was present) and wrote a single value
   with a leftover length; it now scans the control range with
   `db_active_index` and zeroes ALL per-integration-point history slots
   in both versions (`interface7` re-verified against its −1e10 target
   with the accumulated records).

Result: `mohr_coul_direct3` rc=0 (σn −199.99997, node_rhside
±99.99998 — identical to the Professional `.dbs`). Corpus 140 → 141
PASS. `mohr_coul_direct4` remains RUNFAIL (GNU 0.251052 vs target
0.10244): its Professional trajectory needs the `-sigyy` reset to act as
a PERSISTENT pre-stress on the blocks and the interface (the stress-dof
state is held for 100 s), which the GNU's kinematic interface normal
stress cannot sustain — see the developer manual.

## Related

- `group_interface_materi_elasti_stiffness index kn kt,first kt,second` —
  elastic interface stiffness (normal kn, tangential kt).

## Estado de implementación

- **Implementado (Fase 1)**: elastic interface element (2D quadrilateral),
  `group_interface`, `group_interface_materi_elasti_stiffness`.
- **Implementado (Fase 3)**: `group_interface_gap`,
  `group_interface_materi_residual_stiffness`,
  `group_interface_materi_plasti_tension_direct`,
  `group_interface_materi_plasti_mohr_coul_direct` (Mohr-Coulomb acumulativo
  con history `element_interface_force_tang`; `phi_flow` = dilatancia),
  `group_interface_materi_memory` (`-updated_linear`/`-total_linear`).
  Validados con la familia `iface_mc` (13º test de `build_safe.sh`,
  10 runs): fricción alta sostiene la carga tangencial, fricción nula
  desliza libre, tracción abre la interfaz, gap cierra bajo compresión,
  dilatancia (`phi_flow`) abre la interfaz bajo deslizamiento plástico, y
  el clamp del MC reproduce `max_fric = c + kn*strain*tan(phi)` de forma
  numérica exacta.
- **Implementado (Fase 2)**: `control_mesh_convert` — conversion
  automatica de `-bar2` a `-quad4` para interfaces: crea los 2 nodos del
  lado opuesto de la interfaz y reconecta los elementos vecinos del otro
  lado (con `control_mesh_convert_element_group` para los grupos a un
  lado).
- **Pendiente**: `group_interface_materi_memory`, conductivity/groundflow,
  y el post-proceso de interfaz. See `ProjectDocs/DESIGN-INTERFACES.md`.
