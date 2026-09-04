# group_interface

## Implementación

- **Elemento**: `interface_element()` in `interface.cc` (new file, added
  to the makefile as `interface.$(OBJ)`). Invoked from `elem()` in
  `elem.cc`, in the structural-elements branch, when
  `db_active_index( GROUP_INTERFACE, element_group )` is true.
- **Keywords** (data_class GROUP_TYPE) registered in `database.cc`:
  - `group_interface` (INTEGER, length 1, required GROUP_TYPE).
  - `group_interface_materi_elasti_stiffness` (DOUBLE, length 3, required
    GROUP_INTERFACE).
  - `group_interface_materi_memory` (INTEGER, length 1, required
    GROUP_INTERFACE).
- **New enums**: `GROUP_INTERFACE`, `GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS`
  in `tochnog.h` / `tochnog-mod.h` (kept in sync), placed between
  `GROUP_INTEGRATION_POINTS` and `GROUP_MATERI_DAMAGE_MAZARS`.

## Modelo físico

- Interface strains = displacement differences between the two opposite
  sides (not gradients). In 2D the element is a quadrilateral with 4
  nodes; nodes {0,1} form side 1 and nodes {2,3} side 2.
- Elastic law:
  - `stress_normal = kn * strain_normal`
  - `stress_shear = kt * 2 * strain_shear`
  where strain = displacement difference / thickness (thickness absorbed
  in kn/kt, so thickness = 1).
- The displacement dof is `veli_indx` (velocity_integrated) or `dis_indx`
  (displacement), selected via `materi_velocity_integrated`.

## Ley constitutiva (Fase 3)

- **Sign convention (RF-5)**: compression = `strain_normal` POSITIVE,
  traction NEGATIVE. This is the convention used by gap, tension_direct
  and `max_fric` (`kn*strain_normal` grows with compression). The old doc
  said "compression negative" — wrong, corrected 2026-08-14.
  NOTE (2026-09-03, verified against the Professional .dbs of
  mohr_coul_direct3): the orientation flip of the 2D normal (side-2 nodes
  numbered lower than side-1) makes the RECORDED accumulated strain
  NEGATIVE under compression (direct3 strain −1e-6, sigma_n −100 at
  t=1) — the internal sign of `strain_normal` depends on the element node
  numbering; the physical branch decisions (gap/tension/max_fric) only
  compare magnitudes/signs consistently within one element. The yield
  limit uses `max_fric = |c - kn*strain_normal*tan(phi)|` so that the
  normal FORCE is positive under compression regardless of the numbering.
- **Gap** (`group_interface_gap`, CONVERGENCE 2026-09-04): the
  interface is CLOSED when the accumulated normal strain `<= gap` and
  OPEN when `strain_normal > gap` (manual Professional 6.625: "Only
  when the sides displacements are such that the normal strain becomes
  lower then the specified gap value the interface will be closed and
  start to generate stresses"). An opened interface "does not have
  stresses" (6.628): the stress record/rhs are exactly 0 and only the
  residual stiffness (`group_interface_materi_residual_stiffness`,
  default 0.01) stays in the matrix for regularization. A physical gap
  is a NEGATIVE gap value; the default without the record is +1e20
  (always closed — "if you want to allow always tension stresses set
  gap to, by example, 1.e20"). The accumulated strain is stored in
  `ELEMENT_INTERFACE_STRAIN_NORMAL`. NOTE (2026-09-04): the pre-fix code
  inverted the condition (open when `strain <= gap`) with a −1e20
  default — correct for tests without a gap record, but it left
  interfaces with an explicit positive gap (patch1: 1.e20) ALWAYS open
  and the closed phase of a physical gap carried only the residual
  stress.
- **Closed-phase stress accumulation** (`ELEMENT_INTERFACE_FORCE_NORM`,
  new internal history per integration point, registered in
  `database.cc`, allocated in `top.cc` for every interface group):
  `stress,normal` is NOT `kn * strain,normal_total`; it is kn times the
  normal strain of the steps that END CLOSED (the free travel of an
  open gap never builds stress). Per step and per IP:
  `force_norm += kn * (delta_mechanical_strain - thermal_increment)`
  when the step ends closed, `force_norm = 0` when it ends open (an
  opened interface does not have stresses; a later re-closure rebuilds
  the stress from the penetration of the closing step). Verified
  step-by-step against the Professional per-step prints of `interface2`
  (gap 0.1, 200 steps of −1e-3: open with stress 0 during the free
  travel, stress −1 at the closing step 100 and −101 = kn·(−101·1e-3)
  at step 200 — NOT kn·(−0.2)). Without a gap record the interface is
  always closed and `force_norm` equals `kn*strain,normal` (all the
  non-gap corpus tests are bit-identical up to FP round-off).
- **Tension limit** (`group_interface_materi_plasti_tension_direct`): the
  interface opens in TRACTION when the accumulated TOTAL normal force
  `|Fn_total| = |kn*strain_normal| > tension_limit` (requires
  `strain_normal < 0`) and only if it was still closed. On opening the
  state becomes OPEN (stress 0, residual stiffness), the same handling
  as the gap-open state.
- **Residual stiffness** (`group_interface_materi_residual_stiffness`):
  fraction of the original stiffness used as matrix regularization when
  the interface is open (default 0.01; no corpus test sets it).
- **Cumulative Mohr-Coulomb** (`group_interface_materi_plasti_mohr_coul_direct
  phi c phi_flow`, manual Professional 6.632 — angles in RADIANS): active
  by the PRESENCE of the record (D2). With phi=0,c=0 the limit is 0 →
  free sliding; without the record the interface stays purely elastic
  (Fase 1). The limit applies to the TOTAL tangential force, stored in
  the history `ELEMENT_INTERFACE_FORCE_TANG` (spring.cc pattern, GET
  VERSION_NORMAL with GET_IF_EXISTS, PUT VERSION_NEW):
  - `trial = f_t_old + kt1*du_tang` — the stored history is the ELASTIC
    trial `kt*gamma_total` (gamma_total = sum(du_tang) over ALL steps;
    it keeps accumulating even while plastic, so the force can reach the
    cohesion after a large slip — interface_patch -10)
  - `max_fric = |c + Fn*tan(phi)|` with `Fn = -kn*strain_normal`
    (POSITIVE under compression — verified against the Professional
    plateau of mohr_coul_direct4: c + |Fn|*tan(phi) = 1.20271 with
    c=1, phi=0.2, sigma_n = -1; the manual text "Fn negative under
    compression" is a sign slip of the manual)
  - clamp: `f_t = clamp(trial, ±max_fric)`; if plastified →
    `stiff_tang = 0` (consistent tangent, D3, avoids Newton oscillation
    at the elastic/plastic boundary)
  - CONVERGENCE (2026-09-03): the assembled rhs carries the FULL
    accumulated forces (spring.cc pattern), NOT the step increment:
    `stress_normal = kn*strain_normal_acc` and `stress_shear = f_t`
    (elastic: `kt*gamma_total`). With an incremental rhs the node
    reactions of a multi-step run only showed the last increment
    (mohr_coul_direct3: −50 vs the Professional −100 per node =
    w*kn*eps_acc) and a constant load made the interface creep one
    increment per step (interface9: 10 equal steps vs the Professional
    single-step static equilibrium).
- **phi_flow = dilatancy (RF-4, non-associated flow)**: if plastified and
  `phi_flow > 0`, the INCREMENTAL plastic slip of the step opens the
  interface in both sliding directions:
  `strain_normal += -dgamma_inc * tan(phi_flow)` where `dgamma_inc =
  dgamma - gamma_pl_old` (the accumulated-trial return multiplier
  dgamma is the TOTAL plastic slip; the slip already accumulated in the
  past, gamma_pl_old = (|trial_old| - |clamp(trial_old)|)/kt1, is
  subtracted so the ratchet is not re-counted — mohr_coul_direct3
  sigma_n = −200 = kn*(−1e-6 − 1e-6*tan(pi/4)), not −250). Feeds back
  into the normal history → gap/tension/max_fric of the next step.
- **Memory model** (`group_interface_materi_memory`): `-updated_linear`
  (default) computes the interface normal/tangent from `coord` (current
  configuration) each step; `-total_linear` reads the time-0 reference
  geometry (`NODE_START_REFINED` of side-1 nodes) so the interface keeps
  its original orientation even when the mesh deforms. Invalid values
  raise `db_error` (pattern `conspr.cc`). Read in `interface_element()`
  before the geometry block.

## Validación

- **Fase 1 elástica**: test 2 bloques (`/tmp/iface2.dat`) — kn=100→0.044,
  kn=1000→0.0018, kn=1e6→-0.003 (aprox soldado), kn=0.001→≈1.0 libre,
  sin interfaz→-0.0066. Límites físicos correctos.
- **Fase 3 validada (familia `iface_mc`, 13º test de build_safe.sh, 10 runs)**:
  - (a/a') fricción alta sostiene: phi=45°, c=0, Fy=5 < límite ≈ 80 →
    vely(nodo 6) ≈ 0 (0.20 en 10 pasos / 0.006 en 1 paso — invarianza de
    nº de pasos OK, mismos targets).
  - (b/b') fricción nula desliza: phi=0, c=0 con record PRESENTE → límite 0
    → `element_interface_force_tang`(2) = 0 (deslizamiento libre; fallback
    cinemático: vely=5 prescrita en los nodos 5-8).
  - (c) tracción abre: tension_limit=1.0, Fx=+10 → velix(nodo 6) = +3.88
    (strain_normal < 0; el código viejo exigía compresión y nunca abría).
  - (d) gap cierra con compresión: gap=0.001, Fx=-10 → strain 0.138 > gap →
    velix(nodo 6) ≈ 0 (−0.06; el código viejo abría bajo compresión).
  - (memory) `-total_linear` explícito: mismos targets que (a)
    (`iface_mc_mem`).
  - (dil/dil_1step) dilatancia RF-4: phi_flow=30° → strain_normal acumulado
    NEGATIVO (−2.8065 en 10 pasos / −2.82675 en 1 paso) frente a +0.0802
    con phi_flow=0 — la apertura por deslizamiento plástico supera la
    compresión. Diferencia 1-paso vs 10-pasos ~0.7% (dilatancia depende del
    |du_tang| incremental).
  - (num) clamp del MC numérico: c=50, phi=0, vely=1000 prescrita →
    `element_interface_force_tang`(2) = 50.0 exacto en todos los pasos
    (= max_fric = c; sin MC sería elástico ~2e5).
  - Discriminadores vs el código viejo: (b)/(b') elástico puro →
    force_tang ~10⁴ ≫ 1.0 FALLA; (c) nunca abre → velix −13.66 FALLA;
    (d) abre con compresión → velix grande FALLA.

## Ensamblaje

- The stiffness matrix is `[K -K; -K K]` on the displacement dofs of the
  two sides, with K = kn (normal) + kt*2 (tangential) projected onto the
  interface normal/tangent.
- Follows the `spring.cc` pattern: assembly on `vel_indx` (velocity dof).
  The nodal force is `-sign*stress*dir` (principle of virtual work), the
  matrix/lhside from `K*dtime`.
- **Sign convention (2026-08-13)**: the nodal force must be
  `-sign*(stress*dir)`; with the opposite sign the interface pushes
  instead of resisting (kn high increased the displacement).

## Validación

- Test 2 bloques (`/tmp/iface2.dat`): un bloque izquierdo fijo y un
  bloque derecho empujado, conectados por la interfaz. Verificado:
  - kn=100 → velix(nodo6)=0.044
  - kn=1000 → 0.0018
  - kn=1e6 → -0.0032
  - kn=0.001 → ≈1.0 (deslizamiento libre)
  - Sin interfaz (bloques soldados) → -0.0066
  El límite kn→∞ tiende al modelo soldado, kn→0 al deslizamiento libre —
  comportamiento físico correcto.
- La fuerza usa la velocidad relativa entre lados `(v_side2-v_side1)*dtime`
  (incremento de desplazamiento), no `dis_indx` (que es -1 con
  velocity_integrated).

## Pendiente / validación

- El test aislado (1 elemento de interfaz con lados prescritos) no valida
  bien porque los dofs prescritos no dejan que la interfaz frene.
- `control_mesh_convert` (bar2 -> quad4) implementado en la Fase 2 (commit
  `490545b`); el quad4 de interfaz también puede definirse manualmente.
- Solo 2D; el caso 3D (normal en z, 2 tangentes) está esbozado pero no
  validado.

## Detalles

- `normal`/`tangent` se calculan de la geometría: tangent a lo largo del
  lado 1, normal perpendicular (2D).
- La deformación usa el incremento `du_new - du_old` (acumulativo).

## Bugfix 2026-08-24: ddum3 sin inicializar en interface_element (NaN con -O1)

`db(GROUP_INTERFACE_MATERI_ELASTI_STIFFNESS, ..., ddum3, ..., GET_IF_EXISTS)`
no escribe `ddum3` cuando el record no existe, y el buffer local estaba sin
inicializar: en modelos que solo usan `group_interface_groundflow_*` (sin
rigidez elástica) `kn/kt1/kt2` quedaban con basura de stack. Compilado con
gcc 14.2 -O1 la basura contenía no-finitos y la matriz ensamblada daba NaN
(BI-CG "initial error -nan", test groundflow_interface). Con -O0 o con el
gcc del entorno anterior la basura no mordía — UB latente desde 2026-08-20.

Fix: `array_set(ddum3, 0., 3)` antes de la lectura. Regla general para todo
el codebase: inicializar SIEMPRE los buffers pasados a GET_IF_EXISTS.

## Sprint 8 additions (2026-08-24)

Three keywords completing the family (manual Professional 6.625/6.630/6.636):

- `group_interface_condif_conductivity`: q = k*(T1-T2) on the temp dofs
  of the facing node pairs with symmetric tangent — mirrors
  `group_interface_groundflow_permeability` exactly (interface.cc, right
  after the groundflow block; guarded by `condif_temperature`).
- `group_interface_materi_expansion_normal`: eps_th_total =
  alpha*T_avg (sides averaged) subtracted from the accumulated strain
  for the gap/tension/Mohr-Coulomb state (`strain_eff`); the INCREMENT
  d(alpha*T) is subtracted from du_norm so the normal FORCE this step
  is stiff*(du_norm - d_eps_th) — the same incremental eigenstrain
  pattern as stress.cc thermal strains. The stored
  ELEMENT_INTERFACE_STRAIN_NORMAL history stays purely mechanical (no
  thermal re-counting). Empirical: force per node pair = kn*alpha*T
  (test shows 2*kn*alpha*T sigxx for the 2-pair quad4).
- `group_interface_tangential_reference_point` (3D): t1 = perpendicular
  part of (ref - element centroid), t2 = n x t1; falls back to the
  default frame when the perpendicular part degenerates. Memory-model
  consistent (centroid from NODE_START_REFINED with -total_linear).

Also fixed the checklist: `group_interface_ground` never existed
(artifact); `_materi_memory` was already done (643865b);
`_elasti_sti`/`_residual_sti` are OCR truncations of `_stiffness`.

GOTCHA (bit me): ANY enum insertion in tochnog.h requires a FULL CLEAN
BUILD — incremental .o mixing old/new enum numbering corrupts the binary
with phantom check errors ("at least one of materi_velocity_integrated
or materi_displacement..." from mismatched ids). This is the documented
AGENTS.md build rule; now verified the hard way.

## Convergence record 2026-09-03 — interface Mohr-Coulomb direct (mohr_coul_direct3/4)

Changes in `interface_element()` (interface.cc) + `data()` (data.cc),
verified against the Professional binary .dbs (user-supplied 25-10-2023):

1. **Full-force assembly** (interface.cc, constitutive block + assembly):
   the element rhs now carries `stress_normal_ip = stiff*strain_normal_ip`
   (accumulated incl. this step + dilatancy) and the clamped accumulated
   trial `f_t` for the shear (both MC and elastic: `kt*gamma_total`).
   The records use the same accumulated values (rec_stress block
   unchanged, reads the same per-IP arrays). The matrix/tangent and the
   Lobatto-weighted per-pair assembly are unchanged.
2. **Incremental dilatancy**: `dgamma_inc = dgamma - gamma_pl_old` with
   `gamma_pl_old = (|f_el_old| - min(|f_el_old|, max_fric))/kt1`.
3. **Fn sign**: `max_fric_abs = |c - kn*strain_eff*tan(phi)|` (the
   accumulated strain is negative under compression, so this equals
   |c + |Fn|*tan(phi)|).
4. **`control_reset_dof -sigxx/-sigyy/-sigzz` → interface pre-stress**
   (data.cc, inside the CONTROL_RESET_VALUE_CONSTANT branch): after the
   node-dof reset, interface elements whose normal aligns with the reset
   axis (dot > 0.99 over the reference geometry) get
   `ELEMENT_INTERFACE_STRAIN_NORMAL := reset_value/kn` (all IPs, both
   versions). Verified: direct4 reset −sigyy −1 → the .dbs probe shows
   `element_interface_strain_normal 3 -1e-06 -1e-06` == the Professional
   epsilon_n = −1e-6.
5. **`control_reset_interface`/`_interface_strain` FIX** (data.cc): the
   reset block was nested inside `if (max_reset>=0)` (the control_reset_dof
   gate) AND bounded by `db_max_index(CONTROL_RESET_INTERFACE...)` which
   returns −1 although the records are active (measured: record 15 active,
   max index −1) AND wrote one value with a leftover length (only
   VERSION_NORMAL). Now: unconditional scan of the control range with
   `db_active_index` (0..1000), fires at its own control index, zeroes ALL
   ns1 slots of ELEMENT_INTERFACE_STRAIN_NORMAL / FORCE_TANG(_2) in BOTH
   versions. Without this fix interface7 (2 one-step phases separated by a
   reset at index 15) accumulates −2 instead of −1 with the new
   accumulated records.

Blast radius (corpus 363): 141 PASS (baseline 140) / 200 RUNFAIL /
22 PARSE; interface1/7/8/9/12/14/15/patch + conspr1-7 rc=0; direct3 rc=0
(σn −199.99997, node_rhside ±99.99998 — identical to the Professional).
mohr_coul_direct4 = RUNFAIL blocker: GNU 0.251052 vs target 0.10244 at
t=100. The GNU reaches the static plastic equilibrium (u3 = cap/4 = 0.25
with cap ≈ c = 1.0 since the kinematic σn relaxes to ~0 as the blocks
drift) while the Professional's σn = −1 (the −sigyy reset) is held for
the whole 100 s run (its stress-dof state persists; u3 relaxes to
0.10244). Reproducing it requires the reset to act as a PERSISTENT
pre-stress on the block + interface stress-dof mechanics, not just an
initial condition — out of scope of the interface family.

## Convergence record 2026-09-04 — gap multi-step, reset semantics, thermal (interface2/10/expans3)

Changes in `interface_element()` (interface.cc), `data()` (data.cc,
resets), `print_interface_stress.cc` and the new history
`ELEMENT_INTERFACE_FORCE_NORM` (tochnog.h/tochnog-mod.h enum sync,
database.cc registration with the FORCE_TANG pattern: DOUBLE, length 4,
fixed_length 0, version_all 1, print_only 1, class/required ELEMENT;
top.cc db_allocate VERSION_NEW under `any_interface`). Verified against
the Professional per-step prints and `.dbs` (25-10-2023):

1. **Gap condition inverted** (interface.cc): closed when the
   accumulated normal strain `<= gap`, open when `> gap` (manual 6.625);
   default gap +1e20. Pre-fix code opened when `strain <= gap` with a
   −1e20 default (correct only for tests without a gap record). patch1
   (gap 1.e20) was stuck ALWAYS OPEN → now closed.
2. **Closed-phase stress history** `ELEMENT_INTERFACE_FORCE_NORM`: per
   IP, per step: closed → `+= kn*(delta_n - thermal_inc_step)` (delta_n
   measured against the pre-step stored strain, so the dilatancy opening
   of the MC block is included); open → `= 0`. `stress,normal_ip` (record
   + rhs) = this history. Equals `kn*strain,normal` for the always-closed
   tests (bit-identical up to FP round-off) and reproduces the
   Professional `interface2` trajectory exactly (−1 at the closing step
   100, −101 at 200).
3. **`control_reset_interface_strain` semantics** (data.cc): per manual
   6.355 it zeroes the accumulated STRAINS but REMEMBERS the stresses.
   The stress lives in `ELEMENT_INTERFACE_FORCE_NORM`, so the strain
   reset now zeroes `ELEMENT_INTERFACE_STRAIN_NORMAL` only and leaves
   FORCE_NORM untouched (pre-fix the code zeroed the strain history and
   the next step re-compressed the interface — interface10: displacement
   −6e-4 → −1.2e-3, strain −6e-4 instead of 0). `control_reset_interface`
   (full) additionally zeroes FORCE_NORM + FORCE_TANG(_2). The
   `control_reset_dof -sigxx/-sigyy/-sigzz` pre-stress additionally sets
   FORCE_NORM := the reset stress value.
4. **Thermal expansion along the normal + into the stress**
   (interface.cc): the thermal pseudo-load subtracts the increment along
   the interface NORMAL (`du_ip -= alpha*dT*n_hat`, pre-fix: x component
   only → spurious tangential slip −alpha*dT/sqrt(2) on a 45-degree
   interface); the closed-phase stress accumulates `kn*(delta_n −
   alpha*dT)` so a heated constrained interface carries
   `kn*(-alpha*T)` (expans3: sigma_n = −1 = 1*(−1)). The strain history
   stays purely mechanical; `strain_eff = n − alpha*T_total` keeps
   driving the gap/tension/MC state.
5. **Tension status record** (interface.cc output block): CLOSED when
   `strain_eff <= gap`, OPENED otherwise (was inverted).

Blast radius (corpus 363): 150 → 153 PASS (interface2, interface10,
expans3 rc=0; interface13 stays RUNFAIL as a corpus-test bug, see the
user manual — both codes compute 0.67082 while the target demands
1.11803 = |du|/2 and the Professional itself reports "Error detected");
interface1/7/8/9/12/14/15/patch + conspr1-7 rc=0 unchanged.
`patch1` = RUNFAIL blocker sharpened: gap fixed (sigxx 1199.49 →
1200.01, the interface transmits again) but its ±1e-3 ABSOLUTE targets
on a kn=1e11 penalty system need direct-solver precision: the default
Bi-CG leaves ~1e-6 relative solution error, which kn=1e11 amplifies to
hundreds of kPa in sigma_n (GNU 1431 per-intpnt vs Pro 960 uniform),
and the SuperLU sparse path segfaults on this model (so_suplu.c) —
pending solver work, same root as interface_bar2_hex8.
`interface11` = mesh_interface_triangle_* family (manual 6.854/6.855/
6.201): generating interface elements by intersecting a triangulated
plane with a tet4 mesh — not registered, real mesh-generation feature
(out of sprint scope, documented).

## Convergence record 2026-09-04 — the interface element measure (patch1)

The assembled per-pair spring force and stiffness must carry the element
MEASURE: `force = w_i * L * sigma` and `stiffness = w_i * L * kn` with
L = the element length (2D line) or the face area (3D surface). The
historical assembly used the bare Lobatto/even weights (sum 1 = the
UNIT-measure element), which is exact only for the unit-length interfaces
(interface1: length 1) and silently under-integrates every other element.

patch1 of the corpus (two inclined quad6 interfaces of length 1.677/0.559
between loaded quad9 blocks) converged to sigma_n = 1610/1073/536 per
intpnt (mean 1431) instead of the Professional's uniform 960: without the
length the discrete force system loses the load-path moment arm (the
length-weighted centroid of the pairs = the load line only WITH the
measure) and the uniform traction cannot balance the applied edge load.
The fix (2026-09-04): `iface_measure` in `interface_element()` — the 2D
chord between the first and the last side-1 node (reference geometry per
the memory branch: -total_linear → NODE_START_REFINED), the 3D face area
by the triangle fan (the tochnog hex8/quad4 face numbering is a bowtie
order — the crossed-diagonals formula gives 0). With the measure the
per-IP sigma_n = 960.0000000 with the direct solver (rc=0); the
unit-length validations are unchanged (rc=0).

PENDING (same lote): the mpc_linear_quadratic tying (mpc3/4: 0.299/0.350
vs 1/3 — the generated ties constrain the velocity dofs but the mixed
σ-dofs of the tied nodes stay free → non-homogeneous field), the
phreatic-multiple + mechanics coupling (ground8), the materi_dynamic
explicit limit (dynamic1/2/5/8).
