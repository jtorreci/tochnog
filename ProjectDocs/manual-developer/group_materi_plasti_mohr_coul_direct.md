# group_materi_plasti_mohr_coul_direct / tension_direct (+ _normal, _normal_automatic)

## Implementación (dos modos)

The direct records have TWO modes selected by the presence of a plane normal:

1. **Full principal-stress mode** (NO `_normal` / `_normal_automatic`):
   `materi_direct_full_mc()` in `stress.cc`, dispatched from `set_stress()`
   right after the elastic stress and BEFORE the plastic-yield test. It
   implements the manual Professional 6.726 + 6.738:
   - spectral tension cap: `matrix_jacobi` eigenvalues above `sigy` are cut
     to `sigy` (6.738; `sigy` defaults to 0 when `tension_direct` is absent
     but `mohr_coul_direct` is present);
   - Mohr-Coulomb principal-stress-difference cut: sorted eigenvalues
     `w0<=w1<=w2`, `f = 0.5(w2-w0) + 0.5(w2+w0) sin(phi) - c cos(phi)`; when
     `f > 0` a ONE-SHOT return runs along the non-associative flow direction
     `deps = (0.5(1+sin psi), 0, -0.5(1-sin psi))` (psi = phi_flow) projected
     with the isotropic elastic C (`lambda_lame`, `gmod` from the group
     elastic data): `w2 -= f*cd1/denom`, `w0 += ...`, and the middle
     principal gets the C-coupling `-f*cd2/denom` (`lambda_lame*sin psi`),
     with the corner rule `w1 = min(w1, new w2)` (the sigma1 = sigma2 edge
     of the surface). Eigenvector tracking: the eigenvalues are sorted with
     an order[] index so the rebuild `sigma = V diag(w) V^T` keeps the
     eigenvector columns aligned.
   - `_visco tm` relaxes the whole correction with `1-exp(-dt/tm)`; `_wall`
     replaces phi/c/sigy when `plasti_on_boundary`.

2. **Plane traction mode** (a law with `_normal` / `_normal_automatic`):
   `materi_direct_cutoff()` in `stress.cc`. Rewritten 2026-10-02 (C14,
   commits `d17461e` + `77f9282`) to the rule measured on the
   Professional 02-08-2026:
   - **Normals.** `materi()` (`materi.cc`) builds one normal PER LAW with
     `direct_normal_get()` (end of `materi.cc`, prototype in `tochnog.h`
     and `tochnog-mod.h`): explicit `_normal` record, else `ELEMENT_NORMAL`
     when that law's own `_normal_automatic` is `-yes`, else zero. Both go
     through `set_stress()` (`mc_normal`, `ten_normal`); when both are
     zero `set_stress()` calls `materi_direct_full_mc(..., 1)` (mode 1),
     otherwise `materi_direct_cutoff()`. Before C14 a single
     `direct_normal` (MC normal, else tension normal, else element
     normal) served both laws.
   - **Tension step first.** `sigy` = `tension_direct` (wall value when
     `plasti_on_boundary`), else 0 when `mohr_coul_direct` exists (6.738;
     same logic as `materi_direct_full_mc`, including the
     `control_materi_plasti_tension_apply -no` gate, which drops the
     record but keeps the default 0). With a tension normal: cap
     `sig_n = n.sig.n <= sigy` with `sig -= (sig_n - sigy) f_ten n x n`.
     Without it: `materi_direct_principal_cap()` (jacobi + the shared
     helper `direct_cap_principal_values()`, which is also step 1 of
     `materi_direct_full_mc`; no visco factor there, as in mode 1).
   - **MC step.** With an MC normal: traction `t = sig.n`, `sn = n.t`,
     `tau = t - sn n` from the stress AFTER the tension step. If
     `|tau| > max(c - sn tan(phi), 0)`:
     `k = direct_flow_factor(C, n, tau, |tau|) * tan(phi_flow)` (0 when
     phi_flow = 0), `tau_new = max((c - sn tan(phi) + k tan(phi) |tau|) /
     (1 + k tan(phi)), 0)`, `dsn = -k (|tau| - tau_new)`; the stress gets
     the tangential scaling `-(1-scale)(tau x n + n x tau)` with
     `scale = 1 - f_mc (1 - tau_new/|tau|)` and `+ f_mc dsn n x n`. Only sn
     and tau change (the Professional does not update the other
     components). Without an MC normal (tension plane only):
     `materi_direct_full_mc(..., 0)`, the full MC cut without its tension
     cap (verified: Pro probes c21/c21b/c21d equal mode 1 for the MC part).
   - **M/G source.** `direct_flow_factor()` (`stress.cc`) evaluates
     `M = n_i n_j C_ijkl n_k n_l` and
     `G = 0.5 t_i n_j C_ijkl (t_k n_l + n_k t_l)` (t = tau/|tau|; the
     symmetrised form is independent of how C stores the shear minor
     symmetry: `C_matrix()` sets `C_0101 = 2G`, `C_0110 = 0`). `C` is
     built in `set_stress()` at the call site as
     `Cmem + Cuser + Chyper + Chypo` plus `1/compressibility` on the
     `C_iikk` terms: the same linear stiffness that fills `ddsdde` at the
     end of `set_stress()`, so plane stress (`group_materi_membrane`),
     transverse isotropy, `shear_factor`, hardsoil/young_power stiffness
     (computed into Cmem before the cut) are all honoured. Chosen over
     reading young/poisson again (the approach of `materi_direct_full_mc`)
     because it is general and needs no per-law branch. `group_materi_factor`
     scales C after the cut and does not change the ratio.
   - **Hardcoded**: fallback `M/G = 2` (the Professional constant) when
     `M <= 0` or `G <= 0` (umat writes `ddsdde` directly, not a C tensor);
     guard `1 + k tan(phi) < 1e-12 -> k = 0`; visco factors now separate
     per law (`f_mc`, `f_ten`; before C14 one `tm` was shared and the
     tension value won).
   - `group_materi_plasti_bounda/_factor` (6.231/6.232) reduce phi and c
     (not phi_flow) on the wall.
   - **Input**: `GROUP_MATERI_PLASTI_TENSION_DIRECT_NORMAL` no longer
     requires `GROUP_MATERI_PLASTI_TENSION_DIRECT` (`database.cc`,
     `data_required = GROUP_TYPE`; the Professional accepts the normal
     alone and applies nothing without either law, probe c20). The
     `_automatic` switch keeps the prerequisite (the Professional rejects
     it, probe c23).
   - `_normal_automatic` reads the generated `ELEMENT_NORMAL` record
     (2026-09-30, `element_normal.md`). A zero `ELEMENT_NORMAL` in a
     `-yes` group without any explicit `_normal` is an input error raised
     in `element_normal_set()` (`area.cc`, helper
     `element_normal_automatic()`; Pro: assertion). A group switched `-no`
     keeps the law's normal zero (test
     `element_normal2d_column_switch_no1`).
   - The zero-normal check is per law since `b7f3522`:
     `element_normal_automatic()` (area.cc) skips a law only when that
     same law has an explicit `_normal`, so an explicit MC normal plus
     `tension_direct_normal_automatic -yes` on a zero `element_normal` is
     an input error (test `element_normal2d_row_mixed_error1`; Pro:
     assertion).
     Anisotropic coupling terms (normal strain -> shear stress) are not
     included in the return, which uses only M and G.

## Física

- Full mode (6.726): "Principal stress differences higher than allowed by
  the mohr-coulomb criterium are not allowed and will be cut off by
  Tochnog" — the alternative programming of the MC law. The cut does not
  use plastic strains; phi_flow enters the cut DIRECTION (the non-assoc
  flow ratio `(1+sin psi)/(1-sin psi)` of the max/min corrections; for
  psi=0 the cut is the mean-preserving difference cut).
- Plane mode (6.727/6.739): limits the friction/tension stress on the
  specific plane (interface semantics). phi_flow: the slip plastic strain
  `deps = dlambda (sym(t x n) + tan(psi) n x n)` gives `dtau = -G dlambda`
  and `dsn = -M tan(psi) dlambda`, hence `dsn = (M/G) tan(psi) dtau`. The
  Professional uses 2 (= M/G at nu = 0): deliberate divergence, BACKLOG
  C15.

## Differences with the Professional version

- State output (K13, 2026-10-06): `set_stress()` stores the cut as
  plastic strain and the yield value (`direct_plastic_strain.md`);
  shear components in tensor form (Professional: twice).
- phi_flow factor of the plane mode: M/G instead of the Professional
  constant 2 (deliberate, BACKLOG C15; see Física).
- Interface strength (`mohr_coul_direct4`, deliberate, BACKLOG C30,
  project decision 2026-10-06): the GNU returns the interface friction
  at once to the manual 6.632 strength c + |Fn| tan(phi) (1.2027 in that
  test); the Professional 02-08-2026 keeps the elastic trial (2.0) and
  drifts slowly below the strength (0.410 at t = 100). Not imitated: no
  code change, the manual decides.

- **`_visco` direct cut-offs** (known difference, not decided; BACKLOG C33).
  GNU: Relaxation factor uses the previous step's `DTIME` (no cut in the
  first step); no visco on the principal tension cap. Professional
  02-08-2026: 1 - exp(-dt/tm) with the current step (sigxx 1.84914e-3 vs GNU
  1e-3). Status: Scheduled for F4, project owner 2026-10-06. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **3 values in `group_materi_plasti_*_direct_normal` in 2D** (GNU extension
  / input dialect; manual pages only). GNU: Accepts up to 3 values in any
  dimension (`data_length` 3, `fixed_length` 0 in `database.cc`); the use of
  a z value in 2D is not verified. Professional 02-08-2026: Manual 6.739:
  ndim values ("In 1D only normal_x"); binary response to 3 values in 2D NOT
  measured. Status: Input-dialect tolerance; Pro behaviour undetermined.
  Index: [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).

## Validación del modo plano (Professional 02-08-2026, C14)

Probes (single quad4, E 1, homogeneous velocity, 1 step; scratch driver
`c14/run.py`): closed form vs Professional vs old GNU vs new GNU, node 1
stresses (sxx, sxy, syy). The new GNU differs from the closed form by
0.1-0.2 % with the default memory (C16) and reproduces it exactly with
`-total_linear` (dedicated tests).

| probe | closed form | Professional | old GNU | new GNU |
|---|---|---|---|---|
| c01 MC c 1e-4, no ten | -2.5e-4, 1e-4, -2.5e-4 | same | -2.5e-8, 1e-4, 2.75e-7 | -2.49962e-4, 1e-4, -2.49913e-4 |
| c04 MC c 1e-2, syy tension | -2.45145e-4, 2.00971e-4, -1.64757e-4 | same | 3.7e-7, 5.0005e-4, 1.99875e-4 | -2.44970e-4, 2.00930e-4, -1.64807e-4 |
| c13 ten 1e-4, no ten normal | -2e-4, 3e-4, -2e-4 | same | 3.75e-7, 5e-4, -1.25e-7 | -1.99788e-4, 2.99937e-4, -2.00087e-4 |
| c16 ten 0 + ten normal | 3e-4, 0, 0 | same | same | same |
| c17 ten normal only | 3e-4, 0, 0 | same | input error | 3e-4, 0, 0 |
| c19 psi 0.3, nu 0 | 0, 6.60909e-4, -1.20979e-3 | same | 5.9e-7, 5.44483e-4, -9.99590e-4 | 8.2e-7, 6.59212e-4, -1.21022e-3 |
| c15 psi 0.3, nu 0.25 (M/G 3) | -4e-4, 7.04154e-4, -1.28895e-3 | -4e-4, 6.92049e-4, -1.26679e-3 | -3.98691e-4, 6.53888e-4, -1.19971e-3 | -3.98593e-4, 7.02818e-4, -1.28954e-3 |
| c22 beyond apex, nu 0 | 0, 0, 9.38133e-4 | same | 0, 1e-7, 1.00001e-3 | 0, 9.4e-8, 9.38179e-4 |
| c24 ten plane x, MC plane y | 0, 1e-4, 0 | same | input error | 1e-7, 1e-4, 2.7e-7 |
| c21b MC no normal + ten normal x | full MC (Pro -6.69503e-4, 1.69503e-4, -3.30497e-4) | — | plane MC on x | -6.69250e-4, 1.69545e-4, -3.30500e-4 |

Dedicated tests `validation-suite/test-2026/direct_cutoff_plane_*.dat`
(`-total_linear`, targets at 1e-8 from the closed form): shear1,
tension_principal1, tension_value1, tension_normal1, tension_normal2,
two_normals1, phi_flow1, phi_flow2, phi_flow_floor1. Old binary: 8 RED
(tension_normal1 passes on both, non-regression); after `d17461e`: the
three phi_flow tests RED; after `77f9282`: 9/9 GREEN. The Professional
passes all except phi_flow2 (C15).

## Validación (Professional 25-10-2023)

- mohr_coul_direct1 (oedometer, phi=0.4 c=1 psi=0.2): PASS — eptxx target
  1.49e-2 reached (the flow ratio reproduces the dilative lateral strain).
- mohr_coul_direct2 (direct shear, phi=0 c=0 psi=pi/4): PASS — vely = velx.
- mohr_coul_direct6/7 (single-increment 100% shear, tension cap at 0):
  the material map is EXACT (sigma_xy = a*sin(phi), sigma_zz =
  -a*(1-sin(phi)) with a the capped shear) when the element kinematics are
  linear (`group_materi_memory -total_linear`: sigxy = 0.0249584 EXACT for
  direct6). The corpus runs keep the GNU DEFAULT `-updated` element strain
  (incremental polar decomposition U-I), which differs from the
  Professional's linear response at the gamma=1.0 single increment -> the
  two tests stay RUNFAIL on the ELEMENT kinematics (not the material law).

## Keywords

- The `_normal` records are now `fixed_length = 0` (ndim values: 1 in 1D,
  2 in 2D, 3 in 3D — manual 6.727 "In 1d only specify normal_x").

## group_materi_plasti_bounda wall detection (group.cc)

`plasti_on_boundary()` (group.cc) keeps the legacy element-group semantics
AND adds the Professional bounda semantics (6.231): when a listed value
matches an ACTIVE `bounda_dof` record, the element is on the wall when one
of its nodes is bounded (`node_bounded`) on the velocity/displacement
parts.

## Verificación runtime de la reducción de pared (2026-09-07, direct8)

La semántica de pared (6.231/6.232) está ACTIVA y es correcta en runtime:
A/B discriminante con el material directo de plano (normal (0,1)) en
`-total_linear` con COMPRESIÓN + cortante simultáneos (el estado puro de
cortante de direct8 no discrimina: con c=0 y sig_n=0 el límite de
fricción es 0 con o sin factor):

- elástico: sigxy = 2500 (E=1e4, nu=0, vx=0.5), sigyy = -1e4.
- ley directa SIN plasti_bounda: sigxy = 2500 SIN corte — la fricción
  aguanta: max_fric = c - sig_n*tan(phi) = 0 + 1e4*0.4228 = 4228 > 2500
  (la compresión AUMENTA el límite, física correcta).
- ley directa CON `group_materi_plasti_bounda 0 10` (10 = el índice del
  record bounda_dof 10 que fija los nodos 1-2) y
  `group_materi_plasti_bounda_factor 0 0.0`: sigxy = 0 — el factor 0
  anula la fricción de pared (max_fric = 0).

La detección de pared (group.cc `group_materi_plasti_boundary_evaluate`,
mecanismo bounda_dof: el elemento está en la pared cuando un valor del
record coincide con un bounda_dof ACTIVO y uno de sus nodos está bounded
en las partes de velocidad/desplazamiento) y la reducción del factor en
`materi_direct_cutoff` (stress.cc, phi y c por el factor) funcionan.

### mohr_coul_direct8 (corpus): estado y blocker

direct8 parsea y corre; la reducción de pared se aplica (ver A/B); el
test queda RUNFAIL por la CINEMÁTICA del elemento, no por la ley: el .dat
NO lleva `group_materi_memory` → default GNU `-updated` (deformación
finita con rotación polar). Con el incremento único de cortante 100 %
(vx=1 durante 1 s sobre altura 1, γ=1) la rotación incremental (~26.57°)
reintroduce la cizalla espacial: el cutoff anula sigxy en el frame
material, pero al rotar el tensor al frame espacial σxy = 1788.85 ≠ 0
(target 0). Con `-total_linear` (A/B, añadiendo materi_displacement) el
test PASA (sigxy = 0, rc=0). Es el mismo blocker de cinemática de
direct6/7 (γ=1, default -updated vs la respuesta lineal del
Professional). El Professional resuelve este modelo sin rotación (sus
nodos no se mueven y su σyy queda congelada en el reset −30), lo que
apunta a que su default para modelos velocity-only sin desplazamiento no
rota — decisión de cinemática del elemento, fuera del alcance de la ley.

## mohr_coul_direct5 (force_edge sobre bar2 axisimétrico 1D)

direct5 necesita cargas de borde (force_edge) sobre una malla 1D de
elementos bar2 axisimétrica generada por `control_mesh_macro -bar`
(1000 barras radiales entre r=1 y r=41). Ver el manual
force_element_edge_water.md / el cambio en area.cc: los bar2 ahora tienen
2 "lados" degenerados (cada extremo = un nodo), el área del lado = 1 y la
carga nodal del anillo axisimétrico = value*2*pi*r. El parse y el
análisis completo funcionan; los targets quedan RUNFAIL por el
acoplamiento de la mecánica mixta 1D (la tensión transversal σyy del GNU
se relaja vía Poisson hasta −12.6 mientras el Professional la mantiene
congelada en el reset −30, y la integración temporal acumula ~1.58× el
desplazamiento del Pro) — diagnóstico fino en el reporte del sprint.
