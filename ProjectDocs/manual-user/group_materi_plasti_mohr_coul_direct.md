# group_materi_plasti_mohr_coul_direct / tension_direct

## Description

`group_materi_plasti_mohr_coul_direct` and `group_materi_plasti_tension_direct`
are the **direct stress cut-off** variants of the Mohr-Coulomb and tension
plasticity laws (manual Professional 6.726 / 6.738). They are direct cut-offs
of the stresses (no plastic strains, no incremental return mapping), which
makes them very stable.

The record has TWO modes, selected by the presence of the plane-normal
records:

### 1. Full principal-stress cut-off (no `_normal`, no `_normal_automatic -yes`)

`group_materi_plasti_mohr_coul_direct phi c phi_flow` caps the PRINCIPAL
STRESS DIFFERENCES with the Mohr-Coulomb criterion

    f = 0.5(sig1 - sig3) + 0.5(sig1 + sig3) sin(phi) - c cos(phi) <= 0

(sig1 largest, sig3 smallest principal stress, tension-positive). Stresses
outside the surface are cut back onto it. `group_materi_plasti_tension_direct
sigy` caps the principal stresses: "principal stresses higher than sigy are
not allowed and will be cut off by Tochnog" (6.738). If `tension_direct` is
not specified but `mohr_coul_direct` is available, sigy is set to 0 (6.738).

The cut follows the non-associative flow direction given by `phi_flow`
(dilatancy): for phi_flow = 0 the difference cut is mean-preserving
(isochoric), for phi_flow > 0 the pair mean shifts into compression
(dilatant flow). Verified against the Professional binary 25-10-2023:

- the shear-only trial (after the tension cap) maps to
  `(-a(1-sin), -a(1-sin), -a(1+sin))` in the principal frame (a = capped
  shear) — the state lands on the sigma1 = sigma2 edge of the surface;
- the oedometer (mohr_coul_direct1) reproduces the full plastic lateral
  strain with the phi_flow dilatancy;
- the direct shear with phi = 0, c = 0 and phi_flow = pi/4 (mohr_coul_direct2)
  dilates with disy = disx.

The material counterpart of the incremental
`group_materi_plasti_mohr_coul` law: same surface, "alternative programming
of the mohr-coulomb law, which tends to be very stable" (6.726).

### 2. Plane traction cut-off (a law with `_normal` / `_normal_automatic -yes`)

As soon as one of the two laws has a plane normal, each law works on its
OWN plane (manual 6.727/6.739; rule verified against the Professional
02-08-2026 on 2026-10-02, C14). The tension law is applied first, then
the Mohr-Coulomb law:

1. **Tension.** `sigy` is the `tension_direct` value, or 0 when only
   `mohr_coul_direct` is given (6.738). If the tension law has its own
   normal (`group_materi_plasti_tension_direct_normal`, or
   `_tension_direct_normal_automatic -yes`) the normal traction on that
   plane is capped: `sig_n <= sigy`. Otherwise the PRINCIPAL stresses are
   capped at `sigy`, exactly as in the full mode. Consequence: an input
   with `mohr_coul_direct_normal` and no `tension_direct` now cuts the
   principal tension at 0 (before 2026-10-02 the GNU applied no tension
   cut at all in that case).
2. **Mohr-Coulomb.** On the MC plane (`_mohr_coul_direct_normal` or
   `_normal_automatic -yes`), with the normal traction `sn` and the shear
   traction `tau = |tau|` of the stress after the tension step: when
   `tau > max(c - sn*tan(phi), 0)` the state returns along the
   non-associative flow direction (`phi_flow` = psi). Only `sn` and `tau`
   change:

       k       = (M/G) tan(psi)
       tau_new = max( (c - sn tan(phi) + k tan(phi) tau) / (1 + k tan(phi)), 0 )
       sn_new  = sn - k (tau - tau_new)

   `M` is the constrained (oedometric) modulus normal to the plane and `G`
   the shear modulus on the plane, both from the elastic stiffness of the
   element (isotropic: `M = E(1-nu)/((1+nu)(1-2nu))`, `G = E/(2(1+nu))`;
   plane stress and the other elasticity laws are taken from the actual
   stiffness). With `phi_flow = 0` the normal traction is unchanged and the
   shear is simply capped. Beyond the apex (`c - sn tan(phi) < 0`) the
   shear goes to 0 and `sn_new = sn - k tau`. If the MC law has no normal
   while the tension law has one, the MC law is the full principal cut-off
   of mode 1 (without its tension cap, already applied on the plane).

`group_materi_plasti_tension_direct_normal` may be given without
`group_materi_plasti_tension_direct` (as in the Professional): with
`mohr_coul_direct` the tension plane then uses `sigy = 0`; without either
direct law the record has no effect. `_tension_direct_normal_automatic`
still requires `tension_direct` (the Professional rejects it too).

`_normal_automatic -yes` takes the plane normal from the element normal
`element_normal` (inward normal of the mesh boundary, see
`element_normal.md`; before 2026-09-30 the GNU used the cross product of
two element edges, which is out of plane in 2D and made the cut-off
inert). `_normal_automatic -no` means no automatic normal: the group
stays in the full principal-stress mode, as in the Professional (until
2026-10-02 the GNU still used the element normal with `-no`).

An element of a `_normal_automatic -yes` group whose `element_normal` is
zero (interior element, single element, middle element of a single row)
stops the calculation with an input error that names the element and the
group and suggests the remedies: restrict the group to the boundary
interface layer, or give the normal explicitly with
`group_materi_plasti_tension_direct_normal` /
`group_materi_plasti_mohr_coul_direct_normal` (for an element resting on
a supported face, the normal of that face). See `element_normal.md` for
the full message. The Professional also stops on such an element, with an
internal assertion (`!array_null( directions, MDIM*MDIM )`) instead of a
readable message. In 2D give the explicit normal with two values
(`group_materi_plasti_tension_direct_normal 0  0. 1.`): the Professional
rejects a third value in 2D.

## Uso

Place the records in the data part, in the element group definition:

```
group_type 10  -materi
group_materi_elasti_young 10  1000.0
group_materi_elasti_poisson 10  0.3
group_materi_plasti_tension_direct 10  1.0
group_materi_plasti_mohr_coul_direct 10  0.785398  1.0  0.0
```

Without a plane normal the MC direct law needs the tension direct record
(`sigy` defaults to 0 when missing). For the plane variant:

```
group_materi_plasti_mohr_coul_direct 10  0.785398  1.0  0.0
group_materi_plasti_mohr_coul_direct_normal 10  0. 1. 0.
```

Tension and friction on two different planes (2D; tension plane x, MC
plane y; sigy = 0 because `tension_direct` is absent):

```
group_materi_plasti_mohr_coul_direct 10  0.  1.e-4  0.
group_materi_plasti_mohr_coul_direct_normal 10  0. 1.
group_materi_plasti_tension_direct_normal 10  1. 0.
```

The normal record accepts ndim values (1 value in 1D, 2 in 2D, 3 in 3D).

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_plasti_mohr_coul_direct` | `phi c phi_flow` | Friction angle (rad), cohesion, dilatancy angle (rad; in the plane mode it shifts sn by k(tau - tau_new), k = (M/G) tan(phi_flow)). |
| `group_materi_plasti_mohr_coul_direct_normal` | `nx ny nz` | Plane normal (explicit; ndim values). Selects the plane mode. |
| `group_materi_plasti_mohr_coul_direct_normal_automatic` | `-yes`/`-no` | `-yes`: take the normal from the element normal. `-no`: no plane normal (full cut-off). |
| `group_materi_plasti_mohr_coul_direct_visco` | `tm` | Visco relaxation time. |
| `group_materi_plasti_mohr_coul_direct_wall` | `phi c phi_flow` | Values used when the element is attached to a wall. |
| `group_materi_plasti_tension_direct` | `sigy` | Tensile limit on the principal stresses (full mode) or on the plane normal traction (plane mode). |
| `group_materi_plasti_tension_direct_normal` | `nx ny nz` | Plane normal of the tension law (explicit). Allowed without `tension_direct` (sigy = 0 with `mohr_coul_direct`). |
| `group_materi_plasti_tension_direct_normal_automatic` | `-yes`/`-no` | `-yes`: take the normal from the element normal. `-no`: no plane normal (full cut-off). |
| `group_materi_plasti_tension_direct_visco` | `tm` | Visco relaxation time. |
| `group_materi_plasti_tension_direct_wall` | `sigy` | Tensile limit used when the element is attached to a wall. |

## Notas

- Requires `materi_stress` in the initialization part.
- The cut-off is applied on the elastic stress BEFORE the incremental
  plastic-yield test, so it combines with the incremental plasticity models.
- `_visco tm`: viscoplastic relaxation. The stress interpolates between the
  elastic and the fully capped (plastic) response with
  `factor = 1 - exp(-dt/tm)`. `dt << tm` -> elastic, `dt >> tm` -> fully
  capped.
- `_wall`: alternative parameters used when the element is attached to a
  wall (its nodes are bounded by a `bounda_dof` record listed in
  `group_materi_plasti_bounda`, manual 6.231).
- `group_materi_plasti_bounda index bounda_indices...` reduces the friction
  (phi, c, phi_flow) of the elements on the wall by the factor
  `group_materi_plasti_bounda_factor` (default 2/3, manual 6.231/6.232) —
  verified on mohr_coul_direct8 (factor 0 -> zero friction stress on the
  wall elements).

## Differences with the Professional version

- **Plastic strain of the cut (deliberate for the shear components,
  K13).** Since 2026-10-06 the cut is stored in `materi_strain_plasti`
  as C^-1 (sigma_trial - sigma) and `materi_plasti_f` holds the direct-law
  yield value (0 after a cut), as the Professional; the Professional
  writes the shear components as twice the tensor value, the GNU keeps
  the tensor value. See `direct_plastic_strain.md`.
- **phi_flow factor in the plane mode (deliberate, BACKLOG C15).** The
  Professional 02-08-2026 uses the constant `k = 2 tan(phi_flow)`; the GNU
  uses `k = (M/G) tan(phi_flow)`, the consistent elastic return of the
  slip plastic strain. Both are equal for `poisson = 0` (M/G = 2). Example
  (plane strain, E = 1, nu = 0.25, phi 0.5, c 0, phi_flow 0.3, MC normal
  (0,1), trial sxy 8e-4, syy -1.2e-3): GNU sxy 7.04154e-4, syy
  -1.28895e-3; Professional sxy 6.92049e-4, syy -1.26679e-3
  (`validation-suite/test-2026/direct_cutoff_plane_phi_flow2.dat`). When
  no elastic stiffness is available (e.g. a `umat` law) the GNU falls back
  to the Professional constant 2.
- **Stale target `materi_direct_mc` (validation-suite).** Its two free
  y-nodes plus the principal tension cap make the element collapse at
  t = 0.7 (sigxy 0.13 -> 0.06); the Professional collapses at the same
  step with the same constraints (sigxy ~0), so its target 1 +- 0.3 only
  held with the pre-2026-10-02 rule.
- **Interface strength `mohr_coul_direct4` (deliberate, BACKLOG C30).**
  Interface shear test (phi 0.2 rad, c 1, normal stress -1, top displaced
  1 at t = 1 and held to t = 100). The manual 6.632 limits the friction
  to c + |Fn| tan(phi) = 1.2027; the GNU returns to it at once (sigxy
  1.2027 from t = 1, `disx` 0.3007, status plastic). The Professional
  02-08-2026 stays at the elastic 2.0 (status `-plastic`) and drifts down
  slowly (t = 5: 1.991, t = 50: 1.095) to 0.410 `-elastic` at t = 100
  (its target `disx` 0.1024). The drift looks like an artefact of the
  Professional; the GNU follows the manual (project decision 2026-10-06).

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

## Estado de implementación

- **Implementado**: the full principal-stress mode (spectral tension cap +
  MC difference cut with the non-associative flow) and the plane mode
  (`_normal`, `_normal_automatic`, `_visco`, `_wall`; separate tension and
  MC planes, tension first with the default sigy = 0, phi_flow return,
  2026-10-02). Dedicated tests:
  `validation-suite/test-2026/direct_cutoff_plane_*.dat` (9, closed form;
  8 equal the Professional, `_phi_flow2` documents the divergence).
- **Validación**: mohr_coul_direct1/2 pass; the single-increment shear-only
  states (mohr_coul_direct6/7) reproduce the Professional EXACTLY when the
  element kinematics are the linear ones (`group_materi_memory -total_linear`);
  the corpus runs of direct6/7 stay RUNFAIL because their single 100%-shear
  increment exercises the GNU's default `-updated` (polar decomposition)
  element strain, which differs from the Professional's linear response at
  that magnitude (element kinematics, not the material law).
- **Pendiente**: softening via dependency_diagram on
  `materi_strain_total_shear_kappa`; the interface variant
  (`group_interface_materi_plasti_mohr_coul_direct`) reaction magnitude
  (mohr_coul_direct3; mohr_coul_direct4 is a deliberate divergence, see
  above).

## Sprint 10 additions: compression_direct + pressure/coord limits

- `group_materi_plasti_compression_direct index sigy`: principal
  stresses lower than sigy are cut off (direct cut-off, no plastic
  strains; manual Professional 6.694). Spectral implementation:
  sigma = V diag(d) V^T, eigenvalues below sigy pulled up.
- `group_materi_plasti_compression_direct_visco index tm`: relaxation
  time; the cut relaxes with factor 1-exp(-dt/tm) (without the record
  the cut is total).
- `group_materi_plasti_pressure_limit index pressure_limit`: neglect the
  DIRECT plasticity laws when the pressure (positive in compression,
  -tr/3) exceeds the limit — free-surface problems (manual 6.688).
- `group_materi_plasti_coord_limit index coord_limit`: neglect them when
  the vertical coordinate exceeds coord_limit (manual 6.689).
- `group_materi_plasti_heat_generation index factor`: Professional name
  (with underscores) of the legacy `group_materi_plasti_heatgeneration`
  (fraction of plastic energy converted to heat; requires
  condif_temperature).

Tests: mdirect_comp (sigyy -10 capped at -5.0 EXACT), mdirect_gate
(pressure_limit disables the cap -> elastic -10; A/B fails without it).

## Sprint 10 lote 2: Drucker-Prager alias + analytic validation

- `group_materi_plasti_druck_prag index phi c phi_flow` — the
  Professional name (with underscore) of the legacy GNU
  `group_materi_plasti_druckprag` (alias in db_number; same physics,
  layout [phi c phi_flow], legacy tests druckpr1/examp23).
- `group_materi_plasti_bounda` / `_factor` — Professional names of
  `group_materi_plasti_boundary` / `_factor`.
- `group_materi_factor index factor` — multiplication factor for the
  material stresses AND stiffness (unit conversion).
- Registered partials: `group_materi_damping_method`,
  `group_materi_density_groundflow` (pre-existing),
  `group_materi_plasti_visco_exponential_limit/_name/_values`.

### Analytic validation (test mdp_shear)

The GNU implements f = sqrt(J2) + 3·alpha·sigma_m − K (the standard DP
form; NOTE sqrt(J2), not sqrt(3·J2)) with the manual's matching
constants. In PURE SHEAR with phi = 0 (alpha = 0, no pressure term, no
dilatancy) the capped stress is exactly

    sigma_xy = K = 6·c·cos(phi) / ( sqrt(3)·(3 − sin(phi)) ) = 2c/sqrt(3)

With c = 30: sigma_xy = 34.641. Measured at dt = 0.05: 34.636 —
0.014% error. (With phi ≠ 0 the associative flow generates normal
stresses in a displacement-imposed shear rig, so sigma_m ≠ 0; use
phi = 0 for the closed-form check.)

## Sprint 10 lote 3: Mohr-Coulomb CLASSIC (analytic validation)

`group_materi_plasti_mohr_coul index phi c phi_flow` — the classic
Mohr-Coulomb law (NEW implementation; the GNU only had the `_direct`
cut-off and the Professional manual advises `_direct` as "more stable
and fast" — both are now available). Surface (manual 6.724,
tension-positive, sig1 largest / sig3 smallest principal):

    f = 0.5(sig1 - sig3) + 0.5(sig1 + sig3) sin(phi) - c cos(phi) = 0

Flow rule with phi_flow (non-associative when different; phi_flow = 0
avoids dilatancy). Eigenvalues via matrix_eigenvalues; the flow
direction is provided by the standard tochnog driver (central finite
differences), so the law integrates with the existing cutting-plane
return mapping, plasti_kappa, boundary reduction, etc.

### Analytic validation (test mmc_tension)

In UNIAXIAL TENSION (sig1 = sigxx, sig3 = 0) the surface reduces to the
classic MC tensile strength:

    sig_t = 2 c cos(phi) / (1 + sin(phi))

With c = 40, phi = 30 deg: sig_t = 46.19. Measured: 46.19 EXACT (first
run). Sanity phi = 0 (Tresca): sig_t = 2c = 80; converges slowly
(79.55 at dt=0.02) because the Tresca vertex in uniaxial tension is a
singular point for the cutting plane (documented numerical behaviour,
not a law error).

NOTE for test designers: a displacement-imposed "pure shear" rig with
plastic flow does NOT stay in pure shear — the accumulated plastic
strain rotates the state to uniaxial tension at 45 degrees (sig =
tau*[[1,1],[1,1]]), which saturates at the SAME sig_t: measured sigxy
= c·cos(phi)/(1+sin(phi)) = 23.09 for the same parameters (that is how
the surface was first confirmed).
