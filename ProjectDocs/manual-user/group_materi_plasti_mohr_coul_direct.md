# group_materi_plasti_mohr_coul_direct / tension_direct

## Description

`group_materi_plasti_mohr_coul_direct` and `group_materi_plasti_tension_direct`
are the **direct stress cut-off** variants of the Mohr-Coulomb and tension
plasticity laws (manual Professional 6.726 / 6.738). They are direct cut-offs
of the stresses (no plastic strains, no incremental return mapping), which
makes them very stable.

The record has TWO modes, selected by the presence of the plane-normal
records:

### 1. Full principal-stress cut-off (no `_normal` / `_normal_automatic`)

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

### 2. Plane traction cut-off (with `_normal` / `_normal_automatic`)

With a plane normal the records limit the traction on that SPECIFIC plane
(manual 6.727/6.739, the "i.c.w." interface semantics):

- `mohr_coul_direct phi c phi_flow` + `_normal nx ny nz`: the shear traction
  on the plane is capped at `max_fric = max(c - sig_n*tan(phi), 0)` where
  `sig_n` is the normal traction (compression increases the limit).
- `tension_direct sigy` + `_normal nx ny nz`: the normal traction on the
  plane is capped at sigy.

`_normal_automatic -yes` takes the plane normal from the element normal.

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

The normal record accepts ndim values (1 value in 1D, 2 in 2D, 3 in 3D).

## Parámetros

| Record | Parameters | Meaning |
|--------|------------|---------|
| `group_materi_plasti_mohr_coul_direct` | `phi c phi_flow` | Friction angle (rad), cohesion, dilatancy angle (rad). |
| `group_materi_plasti_mohr_coul_direct_normal` | `nx ny nz` | Plane normal (explicit; ndim values). Selects the plane mode. |
| `group_materi_plasti_mohr_coul_direct_normal_automatic` | `-yes` | Take the normal from the element normal. |
| `group_materi_plasti_mohr_coul_direct_visco` | `tm` | Visco relaxation time. |
| `group_materi_plasti_mohr_coul_direct_wall` | `phi c phi_flow` | Values used when the element is attached to a wall. |
| `group_materi_plasti_tension_direct` | `sigy` | Tensile limit on the principal stresses (full mode) or on the plane normal traction (plane mode). |
| `group_materi_plasti_tension_direct_normal` | `nx ny nz` | Plane normal (explicit). |
| `group_materi_plasti_tension_direct_normal_automatic` | `-yes` | Take the normal from the element normal. |
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

## Estado de implementación

- **Implementado**: the full principal-stress mode (spectral tension cap +
  MC difference cut with the non-associative flow) and the plane mode
  (`_normal`, `_normal_automatic`, `_visco`, `_wall`).
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
  (mohr_coul_direct3/4).

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
