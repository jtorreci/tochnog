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

2. **Plane traction mode** (with `_normal` / `_normal_automatic`):
   `materi_direct_cutoff()` in `stress.cc` (the pre-existing plane cut-off:
   `max_fric = max(c - sig_n*tan(phi), 0)` on the plane with the given
   normal). NEW: `group_materi_plasti_bounda/_factor` (Professional 6.231/
   6.232) now reduce phi and c by the factor for the elements on the wall.

## Física

- Full mode (6.726): "Principal stress differences higher than allowed by
  the mohr-coulomb criterium are not allowed and will be cut off by
  Tochnog" — the alternative programming of the MC law. The cut does not
  use plastic strains; phi_flow enters the cut DIRECTION (the non-assoc
  flow ratio `(1+sin psi)/(1-sin psi)` of the max/min corrections; for
  psi=0 the cut is the mean-preserving difference cut).
- Plane mode (6.727/6.739): limits the friction/tension stress on the
  specific plane (interface semantics).

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
