# group_materi_damage_mazars

## Files and functions

- `damage.cc`:
  - `damage()` (line 24) — dispatcher: reads `GROUP_MATERI_DAMAGE_MAZARS` for
    the element group and calls `damage_mazars()`; if the record is absent the
    damage is kept at its old value (`new_damage = old_damage`).
  - `damage_mazars()` (line 45) — the Mazars law itself:
    1. Parameters: `eps0, at, bt, ac, bc, beta` (`materi_damage_mazars[0..5]`).
    2. If `|epe| < 1e-10`: `d = 0`.
    3. Else: computes `gamma` from the compressive stress components
       (`sum(sigma_-^2)` / `sum(sigma_-)`), the principal elastic strains via
       `matrix_jacobi()`, then `epseq`, `alpha`, `dt`, `dc` and the combined
       `new_damage = dt*alpha^beta + dc*(1-alpha)^beta`.
    4. Monotonicity clamp: `if (new_damage<old_damage) new_damage = old_damage`.
- Entry point called from `materi.cc` in the `materi_stress` block:
  - `old_damage`/`new_damage` are read from the `materi_damage` unknown dofs
    (lines 212-215) and passed into the constitutive update (line 302);
  - the stress is scaled by `(1-d)` and the damage residual is assembled into
    the nodal force vector (lines 491-494, `tmp = volume*h*(new-old)/dtime`).
- `database.cc:2245` — registration: `GROUP_MATERI_DAMAGE_MAZARS`,
  `DOUBLE_PRECISION`, length 6.
- `database.cc:3139` — the `materi_damage` unknown; `database.cc:4194` maps
  `-MATERI_DAMAGE` dof labels.
- Enum mirrored in `tochnog.h` and `tochnog-mod.h`.

## Implementation details

- `epseq = gamma * sqrt(sum(pos. princ. strains^2))` is the driving variable;
  `gamma` (compressive triaxiality) is computed from the CURRENT stress tensor
  `new_sig` (compressive components only), so the model distinguishes tension
  and compression without an explicit crack-orientation split.
- **Bug fixed (2026-08-10, commit [P4-A])**: the equivalent-strain accumulation
  at `damage.cc:92` used `epseq += epseq + epsp[idim]*epsp[idim]`, which is
  `epseq = 2*epseq + epsp^2` — a DOUBLING accumulation instead of a plain sum
  of squares. It was present verbatim in the original 2014 GNU sources
  (`Sources-2014.zip`), so it affected every GNU version. The correct line is
  `epseq += epsp[idim]*epsp[idim]`. See the "Bug report and fix" section below.
- `alpha` (brittleness) is clamped to [0,1]; `dt`, `dc` to [0,1].
- `matrix_jacobi()` (miscel.cc) computes the principal strains and `iwork[MDIM]`
  is a scratch eigenvector buffer.
- Regression coverage: `validation-suite/test-2014/damage1.dat` (uniaxial) and
  `damage2.dat` (non-equal biaxial, added with the fix).

## Bug report and fix (2026-08-10)

**Symptom.** In multi-axial tension states (two or three positive principal
strains) the Mazars damage grew faster than the model defines, activating
earlier than it should.

**Root cause.** `damage.cc:92`:

```c
if ( epsp[idim]>0. ) epseq += epseq + epsp[idim]*epsp[idim];  // BUG
```

`epseq += epseq + x` is `epseq = 2*epseq + x`. The intended Mazars equivalent
strain is the square root of the sum of squares of the positive principal
strains, so the accumulation must be a plain sum of squares:

```c
if ( epsp[idim]>0. ) epseq += epsp[idim]*epsp[idim];          // FIX
```

The line directly above (`epssiz += epsp[idim]*epsp[idim]`) already does the
correct sum of squares for the norm, which is the strongest evidence of the
intended form.

**Numerical trace** (reproduced in `/tmp/mazars_trace.c`):

| State | bug (2*epseq+x) | correct (sum of squares) | over-estimation |
|-------|-----------------|--------------------------|-----------------|
| uniaxial `{2e-4, 0, 0}` | 0.0002 | 0.0002 | 0% |
| biaxial `{2e-4, 1e-4, 0}` | 0.0003 | 0.0002236 | +34% |
| triaxial `{2e-4, 1e-4, 5e-5}` | 0.000427 | 0.000229 | +86% |

**Why the regression tests did not detect it.** `damage1.dat` uses a single
`-bar2` element (uniaxial extension), where only ONE principal strain is
non-zero: `epseq = 2*0 + epsp^2 = epsp^2`, identical to the correct form. The
bug only manifests when at least TWO principal strains are positive
simultaneously (biaxial/triaxial tension).

**Verification.** New regression `validation-suite/test-2014/damage2.dat`:
single `-quad4` element under non-equal biaxial extension
(`velx=0.5`, `vely=0.1`), same Mazars parameters as `damage1.dat`.
Target `dam = 0.471973` (tolerance 1e-5):
- fixed binary: passes (exit 0);
- buggy binary: produces `dam = 0.789301` and FAILS the target.

`damage1.dat` and all P4 family tests still pass unchanged with the fix.

## Effective stress and gamma (2026-09-29)

- `damage()` is now called with the EFFECTIVE (undamaged) stress; the
  (1-d) compounding defect of the stress update is fixed, see
  `materi_damage_effective_stress.md`.
- `gamma` reads `new_sig[0..2]` = xx, xy, xz (not principal stresses). It
  is not objective, but the Professional 02-08-2026 is not either and
  matches no candidate form (GNU form, principal-stress gamma, gamma = 1);
  probes and values in `materi_damage_effective_stress.md`. Left unchanged.

## External dependencies

- Core `db()`/`get_group_data()` accessors; `matrix_jacobi()`; the
  `materi_damage` unknown infrastructure in `materi.cc`.
- Globals `materi_damage`, `MDIM`, `new_sig`, `new_epe`.

## Hardcoded parameters / pending refactorings

- The damage is an integration-point state since C68 (F4 T6f,
  2026-10-12): see "Integration-point state (C68)" below. Before, the
  constitutive update read the h-lumped nodal value back, which coupled
  the damage to the mesh interpolation.
- There is no user-facing switch to force a pure elastic run other than omitting
  the record.
- The `epseq` accumulation bug (see "Bug report and fix") is FIXED. The
  upstream 2014 GNU sources still carry the buggy line; if the sources are ever
  re-imported the fix must be re-applied.

## Differences with the Professional version

- **Mazars equivalent strain** (deliberate divergence, class A, project
  decision 2026-10-12; BACKLOG C7). GNU (`damage_mazars()`, `damage.cc`): manual 2.2.5, eps_eq =
  sqrt(sum <eps_i>+^2) of the positive principal strains (hex8
  diag(2e-4,-3e-4,1e-4): dam 0.28610207515). Professional 02-08-2026:
  accumulates them in ascending order as s = 2 s + eps^2 (weights 4, 2, 1;
  the 2014 accumulation defect the GNU fixed on 2026-08-10): dam 0.36409.
  It differs whenever two or three principal strains are positive, also
  through a von Mises + damage path (bar with `vonmises 3e6`: GNU dam
  0.925372, eppxx 1.87863e-4; Professional 0.96946, 1.368e-4). Guard
  tests `test-2026/damage_mazars_eq_two_positive1`,
  `damage_mazars_vonmises1`. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- The triaxiality factor gamma is built from the xx, xy, xz stress
  components in both codes (`damage.cc:84-95`) (not objective: a compressive xx with xy or xz
  shear gives an orientation-dependent `dam`). Same in the Professional,
  so not a divergence; noted in BACKLOG C7, no fix planned.

## Integration-point state (C68, F4 T6f, 2026-10-12)

`materi_damage` is an integration-point history state (Professional manual
4.10), like the hardening parameter kappa (`materi_plasti_kappa.md`).

- `elem.cc` `elem()`: the ELEMENT_DOF read and write masks include
  `(materi_damage && i==dam_indx)`, so `old_damage` in `materi.cc` is the value remembered at the
  integration point, not the h-lumped nodal value interpolated back.
- `materi.cc` `materi()`: one block right after the kappa block (after
  the `stress_recovery_post_pass` skip, outside the `|h| > EPS_H`
  guard) assembles the nodal (output) numerator
  `volume * w * (new - old) / dtime` with `w = strain_recovery_weight()`
  and writes the IP value to `new_unknowns[<indx> + k*nder]` when
  `options_element_dof -yes`. The old h-weighted block was removed.
- `general.cc` `general()`: the lumped denominator of `MATERI_DAMAGE` uses the
  same weight.

Verification: `test-2026/damage_ip_state1` (four tria3, 1-point rule, non-uniform
field, targets = Professional 02-08-2026): previous binary rc 1 (node
reactions 1.9e-1 relative off), fix rc 0 (5e-13); nodal dam within
2e-3.
