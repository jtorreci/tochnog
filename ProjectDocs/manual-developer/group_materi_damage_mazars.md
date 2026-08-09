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
- Note a subtle bug-style pattern at `damage.cc:92`:
  `epseq += epseq + epsp[idim]*epsp[idim]` — the left `epseq` is a stale
  (zero/previous) accumulator term. The professional result is still recovered
  in the regression test because the extra `epseq` terms are zero at the time of
  the first accumulation; it is a latent defect for multi-axial histories.
- `alpha` (brittleness) is clamped to [0,1]; `dt`, `dc` to [0,1].
- `matrix_jacobi()` (miscel.cc) computes the principal strains and `iwork[MDIM]`
  is a scratch eigenvector buffer.
- Regression coverage: `validation-suite/test-2014/damage1.dat`.

## External dependencies

- Core `db()`/`get_group_data()` accessors; `matrix_jacobi()`; the
  `materi_damage` unknown infrastructure in `materi.cc`.
- Globals `materi_damage`, `MDIM`, `new_sig`, `new_epe`.

## Hardcoded parameters / pending refactorings

- **Latent bug at `damage.cc:92`** (`epseq += epseq + ...`): should be
  `epseq += epsp[idim]*epsp[idim]`. Harmless in the current regression test but
  worth fixing in a dedicated work unit.
- The damage is stored as a nodal unknown dof (not element history); the
  element-level assembly in `materi.cc:491-494` uses the shape functions `h`
  to distribute it, which couples damage to the mesh interpolation.
- There is no user-facing switch to force a pure elastic run other than omitting
  the record.
