# control_distribute (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — new enums
  `CONTROL_DISTRIBUTE_CORRELATION_DISTANCE/LENGTH`,
  `CONTROL_DISTRIBUTE_MINIMUM_MAXIMUM`, `CONTROL_DISTRIBUTE_PARAMETERS`,
  `CONTROL_DISTRIBUTE_SEED`, switch enum `LOGNORMAL`.
- `database.cc` — registrations (see control_distribute.md for the user
  view) + switch name `lognormal`.
- `distri.cc` — `distribute()` (called from `step_start`, top.cc): the
  record layout discriminates the path:
  - length 4 → Professional layout (new);
  - length multiple of 3 → GNU legacy triplets, code unchanged
    (`ho_othr4` regression).

## Professional path design

1. Seed: `distribute_ran = -(labs(seed))-1`; a negative idum makes the
   numerical-recipes ran1 reinitialize deterministically
   (`math.cc scalar_ran_uniform`), so equal seeds reproduce fields.
2. Parameters: mean/std OF THE VALUE. Lognormal uses the standard
   transformation to the underlying normal:
   `sigma_ln = sqrt(ln(1+(std/mean)^2))`, `mu_ln = ln(mean)-sigma_ln^2/2`,
   `draw = exp(mu_ln + sigma_ln*z)`.
3. Entities: group targets → elements filtered by `ELEMENT_GROUP`
   (specific index) or all elements whose group has the record (`-all`);
   node-like targets → the record indices (`-all` or one).
4. Spatial correlation: entity coordinates are element centroids
   (group path) or node coordinates (NODE-class items). Standard normals
   z_i are drawn independently, then for non-constant correlation
   `v_i = Σ_j w_ij z_j / sqrt(Σ_j w_ij²)` with
   `w_ij = exp(-d_eff_ij)`, `d_eff` the distance scaled per direction by
   the correlation lengths; pairs with physical distance beyond the
   correlation distance (default 4·L) get zero weight. The
   normalization keeps unit variance. Correlation length > 1.e12 →
   constant field (all entities share z_0).
5. Clamp to [min,max] on the drawn value.
6. Application:
   - group: `delta = draw - record_value` stored in
     `ELEMENT_DISTRIBUTE_VALUES`; the `+=` in `get_group_data`
     (group.cc) yields exactly the drawn effective value. The records
     are refreshed at every `distribute()` call (redraw per step).
   - node-like: the drawn value REPLACES `dval[number]` of the record
     (manual examples: nodal temperatures, y coordinates).

## Verification (suite 70/70)

- `cdist_normal`: two elements of one group, N(1000,50) seed 7 →
  per-element deltas read directly from `element_distribute_values`
  (generic targets): +35.4136 / -21.3726 — distinct draws, exact for the
  seed (reproducibility), group filtering by index exercised.
- `cdist_corr`: same model + `correlation_length 2e12` → both deltas
  IDENTICAL (35.4136) although the elements are 25 length units apart;
  with independent draws they differ (previous test).
- `cdist_clamp`: N(1000, 1e6) + `minimum_maximum 995 1005` → deltas
  exactly +5 / -5 (both draws cut to the bounds).
- `ho_othr4` (GNU layout) unchanged: 70/70 suite.

## Limitations

- Correlation smoothing is the exponential-kernel weighted average
  (O(n²) with distance cutoff), not a spectral/FFT method; adequate for
  the element/node counts of typical GNU models.
- distribute() fires per STEP of the owning control block: the field is
  redrawn at each step (matches GNU behaviour; use one-step blocks for a
  single draw).
