# groundflow_pressure_atmospheric (developer)

Upper cap on the static/total pore pressure computed from the phreatic
level. Keyword inherited from the GNU 2014 sources (not in Professional);
this change only documents and verifies it — no code was modified.

## Files and functions

- `database.cc:2476-2480` — keyword registration: `DOUBLE_PRECISION`,
  `data_length = 1`, `no_index = 1` (single global value), 
  `data_class = GROUNDFLOW`.
- `tochnog.h:547` / `tochnog-mod.h:521` — enum
  `GROUNDFLOW_PRESSURE_ATMOSPHERIC` (the two headers must stay in sync).
- `groundfl.cc` — `groundflow_phreatic_coord()`:
  - Reads the keyword with `GET_IF_EXISTS` (line 226); local default 0.
  - Clamp at lines 323-325, applied to the outgoing values:
    ```c
    if ( static_pressure>=pressure_atmospheric ) static_pressure = pressure_atmospheric;
    if ( total_pressure>=pressure_atmospheric ) total_pressure = pressure_atmospheric;
    total_pressure += addtopressure;
    ```
    (`groundflow_addtopressure` is added AFTER the clamp.)
- `check.cc:345` — requires the `groundflow_pressure` unknown.

## Semantics (sign convention)

With `force_gravity (0,-1)`, static pressure `fg[ndim-1]*dens*(wl-y)` is
NEGATIVE below the phreatic level (compression) and POSITIVE above it
(suction). The clamp `value >= pa -> pa` therefore caps the SUCTION side:

- default `pa = 0`: suction above the phreatic surface is annulled.
- `pa > 0`: suction is kept up to `pa`.
- compression (negative values) is never clamped for `pa >= 0`.

## Consumers of the clamped values

- `groundfl.cc:410-416` — `groundflow_phreatic_level_multiple_static -yes`
  imposes `node_dof[pres] = static_pressure`, i.e. the CLAMPED value.
- `groundfl.cc:354-367` — `groundflow_phreaticlevel_bounda` METHOD1/2 set
  `dens*fg*location` instead: NOT clamped (uses the water level position).
- `calcul.cc:603-630` — post-calculations that print/derive total pressure.
- `area.cc:358,537` — phreatic-dependent edge load computations.
- `materi.cc:418` — water unit weight contribution in the equilibrium.
- `bounda_water` (`bounda.cc`) bypasses `groundflow_phreatic_coord()`
  entirely and computes the direct hydrostatic formula, so the cap does
  NOT apply there (see `bounda_water.md`).

## Verification

Tests `groundflow_pressure_atm` and `groundflow_pressure_atm_def`
(validation-suite/test-2014, registered in scripts/build_safe.sh):

- `groundflow_pressure_atm`: `groundflow_pressure_atmospheric 0.5`, level
  y=1, `_static -yes`. Targets: node y=0 `pres = -1.0` (compression passes),
  node y=2 `pres = +0.5` (suction +1.0 capped to 0.5). PASS.
- `groundflow_pressure_atm_def`: no keyword (default 0). Targets: node y=0
  `-1.0`, node y=2 `0.0` (suction annulled). PASS.
- A/B: the two inputs differ ONLY in the keyword; the suction node differs
  `+0.5` vs `0.0`, proving the causal effect of the clamp.
- Full suite 57/57 OK after the change (build_safe.sh).

## Pending refactorings

- None for this keyword (code inherited). Open question recorded in
  `bounda_water.md`: should `bounda_water` apply the same cap for
  consistency with the unsaturated model?
