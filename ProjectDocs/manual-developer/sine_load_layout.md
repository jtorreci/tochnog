# sine_load_layout (`*_sine` family)

## Files and functions

- `area.cc` — within `area()`:
  - `FORCE_ELEMENT_EDGE_SINE` consumer: `type[itype]==FORCE_ELEMENT_EDGE`
    branch (~457-470).
  - `FORCE_ELEMENT_EDGE_NORMAL_SINE` consumer: `..._NORMAL` branch
    (~490-503).
  - `FORCE_ELEMENT_EDGE_PROJECTED_SINE` consumer: `..._PROJECTED` branch
    (~530-548).
  - `flux_edge_is_master` branch (~565-585): reads
    `flux_edge_companion(type[itype],3)` and therefore serves BOTH
    `GROUNDFLOW_FLUX_EDGE_NORMAL_SINE` (via `groundflow[3]`) and
    `CONDIF_HEAT_EDGE_NORMAL_SINE` (via `condif[3]`).
  - The computed factor is printed under `TOCHNOG_DEBUG=yes` by
    `pri("load", load)` (~610).
  - `time_end` local added to the `double` declaration block (~226).
- `condif.cc` — `CONDIF_HEAT_VOLUME_SINE` consumer in the heat-volume
  block (~130-148); `time_end` local (~68); debug print
  `pri("heat_volume load", load)` (~166).
- `post.cc` — `FORCE_ELEMENT_EDGE_SINE` mirror in
  `post_force_edge_summed_calculate()` (~652-664); `time_end` local
  (~585). No debug print for the factor.

## Implementation details

- Record layout (Professional): `[start_time, end_time, freq_0, amp_0, ...]`.
- Read: `time_start = data[0]; time_end = data[1];`
  `nfreq = ( db_len(...) - 2 ) / 2;`
- Evaluation: `frequency = data[2 + 2*i]; amplitude = data[3 + 2*i];`
  `load += amplitude * sin( 2*PIRAD * frequency * time_total );`
- Gating: `if ( time_total > time_start && time_total <= time_end )`.
  Outside the window `load` stays 0 (the consumer initialises `load=0`
  before the loop). `time_total = time_current + dtime` (end of step),
  consistent with the pre-existing start-time comparison.
- `force_element_volume_sine` (`force.cc`) is intentionally NOT changed:
  its Professional signature (6.507) has no `end_time`, so its
  `nfreq=(len-1)/2`, `data[1+2i]`/`data[2+2i]` layout is correct.

## Audit correction: `condif_heat_edge_normal_sine`

The bounded audit that produced the fix list marked this keyword as
"registered but no consumer" because a grep for the enum name only hits
`database.cc` and `check.cc`. The consumer exists but is reached through
the companion table: `flux_edge_is_master(type[itype])` for
`type[itype]==CONDIF_HEAT_EDGE_NORMAL` calls
`flux_edge_companion(CONDIF_HEAT_EDGE_NORMAL, 3)`, and
`flux_edge_companion()` returns `condif[3] == CONDIF_HEAT_EDGE_NORMAL_SINE`.
The same branch already handled `GROUNDFLOW_FLUX_EDGE_NORMAL_SINE`. The
layout fix at that shared branch therefore also fixes
`condif_heat_edge_normal_sine`; no new consumer was invented.

## Verification

GNU binary built clean (`scripts/build_safe.sh --clean`, suite 16/16).
Professional binary 02-08-2026, `data_ignore=-print_apply`, `LC_ALL=C`,
scratch cwd. Synthetic Pro-dialect inputs run in BOTH binaries; the GNU
`TOCHNOG_DEBUG=yes` factor print is compared with the Professional
response:

- `force_edge_sine`: GNU factor 5 (0 out-of-window), sigxx -25; Pro -25.
- `force_edge_normal_sine`: GNU factor 5 (0 out), sigxx 25; Pro 25.
- `force_edge_projected_sine`: GNU factor 5 (0 out); Pro response ratio
  5.0. The GNU projected-force APPLICATION NaNs on the single-element
  synthetic even without the sine (pre-existing projected-force issue),
  so only the prescribed factor is compared.
- `groundflow_flux_edge_normal_sine`: GNU factor 5 (0 out), pressure
  ratio 5.0; Pro ratio 5.0.
- `condif_heat_volume_sine`: GNU `heat_volume load` 5 (0 out),
  T=1.25; Pro T=1.25.
- `condif_heat_edge_normal_sine`: GNU factor 5 (0 out), temperature
  ratio 5.0; Pro ratio 5.0.
- post.cc mirror: GNU `post_force_edge_summed` (25, 0); Pro (25, 0).

Evidence: `/tmp/opencode/dyn8-sine-fix-evidence/`.

## Pending

- Frequency-0 static value (manual 6.34/6.465 special option) is not
  implemented; a zero frequency contributes 0.
- The projected-force application NaN on the single-element synthetic is
  a separate pre-existing issue (out of scope for the layout fix).

## External dependencies

None. `sin` from the C library, `PIRAD` the global pi constant.
