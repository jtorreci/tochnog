# groundflow_velocity_sign (developer)

Sign convention of the groundflow (Darcy) velocity. This is a bug fix, not
a new keyword: the record `groundflow_velocity` already existed and the
element equation had the wrong sign.

## Files and functions

- `groundfl.cc`, `groundflow()` — the velocity equation assembled for the
  `-GROUNDFLOW_VELOCITY` dof:

  ```c
  // fixed (Professional sign)
  tmp = h[inol] * ( pe[jdim]*grad_new_unknowns[jdim*nuknwn+pres_indx] -
    old_unknowns[iuknwn] ) / dtime;
  element_rhside[indx] += volume * tmp;
  ```

  The previous code had `-pe[jdim]` (the 2013 F. Lorenzo change, see the
  header comment of the file). The file header comment about the 2013
  change is now superseded by the inline comment at the velocity equation.

- `general.cc`, `general()` — consumes the velocity for the coupled heat
  convection: for a `-CONDIF_TEMPERATURE` dof, `vel +=
  new_dof[gvel_indx + jdim*nder]` when `groundflow_velocity` is
  initialized (manual 8.3).

## Semantics

Professional manual 8.2 ("Groundflow velocities"): `v_i^g = k_ij dh/dx_j`.
The sign is the positive head gradient; the GNU sources had it negated.
The velocity is a solved dof (lumped inertia `1/dtime`), not a
post-processed quantity.

Measured reference (Professional 02-08-2026, `condif5.dat`):
node dofs `temp/gvelx/pres` = `1,10,-2` / `1,10,-1` / `0,10,0`.
GNU after the fix matches exactly; before it the velocity was `-10` and
the middle temperature `1.1e-16`.

## Verification

- `condif5` (sfnet corpus): before `temp = 1.11022e-16` (RUNFAIL), after
  `rc = 0` with `node_dof 2 temp = 1.0` = target, and the full node dof
  vector equals the Professional binary.
- Regression A/B on the sfnet 197 corpus (`scripts/final-corpus-run.sh`):
  `173 PASS / 24 RUNFAIL / 0 PARSE` -> `174 PASS / 23 RUNFAIL / 0 PARSE`;
  the per-test diff is ONLY `condif5 RUNFAIL -> PASS`, no regressions.
- Tests that also initialize `groundflow_velocity` (`ground1`, `ground3`)
  still pass: their targets are pressure-based, so the velocity sign does
  not change their classification.
- Internal suite (`scripts/build_safe.sh`): 16/16 runs OK.

## Notes / pending

- The Professional appears to compute `groundflow_velocity` only when a
  coupled process needs it (in a pure `-groundflow` model such as
  `ground1` the Professional reports `gvel = 0`, while GNU solves it
  always). This is a difference in the laziness of the computation, not
  in the sign; it does not affect any corpus target (the groundflow-only
  targets are pressures).
- The `groundflow_velocity` field is also used by `post_calcul` and the
  seepage/consolidation paths; the sign change only affects the direction
  of the field, and the A/B shows no regression on those tests.

## Differences with the Professional version

- **groundflow velocity output (`ground1`)** (deliberate divergence decided by
  the project owner on 2026-10-07; AUDIT K10; SEGUIMIENTO Descarte). GNU: v = k dh/dy =
  -0.05 (manual), at the nodes and at the integration points. Professional
  02-08-2026: 0 at the nodes and at the integration points (post-processing
  only; the pressure field agrees). Status: deliberate divergence (project
  decision 2026-10-07): the GNU follows the manual; the Professional
  behaviour is also Pro-suspect. Root cause measured 2026-10-07: the
  Professional zeroes the groundflow velocity where the calculated total
  pressure reaches `groundflow_total_pressure_limit`, whose default is 0
  (manual 6.587); with `groundflow_total_pressure_limit 1.e10` the
  Professional gives -0.05 = GNU (see `groundflow_total_pressure_limit.md`,
  DIVERGENCES TPL, pending owner decision). Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **`gvelx` sign in the interface-coupled groundflow path
  (`groundflow_interface`)** (known difference, not decided; BACKLOG C2).
  GNU: gvelx -1.99962e-6. Professional 02-08-2026: +1.99960e-6 (pressure
  field identical). Status: Open, LOW; root cause not established. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).
- **groundflow `pres` dof datum** (known difference, not decided; AUDIT
  K09). GNU: Stores the pressure (`ground4` pres 0). Professional
  02-08-2026: Stores the head h (500), as its manual says. Status:
  Convention difference (audit BOTH-plausible); total pore pressure and
  stresses identical. No decision. Index:
  [DIVERGENCES-PROFESSIONAL.md](../DIVERGENCES-PROFESSIONAL.md).

- **Decision 2026-10-07 (project owner):** align the GNU with manual 6.587 and the Professional on `groundflow_total_pressure_limit` (default 0, cap the calculated total pressure, zero the groundflow velocity where the limit is reached). This reverses the earlier K10 decision. Implementation scheduled for F4; see `ProjectDocs/DIVERGENCES-PROFESSIONAL.md` (TPL, K10).
