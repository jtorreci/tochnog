# post_node_rhside_ratio_method

## Implementation (F5 T2 V9, C110)

- **Record:** INTEGER, length 1, no index; only `-post_node_rhside_free` is accepted (`db_error` otherwise).
- **Consumer:** `post.cc`, `post_node_rhside_fixed_free()`. One node loop over the active `NODE_RHSIDE` records accumulates, per principal dof (`dof_principal>=0`) and per `NODE_BOUNDED` state, the signed sum (stored as `POST_NODE_RHSIDE_FIXED` / `_FREE`), the sum of absolute values and the node count. The ratio is computed per dof from the means: `mean_free/mean_fixed`, or `mean_free` when `mean_fixed<EPS_RHSIDE_FIXED` (1.e-10, manual 6.950) or the method is `-post_node_rhside_free`; the largest dof value is stored. `POST_NODE_RHSIDE_RATIO_UNKNOWNTYPES` filters the dofs by `dof_type`.
- **Callers:** `top.cc`, after every equilibrium solve (the ratio drives `control_timestep_iterations_automatic` and `control_timestep_size_automatic_decrease`), and once more after the post-solve material pass (C110), so the stored records are those of the final residual. The local `post_node_rhside_ratio` of `top.cc` keeps the last in-loop value.
- **Hardcoded:** `EPS_RHSIDE_FIXED` 1.e-10 (manual 6.950; was 1.e-6 with no fallback before C110).

## Verification

Professional 02-08-2026 probes (v9t t3 input, scratch `c110pro`): `post_node_rhside_ratio1` (default ratio 12.9442, signed sums) and `post_node_rhside_ratio_method1` (1.23884), checked in `scripts/build_safe.sh` with a relative tolerance of 2e-4. Pro `other/slide4`: fixed -5.770607e-3 / 9.995002e-3 = Professional (C63).

## Pending

- The manual default `post_node_rhside_ratio_dof_type -materi_velocity` (6.951) when the record is absent is not applied: without it every principal dof enters the ratio (not measured on the Professional).
