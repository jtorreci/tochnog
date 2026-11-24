# dof_limit (developer)

- Enum `DOF_LIMIT` (`tochnog.h` / `tochnog-mod.h`, before `LAST_DUMMY`),
  registered in `database.cc` next to `DOF_TYPE`: DOUBLE, length
  `2*npuknwn`, `no_index`, class DOF, fixed length (the generic parser
  error "Problem reading : dof_limit ... Number of data values expected"
  for another count, the Professional text).
- `dof.cc`: `static dof_limit_clamp( node_dof_new, limit )` clamps slot
  `ipuknwn*nder` of every primary dof to `[limit[2i], limit[2i+1]]`.
  - `parallel_new_dof_diagonal()` reads the record once
    (`GET_IF_EXISTS`, buffer `2*MPUKNWN`) and clamps each node twice: after
    the solved / non-principal update loop (so the displacement, velocity
    integrated, acceleration and derivative updates use the clamped
    velocity) and at the end of the node (the integrated dofs).
  - `dof_limit_apply()` (serial, every active node, VERSION_NEW) is called
    from `top.cc` at the end of the equilibrium block of each time step,
    after the hoisted stress sync (which assigns stresses outside
    `parallel_new_dof_diagonal`) and after `node_dof_empty_apply()`.
- Prescribed nodes are clamped too; `bounda()` re-imposes their values at
  the next iteration start.

Verification: Professional 02-08-2026 probes f5/v5b `b_lim` (temperature,
free and prescribed nodes clamped, free node unaffected), `b_lim_hi`
(upper only), `b_lim_short` / `b_m4` / `b_m6` / `b_m24` (wrong counts
rejected: one pair per `dof_label` entry, not per principal dof),
`b_m22` (velx clamped, disx follows, element_dof / element_intpnt_dof not
clamped), `f_s1` (averaged stress clamped in node_dof only). Tests
`dof_limit1`, `dof_limit2`, `dof_limit_error1`.

Hardcoded: none. Pending: the Professional leaves the prescribed values
unclamped in its element pass; with `derivatives` one pair per primary dof
(assumption, not compared).
