# group_materi_plasti_cap1

## Implementación

- **Ley**: new block in `plasti_rule()` in `plasti.cc` (CRLF file),
  right after the legacy `GROUP_MATERI_PLASTI_CAP` block. It reads the
  group record `GROUP_MATERI_PLASTI_CAP1` (DOUBLE, fixed length 8) and
  implements both `GET_YIELD_RULE` and `GET_FLOW_RULE` (associative:
  `f_flow = f_yield`) for `plasti_type == -NONE` and for the explicit
  type. The flow gradient is provided by the standard driver via central
  finite differences — exactly like the other blocks.
- **Keywords** (data_class MATERI) registered in `database.cc`:
  - `group_materi_plasti_cap1` (DOUBLE, length 8, required GROUP_TYPE):
    `phi c M lambda_star kappa_star K_ref p_ref m`.
- **New enum**: `GROUP_MATERI_PLASTI_CAP1` in `tochnog.h` /
  `tochnog-mod.h` (kept in sync, same order; after
  `GROUP_MATERI_PLASTI_CAP`, before `GROUP_MATERI_PLASTI_COMPRESSION`).
- **History dof**: `pc` is read from `new_unknowns[cap1_indx]` (the
  `materi_plasti_cap1_history` initia); the plasti block exits with an
  error if the initia is missing.
- **check.cc**: `GROUP_MATERI_PLASTI_CAP1` requires `materi_stress`,
  `materi_strain_plasti` and `materi_plasti_cap1_history`.

## Física

The manual invariants (theory cap1):

```
p = -(sig11 + sig22 + sig33)/3                     (compression positive)
q = sqrt(0.5*[(s11-s22)^2 + (s22-s33)^2 + (s33-s11)^2] + 3*(s12^2+s23^2+s31^2))
f = q^2/M^2 + p*(p* - p*c),   p* = p + c*cot(phi),  p*c = pc + c*cot(phi)
```

In plasti.cc the invariants are computed from the raw tension-positive
`sig[]` array (MDIM*MDIM layout; the diagonal is `sig[0]`, `sig[4]`,
`sig[8]`, the shear components `sig[1]`, `sig[2]`, `sig[5]`) — the q
formula is invariant to the sign convention (only squared differences
and squares). `phi` is in **radians**; `M` is READ from the record (the
manual's `M = 6*sin(phi)/(3 - sin(phi))` is only a suggestion).

Singularities: `cot(phi) = cos(phi)/sin(phi)` diverges at `phi = 0`
(the manual expects `phi > 0`); `sin(phi) == 0` or `cos(phi) == 0`
raises `db_error`, as does `M == 0`.

## Cómo se actualiza pc (decisión de diseño)

The hardening is implemented with the **kappa pattern** (the same
machinery as `materi_plasti_kappa`), NOT inside the plasti block:

1. `set_stress()` (stress.cc) — after the plastic iterations converge
   (both the incremental branch and the cutting-plane branch, next to
   the `new_kappa` update), with the **converged `inc_epp`** of the
   step:

```
deps_p_cv = -(inc_epp[0] + inc_epp[4] + inc_epp[8]);   // compression positive
if ( deps_p_cv < 0. ) deps_p_cv = 0.;                  // no softening
new_cap1pc = old_cap1pc + deps_p_cv * K_ref/(lambda*/kappa* - 1)
                          * pow((old_cap1pc + c*cot(phi))/p_ref, m);
```

   The rate form of the manual is inverted and integrated explicitly
   with the OLD pc (the `(p*c/p_ref)^m` factor uses the step-start pc).
   The parameters are read from the element's `GROUP_MATERI_PLASTI_CAP1`
   record (`get_group_data ... GET_IF_EXISTS`); without the record
   `new_cap1pc = old_cap1pc` (the dof freezes). The element group guard
   is needed because the dof is global while the cap1 record is
   per-group.
2. `materi()` (materi.cc) assembles the evolution equation of the dof
   (same as kappa):
   `element_rhside += volume*h*(new_cap1pc - old_cap1pc)/dtime`; the
   solver resolves it with the lumped inertia from `general()` (dof type
   `MATERI_PLASTI_CAP1_HISTORY`, `inertia = 1`, `conv_part = 1`).
3. `dof.cc` clamps the dof to >= 0 after each step (kappa pattern).

Why this design: the plasti block is called 3 times per iteration
(yield, flow, flow-grad) and mutating `new_unknowns[cap1_indx]` there
would be order-dependent and inconsistent with the RHS assembly; the
kappa pattern computes the target value once per step from the
converged plastic strain and lets the RHS drive the dof to it — the
same weak-coupling behaviour as kappa (the yield surface within a step
uses the previous iterate of pc; with `control_timestep_iterations` the
coupled system converges, see the discrete fixed point below).

**Discrete fixed point**: per plastic step the cutting-plane returns the
stress to the CURRENT pc dof, so the hardening increment satisfies
`dpc = (deps_p_el - dpc)/(lambda*/kappa* - 1)` (with the elastic
pressure jump `deps_p_el = K*deps_vol`); for `lambda*/kappa* = 10` and
`deps_p_el = 5.0`: `dpc = 0.5/step` (continuum bound `0.5556/step`).
Measured `pc = 110.0` exactly after 20 plastic steps.

## Validation

- `mcap1.dat`: hex8 isotropic compression (`phi = 0.5`, `c = 10`,
  `M = 1.2`, `lambda* = 0.2`, `kappa* = 0.02`, `K_ref = 833.3`,
  `p_ref = 100`, `m = 0`, `pc_0 = 100`, 40 load + 5 unload steps of
  0.05): cap activates at `p = pc = 100`, `pc` grows 0.5/step to
  `109.9996`, elastic unload ends at `sigxx = -84.9996`
  (analytic `-85.0`). Targets `-sigxx -85.0 ± 4`, `-pc 110.4 ± 2.5`.
- `mcap1_elast.dat`: elastic twin -> `sigxx = -175.0000` (exact).
- `mcap1_comb.dat`: cap1 + druck_prag -> identical to mcap1 (max-f
  selection; `f_dp < 0` on the isotropic path).

## Pendiente / gotchas

- GOTCHA (test, not code): `bounda_time_increment` is a LOCAL in
  `bounda()` that persists across bounda records (`db` with
  GET_IF_EXISTS does not reset it). A pairs record AFTER an increment
  record is misinterpreted as load-only (its time values become loads —
  observed `velx = +100` on the fixed face at step 2). Fix in the tests:
  use the increment format on ALL bounda records (zeros on the fixed
  faces). msph works only because its increment record is the LAST one.
- GOTCHA (test): `k = floor(t/increment)`, so the j-th value applies at
  `t = (j-1)*dt`; 41 load + 5 unload values give 40 loading steps + 5
  unloading steps.
- The hardening uses the TOTAL plastic volumetric strain of the step
  (trace of inc_epp). When combined with a shear law with dilatancy,
  the shear volumetric strain also contributes to pc. A cleaner (but
  more invasive) design would gate the update on
  `plasti_type == GROUP_MATERI_PLASTI_CAP1` at convergence; the current
  version matches the manual for cap-dominated compression and is
  documented.
- CRLF: `plasti.cc` is CRLF; the block was added preserving CRLF.
