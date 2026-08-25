# control_materi_plasti_hardsoil_gammap_initial

## Implementación

- **Record**: `control_materi_plasti_hardsoil_gammap_initial` (INTEGER,
  length 1, data_class CONTROL) — the enum and the registration already
  existed (Sprint 9 stub); the behavior is new.
- **Read**: in the elastic HS block of `set_stress()` (`stress.cc`),
  `db( ..., icontrol, &swit, ..., GET_IF_EXISTS )` with the current
  ICONTROL (the index selects the matching `control_timestep` block,
  manual 6.146).
- **First-timestep initialization** (`stress.cc`, in the elastic HS
  block): if the switch is `-yes` and the element has the HS plastic
  group, the record
  `ELEMENT_INTPNT_MATERI_PLASTI_HARDSOIL_GAMMAP_INITIAL` is checked
  (GET_IF_EXISTS): a value `< 0` means "not initialized yet" (the
  sentinel). On the first timestep the block computes
  `gamma_p_extra = q/(E50*(1-q/qa)) - 2*q/Eur` at the INITIAL stress
  state (`new_sig` at the start of the first step) and PUTs it. If the
  initial `q >= qa` (beyond the asymptote) it clamps to `1.e10` (no
  finite gamma_p brings f to zero, documented); if `gamma_p_extra < 0`
  (the initial state is already inside) it clamps to 0.
- **Pre-allocation**: `top.cc` (before the time loop, pattern of
  `NONLOCAL_ELEMENT_INFO`) allocates the element record for every
  element with the sentinel `-1.` when the control record exists — the
  first-timestep PUT runs inside the PARALLEL element loop, where db
  allocation is forbidden ("Data is allocated in a parallel loop").
- **Consumption**: the hardsoil plastic block in `plasti_rule()`
  (`plasti.cc`) reads the record (`GET_IF_EXISTS`, buffer initialized to
  0, negative values treated as 0) and adds it to gamma_p in the yield
  function from then on.
- **Record**: `element_intpnt_materi_plasti_hardsoil_gammap_initial`
  (DOUBLE, length 1, data_class ELEMENT) registered in `database.cc`
  (the manual stores one value per integration point; this port stores
  ONE value per element — the last integration point's value — see
  PENDING).

## Física

Manual (6.146): "add an initial contribution to gamma_p exactly such
that the yield function is zero-valued. This is convenient to start the
calculation with hardsoil with deviatoric stresses which would have been
outside the yield surface without this extra contribution." With the
extra, `f(initial) = 0` and no initial plastic return occurs; the extra
stays in the record and is added to gamma_p forever (the surface starts
at `gamma_p_extra`).

## PENDING / simplifications (documented)

- **Per-integration-point storage**: the manual saves
  `gammap_initial_integration_point_0, _1, ...` per integration point.
  This port stores ONE value per element (the last integration point's
  value): for uniform stress states (all tests) it is exact; for
  non-uniform initial stresses the per-integration-point array is
  pending (would need the ipoint index in `set_stress` and a
  `data_length = npoint` record).
- The initialization is keyed on the element record sentinel, so it runs
  in the FIRST timestep of the FIRST control block that processes the
  element; with several `control_timestep` records the manual's
  per-index behavior is approximated (documented).
- `gamma_p_extra` is evaluated with the FIRST-LOADING moduli (E50/Eur at
  the initial state, consistent with the yield function which uses both).

## Validation

- `mhardsoil_gp0`: `gamma_p_extra = 2/(1000*(1-2/38.49)) - 4/3000 =
  0.0007763` EXACT in the element record (target on the record
  directly); the initial stress stays `sigxx = -2.0` EXACT and `kappa`
  stays 0 (f = 0 at the start, no return).
- `mhardsoil_gp0_off` (A/B without the control): the same state is
  outside the surface, the return relaxes the deviatoric stress
  (`sigxx -> -1.225`) and `kappa` grows to `0.000293`.
