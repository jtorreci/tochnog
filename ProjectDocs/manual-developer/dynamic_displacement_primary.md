# Displacement-primary dynamics for `materi_dynamic < 1`

## Checkpoint: 2026-09-08

The gated path in `so.cc` is restricted to the plain elastic mechanical
transient (`materi_velocity` plus `materi_displacement`, `inertia_apply -yes`,
Bi-CG, no plasticity, damping, viscosity, MPC, or multi-physics records).
For its `f=0` normal pass it holds the truss and continuum operator in the
reference configuration. During the auxiliary `f=1` reassembly,
`disp_primary_reference` preserves that selection.

The implementation uses the measured Professional map
`a = 2*w1 - w2`, where each `wi` is solved with the same
`M + dt^2 K` operator and the second right hand side is `M*w1`. It closes
`dynamic1` and `dynamic2`: both exit rc=0 at
`disx=0.0958890141894` after 100 s.

It is not yet a complete continuum closure. `dynamic5` remains rc=1 at
`-8.500245769e-05` (Professional: `-8.587735696826e-05`) and `dynamic8`
remains rc=1 at `-3.38705737383e-04` (Professional:
`-9.985253775e-04`). In dynamic8 the second Bi-CG solve repeatedly breaks
down and the code falls back to `x2=0`; investigate that auxiliary continuum
system before broadening the gate.

## Factor resolution (F4 T6d / C11, 2026-10-11)

The gate and `disp_primary_f` take the factor from `general.cc`
`materi_dynamic_effective()`: `control_materi_dynamic` at the current
ICONTROL overrules `materi_dynamic` (manual 6.141). Before, the gate read
only the global record and refused any input with
`control_materi_dynamic` (dynamic1 with only `control_materi_dynamic 10
0.`: 6.4477, Professional 0.0958890; now equal). The auxiliary `f=1`
reassembly no longer PUTs `MATERI_DYNAMIC = 1` and restores it: the
helper returns 1 while `disp_primary_reference == -YES`.
