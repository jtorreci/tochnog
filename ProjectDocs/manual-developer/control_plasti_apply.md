# control_plasti_apply — developer

Per-timestep form of `plasti_apply`: enum `CONTROL_PLASTI_APPLY` (tochnog.h /
tochnog-mod.h, before `PLASTI_APPLY`), record INTEGER 1 class CONTROL
(indexed) in `database.cc`, consumed through
`materi_apply_gate_off( PLASTI_APPLY, CONTROL_PLASTI_APPLY )` in
`set_stress()` (stress.cc) and `slide_spring()` (slide.cc). The record at the
current `ICONTROL` overrules `plasti_apply`. Implementation details, the
Professional measurements and the pending items are in `plasti_apply.md`.
