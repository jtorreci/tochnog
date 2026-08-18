# strain_settlement_parameters

## Implementación

- **Función**: `strain_settlement_creep()` in `materi.cc`. Called at the end
  of `set_deften_etc()`: adds the creep strain increment to `inc_ept` and
  updates `new_ept = old_ept + inc_ept`. Vertical component (y in 2D = index
  [4], z in 3D = index [8]); horizontal components scaled by `lateral_factor`.
- **Keywords** (data_class CONTROL, index 0) in `database.cc`:
  `strain_settlement_parameters` (DOUBLE 6), `_element_group` (INTEGER
  variable), `strain_settlement_diagram` (DOUBLE variable),
  `_diagram_dof` (INTEGER 1), `_diagram_number` (INTEGER 1). Enums in
  `tochnog.h`/`tochnog-mod.h` (between `STEP` and `STRESS`).
- **Ley de creep** (decisión 2026-08-18, el OCR del manual es ambiguo):
  `eps_zz(t) = Ar*(t/t_ref)^n/(t_plus+(t/t_ref)^n)`, saturada. `t` = tiempo
  desde la activación del elemento (`mesh_activate_gravity_time` o el global
  start). El incremento por paso es `eps(t_new)-eps(t_old)`.
- **Diagrama**: `strain_settlement_diagram` + `_dof` + `_number` hacen que un
  parámetro (1=time_plus, 2=Ar, 3=t_ref, 4=n, 5=lateral) dependa de un dof
  del primer nodo del elemento, vía `table_xy`.

## Validación

- `strain_settle`: quad4 con creep (Ar=0.01, t_plus=1, t_ref=1, n=1) → disy
  del nodo superior ≈ -0.000263 (el creep comprime la columna; sin el record
  disy=0).
- `strain_settle_diag`: `_diagram` hace Ar=0.02 (el doble) → disy ≈ -0.000526.
