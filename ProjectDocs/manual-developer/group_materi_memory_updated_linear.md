# group_materi_memory -updated_linear

## Implementation

The `-UPDATED_LINEAR` memory type is handled in three places:

1. **materi.cc, `materi()`**: the UPDATED branch (rotated old stress
   from old_unknowns) now includes `memory==-UPDATED_LINEAR`:
   ```c
   else if ( memory==-UPDATED || memory==-TOTAL_LINEAR ||
       memory==-UPDATED_WITHOUT_ROTATION || memory==-UPDATED_LINEAR ) {
     if ( materi_stress )
       array_move( old_sig, rotated_old_sig, MDIM*MDIM );
     ...
   }
   ```
   and the stress-rotation block treats `-UPDATED_LINEAR` like
   `-UPDATED` (`inc_rot` applied).

2. **materi.cc, `set_deften_etc()`**: `-UPDATED_LINEAR` gets the
   identity rotation matrices and the LINEAR engineering strains
   (same as `-UPDATED_WITHOUT_ROTATION`):
   ```c
   if ( memory==-UPDATED_WITHOUT_ROTATION || memory==-UPDATED_LINEAR ||
        memory==-TOTAL_LINEAR ) { ... identity rotations ... }
   ...
   if ( memory==-UPDATED_WITHOUT_ROTATION || memory==-UPDATED_LINEAR ) {
     // linear engineering strains
   ```
   The velocity-gradient term that builds `inc_deften` (and hence
   `inc_ept`) is REFERENCE-geometry based for `-UPDATED_LINEAR` only;
   see section 4.

3. **stress.cc, `set_stress()`**: the compressibility block treats
   `-UPDATED_LINEAR` like `-TOTAL_LINEAR` (linear strain measures):
   ```c
   if ( memory==-TOTAL_LINEAR || memory==-UPDATED_LINEAR ) { ... }
   ```

Additionally the UPDATED memory self-check (stress.cc) only rejects
`materi_displacement` when the updated formulation is EXPLICITLY
requested (`db_active_index(GROUP_MATERI_MEMORY, gr, ...)`): with the
default memory (record absent) the Professional accepts
velocity+displacement together (the Masin clay corpus tests hypo12/13
run without a memory record).

4. **Reference-geometry strain increments (R2, 2026-09-13)**. The
   `-UPDATED_LINEAR` total must be the LINEAR engineering strain referred
   to the UNDEFORMED mesh, not to the current configuration. With the
   default velocity-integrated `-FOLLOW_MATERIAL` mesh (`top.cc:105-109`)
   the nodes are advanced by `vel*dtime` in `locate()` (`locate.cc:24-54`,
   called from `top.cc:550`), so `NODE VERSION_NORMAL` read in `elem()`
   (`elem.cc:452`) is the CURRENT geometry; `pol()` then builds `old_d`
   from it (`polynom.cc:655-661`) and `inc_deften = I + grad_current(v)*dt`
   (`materi.cc:1472-1479`). Accumulating that gives the LOGARITHMIC strain
   (`eptyy = -ln(0.7) = -0.356567822765` for hypo7), not the engineering
   `-0.3` the Professional reports.

   Fix: for `-UPDATED_LINEAR` only, `elem()` additionally loads the
   reference coordinates from `NODE_START_REFINED` (populated in
   `top.cc:272-274` before the time loop and never advanced by `locate()`),
   and `pol()` produces an extra gradient `ref_d` from them (a separate
   pass, so the current-geometry `detj`/`volume` is untouched). `elem()`
   forms `old_grad_ref_new_unknowns = ref_d @ new_dof` and threads it
   through `materi()` -> `set_deften_etc()`. There, when
   `memory==-UPDATED_LINEAR && old_grad_ref_new_unknowns!=NULL`, the
   velocity gradient uses the reference derivatives; every other memory
   option keeps `old_grad_new_unknowns` (the current geometry). The
   axisymmetric hoop increment (`inc_deften[8]`) uses the reference radius
   `R_ref` in the same gated way, so the total is `sym(grad_ref(u))` (the
   engineering tensor). `old_deften`/`new_deften` are unaffected.

   Signature changes (both `tochnog.h` and `tochnog-mod.h` must stay in
   sync): `pol()` gains `double ref_coord[]` (in) and `double ref_d[]`
   (out); `materi()` gains `double ref_coord[]` and
   `double old_grad_ref_new_unknowns[]`; `set_deften_etc()` gains
   `double ref_coord[]` and `double old_grad_ref_new_unknowns[]`.
   `elem()` passes `NULL` for all of them when the memory option is not
   `-UPDATED_LINEAR`, so the reference path is fully inert otherwise
   (verified: hypo1-4 `.dbs` byte-identical).

## Verification

Measured 2026-09-13 on `dev/masin-linear-strain` (R2) against a clean
`a02e7c2` baseline build and the Professional binary
(`tochnog_version_02-08-2026`). `eptyy`/`sigyy` from the last
`post_point_dof` record; `hyhis0` is the void-ratio history variable.

| test | route | eptyy | sigyy | hyhis0 | Pro target sigyy |
|------|-------|-------|-------|--------|------------------|
| hypo7 | a02e7c2 baseline | -0.356567823 | -230.793 | 1.110389 | -231.80 |
| hypo7 | R2 | **-0.300000000** | -235.580 | 1.099038 | -231.80 |
| hypo7 | Professional | -0.300000000 | -231.809 | 1.110746 | |
| hypo8 | a02e7c2 baseline | -0.356567823 | -217.671 | 1.302444 | -230.58 |
| hypo8 | R2 | **-0.300000000** | -219.621 | 1.303106 | -230.58 |
| hypo8 | Professional | -0.300000000 | -230.588 | 1.284286 | |
| hypo9 | a02e7c2 baseline | -0.356567823 | -230.793 | -5.77e-05 | -231.81 |
| hypo9 | R2 | **-0.300000000** | -235.580 | -5.77e-05 | -231.81 |
| hypo9 | Professional | -0.300000000 | -231.811 | 4.23e-05 | |

Reading: R2 makes the REPORTED strain measure exactly the Professional's
LINEAR engineering strain for all three tests. The stress/void-ratio
convergence is NOT systematically improved: hypo7/hypo9 move away from the
target (-235.6 vs -231.8; baseline was -230.8) and hypo8 improves only
marginally (-219.6 vs -230.6; baseline -217.7). The Professional's void
ratio for hypo7 (1.110746) sits much closer to the baseline LOG run
(1.110389) than to the R2 engineering run (1.099038), so the Professional
appears to report engineering strain while its kernel state follows the
logarithmic measure - i.e. the "log vs linear" hypothesis for the stress
residual is NOT supported by this experiment (it confirms the refutation
already recorded in SEGUIMIENTO 2026-09-08, DEV-C kernels masin). R2 is
therefore a REPORTED-STRAIN-MEASURE correction, not the stress-residual
fix.

- hypo1-4 (`-updated_without_rotation` + `-fixed_in_space`): `.dbs`
  byte-identical to the a02e7c2 baseline (the gate is exactly
  `-UPDATED_LINEAR`).
- sfnet 197-test corpus: identical PASS/RUNFAIL/PARSE classification as
  the a02e7c2 baseline (160 PASS / 37 RUNFAIL / 0 PARSE); validation-suite
  internal suite 16/16.
- Parked basic-Masin candidate + R2: still diverges
  (hypo7 `-737.68` vs `-231.81`); feeding the engineering increment does
  NOT rescue the basic route. R2 reproduces the `-fixed_in_space` strain
  measure exactly (`-737.681616923` vs the R1 `-737.682713253`), which
  validates the kinematics.
- hypo12/13 (no memory record, default): the displacement restriction no
  longer fires; the models run (the masin visco kernel overflows the
  void ratio - pending kernel calibration).
