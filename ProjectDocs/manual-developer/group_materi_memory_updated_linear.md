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
   The stress-rotation block ("rotate to new configuration") does NOT
   rotate `-UPDATED_LINEAR` (fixed 2026-09-29, audit cluster K01): only
   `-UPDATED` (by `inc_rot`) and `-TOTAL`/`-TOTAL_PIOLA` (by `new_rot`)
   rotate. Commit `22801c3` (2026-09-01) had added `-UPDATED_LINEAR` to
   that list on the assumption that `inc_rot` is the identity for it; it
   is not, because `set_deften_etc()` overwrites the identity with
   `set_deften_u_rot(inc_deften, inc_u, inc_rot)` (called for every
   memory). The stress was rotated by the real incremental rotation. The
   Professional manual 2.2.1 defines `-updated_linear` as "add the extra
   stresses to the stresses of the previous time point" with no rotation
   step (only `-updated` rotates by R); the 2014 source also rotated only
   `-updated`/`-total`/`-total_piola`. See section 5.

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
   if ( memory==-TOTAL_LINEAR || memory==-UPDATED_LINEAR ||
        memory==-UPDATED_AREA ) {
     tmp_new = ( new_epe[0] + new_epe[4] + new_epe[8] ) / compressibility;
     tmp_inc = ( inc_epe[0] + inc_epe[4] + inc_epe[8] ) / compressibility;
   }
   ```
   Until 2026-09-29 `tmp_inc` was `tr(new_epe - old_epe)/compressibility`.
   `old_epe` is loaded only when `materi_strain_elasti` is initialized
   (materi.cc, the `old_epe` load from the dof); without it `old_epe` is 0
   and the TOTAL elastic strain was added to the stress every step, so a
   compressibility-only material became far too stiff (Pro `interac1`:
   element 2 strain ~0 with stress -0.01 instead of -0.005 / -0.005).
   Legacy from the 2014 source. `tmp_new` is only used with
   `formulation==TOTAL`.

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

5. **Audit cluster K01 fix (2026-09-29): no rotation + compressibility
   increment.** Both changes above (materi.cc rotation list, stress.cc
   `tmp_inc`). Evidence against the Professional 02-08-2026 (all 13 K01
   tests re-run, every variable compared; table in
   `ProjectDocs/AUDIT-ALL-VARIABLES-2026.md`):

   | model | variable | Pro | GNU before | GNU after |
   |---|---|---|---|---|
   | Pro `interac1` (2 bars, compressibility, analytic -0.005) | sigxx n1 / n3, eptxx n3 | -0.005 / -0.005 / -0.005 | -0.01 / -0.01 / ~0 | -0.005 / -0.005 / -0.005 |
   | sfnet `druckpr1` (simple shear gamma 1, DP c 0.1, analytic 0.11547) | sigxx / sigxy / sigyy | 0 / 0.11547 / 0 | 0.0924 / 0.0693 / -0.0924 | 0 / 0.11547 / 0 |
   | same, 10 steps | sigxx / sigxy | 0 / 0.11547 | 0.0365 / 0.1096 | 0 / 0.11547 |

   The pre-fix `druckpr1` answer is the correct stress rotated by
   atan(gamma/2) = 26.57 deg; a smaller time step did not remove it. With
   `-updated_without_rotation` the pre-fix GNU was already exact.
   Self-checking inputs: `validation-suite/test-2026/updated_linear_compress1.dat`,
   `updated_linear_shear1.dat` (pass on the fixed GNU and the Professional,
   rc=1 on the pre-fix binary).

6. **Reference vs current gradient: what the Professional does
   (measured 2026-09-29, NOT changed).** The manual says "the linear
   engineering strains in the DEFORMED configuration are used as incremental
   strain", while section 4 (commit `c945833`) builds the `-UPDATED_LINEAR`
   increment with the REFERENCE gradient. Discriminating probe: one `bar2`,
   E = 1, left end fixed, right end velx = 1 for t in [0, 0.5], dt 0.01
   (50% extension; reference gradient gives eps = 0.5, current gradient
   gives sum(du/l) ~ ln(1.5) = 0.405):

   | mesh record | memory | Pro | GNU |
   |---|---|---|---|
   | none (Pro default `-fixed_in_space`, GNU default `-follow_material`) | `-updated_linear` | 0.5 | 0.5 |
   | none | `-updated` | 0.5 | 0.40714 |
   | `mesh -follow_material` | `-updated_linear` | **0.40714** | **0.5** |
   | `mesh -follow_material` | `-updated` | 0.40714 | 0.40714 |
   | `mesh -fixed_in_space` | `-updated_linear` | 0.5 | 0.50003 |

   Reading: the Professional uses the CURRENT-geometry gradient for
   `-updated_linear` (as its manual says); on a mesh that follows the
   material it gives exactly the GNU `-updated` value. Its engineering
   0.5 without a mesh record comes from its DEFAULT mesh, which is
   `-fixed_in_space` (manual 6.806: "Default each specifier is set to
   -fixed_in_space"), so the geometry never moves. The GNU defaults to
   `-FOLLOW_MATERIAL` whenever `materi_velocity` is initialized without
   `materi_displacement` (`top.cc`, the `options_mesh` default), and the
   reference-gradient path of `c945833` emulates the Professional's fixed
   mesh for `-updated_linear` only. Consequence: on an explicit
   `-follow_material` mesh the GNU `-updated_linear` differs from the
   Professional (0.5 vs 0.40714). The underlying divergence is the mesh
   default, not the memory; tracked as BACKLOG C8. The same mesh default
   explains the `interac1` `-updated` residual (GNU -0.0050117, Pro
   -0.005): with an explicit `mesh -follow_material` both give -0.0050117,
   with `mesh -fixed_in_space` Pro -0.005 and GNU -0.0050000474. The GNU
   `-updated` value is the one consistent with the manual's `-updated`
   definition (incremental stretch in the current configuration, so the
   accumulated strain is logarithmic-like: ln(1/(1-0.005)) = 0.0050125).
