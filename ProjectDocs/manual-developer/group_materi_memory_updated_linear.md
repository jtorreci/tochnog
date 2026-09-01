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

## Verification

- hypo7/hypo8 (explicit `-updated_linear`): run and converge to within
  3-6% of the Professional targets (masin kernel accuracy pending).
- hypo12/13 (no memory record, default): the displacement restriction no
  longer fires; the models run (the masin visco kernel overflows the
  void ratio - pending kernel calibration).
