# control_print_dof_line_move

## Implementación

- Stored as a CONTROL INTEGER record (`-yes`/`-no`). Applied in
  `print_dof_line_point()` of `print_dl.cc` AFTER the files of the call
  are written (print at the current position, THEN move):
  ```
  if ( is_line && move==-YES && materi_velocity && nfound>0 ) {
    if ( db( DTIME, 0, idum, &dtime, ldum, VERSION_NEW, GET_IF_EXISTS ) ||
         db( DTIME, 0, idum, &dtime, ldum, VERSION_NORMAL, GET_IF_EXISTS ) ) {
      for ( ipoint=0; ipoint<npoints; ipoint++ ) {
        if ( !point_found[ipoint] ) continue;
        for ( idim=0; idim<ndim; idim++ )
          coordinates[ipoint*ndim+idim] +=
            point_dof[ipoint*MUKNWN+vel_indx+idim*nder_]*dtime;
      }
      db( CONTROL_PRINT_DOF_LINE_COORDINATES, icontrol, idum, coordinates,
        ncoord, VERSION_NORMAL, PUT );
    }
  }
  ```

## Diseño / decisiones

- Pattern of `POST_POINT_MOVE` (post.cc:59-66): the displacement is the
  INTERPOLATED velocity at the point times `DTIME`. `DTIME` is looked up
  in VERSION_NEW first, then VERSION_NORMAL (same as post.cc).
- Only meaningful when `materi_velocity` is initialized (manual 6.278);
  when it is not, the move is SKIPPED silently (decision; no error).
- The updated coordinates are stored back into the
  `control_print_dof_line_coordinates` record (the record then holds the
  CURRENT particle positions, so the move persists across calls/steps).
- Only the LINE family has `_move` (the Professional manual defines none
  for the point family).

## Detalles

- Points not found in the current call are not moved (their position is
  undefined).
- GOTCHA of the test (dpline_move): the two `control_timestep` blocks of
  one .dat run in SERIES, so the no-move block's prints happen at later
  times; the discriminator is the y-coordinate of the moved point, not
  the value.

## Pendiente

- None.
