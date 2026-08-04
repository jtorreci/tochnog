# post_point_move

## Archivos y funciones

- `post.cc` → `post()` (`post.cc:30`) — the per-post-point loop that moves
  the point (`post.cc:48-69`).
- `post.cc` → `parallel_post_point()` (`post.cc:234`) — locates the element
  containing the post point and interpolates the nodal DOFs to the point
  (used by the move logic to obtain the velocity at the point).
- `database.cc:3677-3681` — keyword registration:
  `strcpy(name[POST_POINT_MOVE], "post_point_move")`, type `INTEGER`,
  `data_length = 1`, `no_index = 1`, class `POST`.
- Enum `POST_POINT_MOVE` in `tochnog.h` / `tochnog-mod.h`.

## Detalles de implementación

- `post()` reads the switch once with
  `db( POST_POINT_MOVE, 0, &post_point_move, ddum, ldum, VERSION_NORMAL,
  GET_IF_EXISTS )` (`post.cc:44-45`).
- For each active `POST_POINT`, it resets `post_point_dof` to 0,
  `post_found = 0`, and runs `parallel_sys_routine( &parallel_post_point )`
  (`post.cc:51-53`). `parallel_post_point` walks the elements with
  `point_el()`, locks, and accumulates
  `post_point_dof += weight[inol] * node_dof` (`post.cc:259-273`).
- If the point was found, the interpolated DOFs are stored with
  `db( POST_POINT_DOF, ipost, ... )` (`post.cc:55-56`).
- Move logic (`post.cc:59-66`): only when
  `post_point_move == -YES && materi_velocity`:
  - read `DTIME` preferring `VERSION_NEW`, falling back to
    `VERSION_NORMAL` (`post.cc:60-61`),
  - for each spatial axis `i`: `post_point[i] += post_point_dof[vel_indx
    + i*nder] * dtime` — the velocity is the `vel_indx` DOF of the
    interpolated vector,
  - persist with `db( POST_POINT, ipost, idum, post_point, ndim,
    VERSION_NORMAL, PUT )`.
- `vel_indx` is a global initialized to `-1` in `initia.cc:52` and set
  during unknown setup.

## Dependencias externas

- `point_el()` (`point_el.cc`) for element location, the parallel system
  (`parallel_sys_routine`, `parallel_sys_lock`/`unlock`), core database
  accessors, and the global `materi_velocity`, `vel_indx`, `nder`, `dtime`.

## Parámetros hardcodeados / refactorizaciones pendientes

- IMPORTANT: `point_el()` has pre-existing bugs for truss/beam elements and
  for some quad4 positions; a post point (and therefore `post_point_move`)
  in those elements may be mislocated. Fixing `point_el()` benefits both
  this feature and `force_point`.
- The `VERSION_NEW`/`VERSION_NORMAL` fallback for `DTIME` is duplicated in
  several routines; a small `get_dtime()` helper would remove the pattern.
- `post_point_move` is re-read for every `post()` call; it could be cached
  in a global like `check_used` if repeated reads matter.
- The condition `post_point_move == -YES` uses the global `-YES` constant;
  an explicit bool would be clearer.
