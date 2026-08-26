# control_print_number_iterations

## Implementación

- **Keyword**: `control_print_number_iterations` (INTEGER, data_length
  1, data_class CONTROL) registrado en `database.cc`. `ival[0]` es el
  switch (`-yes` / `-no`).
- **Lógica**: en `top()` (`top.cc`), justo ANTES del bucle de
  iteraciones de equilibrio (for de `iteration`), se lee el record del
  `icontrol` actual con `GET_IF_EXISTS` en una variable local
  INICIALIZADA a `-NO` (regla GET_IF_EXISTS):
  ```
  long int print_number_iterations = -NO;
  db( CONTROL_PRINT_NUMBER_ITERATIONS, icontrol,
    &print_number_iterations, ddum, ldum, VERSION_NORMAL, GET_IF_EXISTS );
  for ( iteration=1; ...; iteration++ ) {
    if ( print_number_iterations==-YES )
      cout << "control_print_number_iterations: time " << time_current
           << " iteration " << iteration << "\n";
    ...
  }
  ```
  Se imprime UNA línea por iteración de equilibrio (monitor de consola,
  manual 6.336: "print the iteration number while doing equilibrium
  iterations in a time step").
- **No confundir con `-inverse_iteration_number`**: ese data item
  (`INVERSE_ITERATION_NUMBER`, top.cc:221) guarda el nº de iteraciones
  del paso YA terminado para el análisis inverso. El monitor imprime el
  contador VIVO durante las iteraciones.
- **Dónde vive el bucle**: el for de equilibrio está dentro del
  `while ( time_current < time_new )` del bloque `control_timestep`
  (top.cc ~365); `icontrol` es el índice del bloque (record
  `ICONTROL`). El read del record se hace una sola vez por paso (no por
  iteración).

## Enums nuevos

- `CONTROL_PRINT_NUMBER_ITERATIONS` (bloque CONTROL_PRINT, tras
  `CONTROL_PRINT_MATLAB`) en `tochnog.h` / `tochnog-mod.h` (sync).

## Pendiente

- El monitor imprime por iteración de equilibrio; los bucles de
  iteraciones inversas (`inverse_iter`) o de timestep no están
  monitorizados.
