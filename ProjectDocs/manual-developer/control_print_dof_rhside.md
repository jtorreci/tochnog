# control_print_dof_rhside

## Implementación

- **Alias puro**: `control_print_dof_rhside` es el nombre Professional
  del GNU `control_print_unknownsrhside`. Se registra ÚNICAMENTE como
  traducción en `db_number()` (`database.cc`):
  ```
  else if ( !strcmp( str, "control_print_dof_rhside" ) )
    return CONTROL_PRINT_UNKNOWNSRHSIDE;
  ```
  Sin enum nuevo, sin dispatch nuevo (top.cc:905-908 ya despacha
  `CONTROL_PRINT_UNKNOWNSRHSIDE` -> `print_unknownsrhside()`).
- **Por qué en db_number**: el detector de fin-de-valores de records
  variable-length de input.cc (`db_number(str) >= 0`) también llama
  `db_number`; si la traducción viviera solo en el punto del keyword,
  un input con `control_print_dof_rhside` después de un record
  variable-length rompería el parseo (gotcha Sprint 9, fila f178849).
- **Output**: `print_unknownsrhside()` (`print_un.cc`) escribe
  `<dof>_rhside.<index>` (e.g. `temp_rhside.20` para `-temp`, como el
  manual 6.284) con `x y z <rhs>` por nodo. El nombre de archivo sale de
  `db_name(dof_label[iuknwn])` — coincide con el manual.
- El `check` de `CONTROL_PRINT_UNKNOWNSRHSIDE` (check.cc:
  `check_unknowns_are_specified`) aplica automáticamente por enum.

## Pendiente

- Nada específico: el alias comparte el comportamiento (y las
  limitaciones) del GNU `control_print_unknownsrhside`.
