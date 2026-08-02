Plan de trabajo (auditoría y endurecimiento)

Objetivo: Reducir riesgos de memoria, desbordes de cadena y condiciones de carrera; facilitar auditorías futuras.

1) Endurecer manejo de cadenas
- Sustituir strcpy/strcat por snprintf/strncat en: bounda.cc, calcul.cc, print_*.cc, database.cc (nombres). Añadir checks de longitud (MCHAR).

2) Sanitizers y warnings
- Añadir objetivos make: audit/asan/ubsan con: -fsanitize=address,undefined -fno-omit-frame-pointer -g -O1, y -Wall -Wextra -Wpedantic.
- Documentar ejecución y exclusiones cuando faltan dependencias opcionales.

3) Paralelismo
- Revisar funciones invocadas dentro de parallel_sys_routine: calcul, contact, elem, post, map, geometry. Catalogar escrituras compartidas y envolver en locks o usar buffers por hilo si aplica.

4) Error/Salida
- Añadir opción DEBUG para llamar db_close() en exit_tn_on_error() y facilitar análisis con herramientas (sin cambiar comportamiento por defecto).

5) Validación
- Correr ejemplos representativos con ASan/Valgrind y registrar reportes en ProjectDocs/ (logs y resúmenes).
- Medir overhead y evaluar trade-offs.

Hitos
- Semana 1: strings + objetivo make audit.
- Semana 2: revisión de paralelismo crítico (contact, elem, post) y fixes.
- Semana 3: validación con sanitizers y reporte final.

