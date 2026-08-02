Decisiones y criterios futuros

Memoria
- Mantener wrappers get_new_* con validaciones y delete[] explícitos. Evitar mezclar malloc/free con new/delete.
- Mantener db_close() en rutas nominales. En errores, conservar salida inmediata; habilitar cierre opcional bajo DEBUG para herramientas.

Cadenas
- Estandarizar uso de snprintf/strncat con límites MCHAR-1; prohibir strcpy/strcat en nuevas contribuciones.
- Verificar tamaño requerido antes de concatenar; fallar con mensaje claro si excede.

Paralelismo
- No asignar memoria ni modificar metadatos de DB durante parallel_active; usar buffers locales por hilo. Toda escritura compartida protegida por mutex.

Build/Calidad
- Objetivos con sanitizers y warnings estrictos obligatorios en CI. Permitir desactivar via flags cuando dependencias de terceros no soportan.

Registro
- Mantener informes de auditoría bajo ProjectDocs/ con fecha y alcance.

