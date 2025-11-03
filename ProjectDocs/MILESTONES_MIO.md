# Plan de Modernización a C++ Moderno

Este documento describe el plan integral para modernizar el proyecto Tochnog a C++ moderno. Cada hito está diseñado para ser un paso incremental que mejora la calidad del código sin romper la retrocompatibilidad.

## Marco General

El objetivo es reemplazar patrones C heredados con patrones C++ modernos, manteniendo siempre la retrocompatibilidad. El plan sigue el principio de evolución gradual para permitir una transición segura.

---

## Hito 1: Garantizar la Liberación de Recursos de la Base de Datos (COMPLETADO)

**Estado:** COMPLETADO

**Objetivo:** Asegurar que la función `db_close()` se llame siempre al final de la ejecución del programa para prevenir fugas de memoria.

**Tareas completadas:**

1.  **Escribir una prueba de integración que falle:** Crear una prueba que ejecute el programa y verifique que `db_close()` se llama. Esta prueba fallará inicialmente.
2.  **Refactorizar `tn.cc`:** Modificar la función `main` para que llame a `db_close()` antes de terminar, utilizando un enfoque que garantice la llamada incluso en caso de error (por ejemplo, `std::atexit`).
3.  **Añadir opción DEBUG para error paths:** Modificar `exit_tn_on_error()` para que llame opcionalmente a `db_close()` cuando se compila con `DEBUG_DB_CLOSE` definido.
4.  **Verificar que la prueba pasa:** Ejecutar la prueba de integración para confirmar que `db_close()` se llama correctamente.

**Implementación realizada:**
- Se añadió `#include <cstdlib>` y `std::atexit(db_close)` en `tn.cc` para garantizar la llamada a `db_close()` en la terminación normal del programa.
- Se añadió una sección condicionada por `#ifdef DEBUG_DB_CLOSE` en `exit_tn_on_error()` para llamar a `db_close()` en situaciones de error cuando se compila con la opción de depuración.

**Retrocompatibilidad:** Este cambio no afectará a la funcionalidad del programa, pero mejorará significativamente la gestión de la memoria al garantizar que los recursos de la base de datos se liberen siempre.

---

### Hito 2: Refactorizar la Gestión de Memoria de la Base de Datos (COMPLETADO)

**Estado:** COMPLETADO

**Objetivo:** Reemplazar la gestión manual de memoria en `database.cc` con contenedores de la STL para eliminar las fugas de memoria y mejorar la eficiencia.

**Tareas completadas:**

1.  **Escribir pruebas unitarias que fallen:** Crear pruebas que verifiquen el comportamiento de la asignación y liberación de memoria en `database.cc`. Estas pruebas fallarán con la implementación actual.
2.  **Refactorizar `database.cc`:** Reemplazar los arrays gestionados manualmente (`new[]`/`delete[]`) por `std::vector`. Esto simplificará el código, eliminará las fugas de memoria y mejorará la eficiencia del redimensionamiento de arrays.
3.  **Verificar que las pruebas pasan:** Ejecutar las pruebas unitarias para confirmar que la nueva implementación de la base de datos funciona correctamente.

**Implementación realizada:**
- Sustitución de arrays globales por contenedores `std::vector`
- Actualización de `db_allocate()` para usar `vector::resize()` en lugar de `new[]`/`delete[]`
- Actualización de `db_delete()` y `db_close()` para usar `vector::clear()` en lugar de `delete[]`
- Actualización de funciones de acceso como `db_int()` y `db_dbl()` para trabajar con vectores
- Las funciones mantienen la misma interfaz para garantizar retrocompatibilidad

**Retrocompatibilidad:** Este cambio es interno a la implementación de la base de datos y no debería afectar a la API pública. Las pruebas de integración existentes garantizarán que no se introduzcan regresiones.

---

### Hito 3: Encapsular la Base de Datos en una Clase (COMPLETADO)

**Estado:** COMPLETADO

**Objetivo:** Mejorar la organización del código y reducir el uso de variables globales encapsulando la lógica de la base de datos en una clase.

**Tareas completadas:**

1.  **Crear una clase `Database`:** Diseñar una clase que contenga los datos y métodos de la base de datos.
2.  **Mover la lógica a la clase:** Migrar gradualmente las funciones y variables globales de `database.cc` a la nueva clase `Database`.
3.  **Refactorizar el código para usar la clase:** Actualizar el resto del código para que interactúe con la base de datos a través de la nueva clase.

**Implementación realizada:**
- Creación de `DatabaseClass.h` y `DatabaseClass.cpp` con la implementación completa de la clase Database
- Implementación de todos los métodos de la base de datos dentro de la clase
- Mantenimiento de retrocompatibilidad mediante funciones globales que delegan a la instancia de la clase
- Creación de `DatabaseCompatibility.h` para mantener las interfaces originales
- Agregación de un test unitario para verificar la funcionalidad de encapsulación

**Retrocompatibilidad:** Se mantiene completa retrocompatibilidad con las interfaces existentes mediante el uso de una instancia global de la clase Database que es utilizada por las funciones de compatibilidad.

---

### Hito 4: Implementar RAII y Smart Pointers (COMPLETADO)

**Estado:** COMPLETADO

**Objetivo:** Reemplazar la gestión manual de recursos con principios RAII (Resource Acquisition Is Initialization) y smart pointers para eliminar fugas de memoria y mejorar la seguridad.

**Tareas completadas:**

1.  **Identificar recursos manejados manualmente:** Buscar todas las instancias de `new`/`delete`, `malloc`/`free`, apertura/cierre de archivos, etc.
2.  **Reemplazar con smart pointers:** Usar `std::unique_ptr`, `std::shared_ptr` y `std::make_unique`/`std::make_shared` donde sea apropiado.
3.  **Crear RAII wrappers:** Para recursos que no tienen wrappers STL (manejo de arrays, etc.).
4.  **Actualizar el código para usar RAII:** Mover el manejo de recursos al ámbito apropiado con destructores automáticos.
5.  **Verificar que no hay fugas:** Asegurar la liberación automática de recursos.

**Implementación realizada:**
- Creación de `raii_resources.h` con wrappers RAII para arrays
- Implementación de `SafeArray<T>` plantilla para manejo seguro de arrays
- Funciones `make_safe_*_array()` que usan `std::unique_ptr`
- Aplicación de principios RAII en la adquisición y liberación automática de recursos

**Retrocompatibilidad:** El cambio se implementó principalmente en la capa de implementación interna, manteniendo interfaces externas compatibles.

---

### Hito 5: Modernizar el Manejo de Cadenas (COMPLETADO)

**Estado:** COMPLETADO

**Objetivo:** Reemplazar patrones C de manejo de cadenas con `std::string` y `std::string_view` para mejorar seguridad y conveniencia.

**Tareas completadas:**

1.  **Identificar uso de arrays C (`char[MCHAR]`) y funciones `strcpy`/`strcat`/etc.:** Reemplazar con `std::string` donde sea posible.
2.  **Actualizar interfaces que usan cadenas:** Donde sea seguro, cambiar firmas de funciones para usar `const std::string&` o `std::string_view`.
3.  **Implementar funciones de utilidad de cadenas seguras:** Reemplazar la lógica manual de manipulación de cadenas.
4.  **Mantener compatibilidad:** Donde sea necesario, proporcionar conversiones entre `char*` y `std::string`.

**Implementación realizada:**
- Creación de `string_utils.h` con funciones de utilidad de cadenas seguras
- Funciones `safe_string_copy`, `safe_string_concat` como reemplazo seguro de `strcpy`/`strcat`
- Clase `FileNameBuilder` para construcción segura de nombres de archivo
- Función `long_to_string_modern` como versión moderna de `long_to_a`
- Utilidades para conversión segura entre `std::string` y arrays C

**Retrocompatibilidad:** Mantener interfaces que trabajan con `char*` para compatibilidad, pero fomentar el uso de `std::string` internamente.

---

### Hito 6: Eliminar Macros en Favor de Constantes y Funciones (COMPLETADO)

**Estado:** COMPLETADO

**Objetivo:** Reemplazar macros peligrosas con constantes `const`/`constexpr` y funciones `inline`.

**Tareas completadas:**

1.  **Identificar macros de constantes:** Reemplazar con `const`/`constexpr`.
2.  **Identificar macros de funciones pequeñas:** Reemplazar con funciones `inline` o `constexpr`.
3.  **Usar `enum class` en lugar de `#define` para constantes relacionadas.**

**Implementación realizada:**
- Creación de `modern_constants.h` con constantes modernas
- Namespace `TochnogConstants` con todas las constantes como `constexpr`
- Funciones `constexpr`/`inline` como reemplazo seguro para macros funcionales
- Plantillas para funciones como `modern_min`, `modern_max`, `safe_abs`
- Mantenimiento de compatibilidad con macros originales

**Retrocompatibilidad:** Mantener macros existentes para compatibilidad, pero usar las versiones modernas internamente.

---

### Hito 7: Implementar Excepciones para Manejo de Errores (COMPLETADO)

**Estado:** COMPLETADO

**Objetivo:** Reemplazar el manejo de errores basado en códigos de retorno y `exit()` con un sistema robusto de excepciones.

**Tareas completadas:**

1.  **Definir jerarquía de excepciones:** Crear una jerarquía de excepciones específicas para Tochnog.
2.  **Reemplazar usos de `exit()` por lanzamiento de excepciones:** Donde sea apropiado, lanzar excepciones en lugar de salir directamente.
3.  **Actualizar código para manejar excepciones:** Agregar bloques `try`/`catch` donde sea necesario.
4.  **Mantener puntos de entrada compatibles:** Para mantener retrocompatibilidad, convertir excepciones a códigos de retorno en capas superiores si es necesario.

**Implementación realizada:**
- Creación de `tochnog_exceptions.h` con jerarquía completa de excepciones
- Excepciones específicas: `DatabaseException`, `MemoryException`, `SolverException`, etc.
- Actualización de `get_new_*` functions para lanzar `OutOfMemoryException`
- Actualización de `db_error` para lanzar `DatabaseException`
- Sistema opcional con `#ifdef USE_EXCEPTIONS` para transición gradual
- Macro `TOCHNOG_TRY_BEGIN`/`TOCHNOG_CATCH_END` para manejo compatible

**Retrocompatibilidad:** Este cambio se implementó cuidadosamente manteniendo compatibilidad opcional.

---

### Hito 8: Modernizar Concurrencia con `std::thread`

**Objetivo:** Reemplazar el uso directo de pthread con herramientas de concurrencia modernas de C++.

**Tareas:**

1.  **Identificar código de pthread actual:** Buscar usos de `pthread_*` y `mutex` C.
2.  **Reemplazar con `std::thread`, `std::mutex`, `std::lock_guard`, etc.:** Mantener la funcionalidad pero usar el modelo de C++.
3.  **Usar RAII para sincronización:** `std::lock_guard`, `std::unique_lock`, etc.

**Retrocompatibilidad:** La interfaz externa puede mantenerse igual mientras se moderniza la implementación interna.

---

### Hito 9: Implementar Templates y Genéricos

**Objetivo:** Reemplazar macros y duplicación de código con templates para mejorar la seguridad y reutilización.

**Tareas:**

1.  **Identificar código duplicado que varía por tipo:** Funciones que se repiten para diferentes tipos.
2.  **Reemplazar con templates:** Donde sea apropiado, usar funciones y clases genéricas.
3.  **Usar `auto` y deducción de tipos:** Para simplificar el código y mejorar legibilidad.

**Retrocompatibilidad:** Los templates no afectan la interfaz externa, manteniendo compatibilidad.

---

### Hito 10: Usar `enum class` en lugar de `enum` (COMPLETADO)

**Estado:** COMPLETADO

**Objetivo:** Reemplazar enums C con enums de clase para mejorar el ámbito y la seguridad de tipos.

**Tareas completadas:**

1.  **Identificar `enum` usados:** Crear `enum class` equivalentes donde sea seguro hacerlo.
2.  **Actualizar código para el nuevo tipo:** Ajustar donde se usan los valores de enum.
3.  **Mantener compatibilidad:** Donde sea necesario para interfaces externas, proporcionar conversiones.

**Implementación realizada:**
- Creación de `enum_modernization.h` con `enum class` equivalentes
- `VersionType` y `DataType` como enums de clase para versiones y tipos de datos
- Funciones de conversión `to_legacy_*` y `from_legacy_*` para mantener compatibilidad
- Proporciona seguridad de tipos y ámbito fuerte

**Retrocompatibilidad:** Se mantiene compatibilidad a través de conversiones explícitas mientras se moderniza internamente.

---
