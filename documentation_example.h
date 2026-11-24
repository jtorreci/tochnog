/*
 * Ejemplo de documentación mejorada para funciones complejas
 * Este archivo demuestra cómo se debería documentar la función db()
 * que tiene una firma muy compleja
 */

#include "tochnog.h"

/**
 * @brief Función principal de acceso a la base de datos
 * 
 * Esta es la función central del sistema de base de datos de Tochnog.
 * Proporciona acceso uniforme para obtener, colocar y gestionar datos
 * en diferentes categorías de la base de datos.
 * 
 * @param idat Identificador del tipo de dato en la base de datos
 *             Valores típicos: NODE, ELEMENT, DOF_LABEL, etc. (usar -valor para lectura)
 * @param index Índice dentro de la categoría especificada
 *             Para datos sin índice, usar 0. Rango típico: 0 a número máximo de elementos/nodos
 * @param int_arr Array para almacenar datos enteros (de entrada o salida según acción)
 *                Debe tener tamaño suficiente para 'length' elementos
 *                Si action es GET: aquí se devuelve la información
 *                Si action es PUT: aquí se proporciona la información a almacenar
 * @param dbl_arr Array para almacenar datos dobles (de entrada o salida según acción)  
 *                Similar a int_arr pero para valores double
 *                Debe tener tamaño suficiente para 'length' elementos
 * @param length Referencia a la longitud del registro
 *               Para PUT: debe contener la longitud de los datos a almacenar
 *               Para GET: se actualiza con la longitud real del registro obtenido
 * @param version Versión temporal de los datos a acceder
 *                Valores típicos: VERSION_NORMAL, VERSION_START, VERSION_NEW
 *                Permite manejar diferentes estados temporales de la malla
 * @param action Acción a realizar con los datos
 *               Valores posibles: GET, PUT, GET_IF_EXISTS, GET_AND_CHECK
 *               GET: Obtener datos existentes
 *               PUT: Almacenar nuevos datos (crea el registro si no existe)
 *               GET_IF_EXISTS: Obtener si existe, no error si no existe
 *               GET_AND_CHECK: Obtener y verificar que la longitud coincida
 * 
 * @return long int Código de resultado
 *                1: Operación exitosa
 *                0: Operación fallida (en modo GET_IF_EXISTS si no existe)
 * 
 * @throws DatabaseException Si se compila con USE_EXCEPTIONS y ocurre un error
 * 
 * @note Esta es la función de acceso más utilizada en Tochnog, llamada miles de veces
 * @note La función maneja automáticamente la creación de espacio si action es PUT
 * @warning Asegurarse de que los arrays tengan tamaño suficiente para 'length' elementos
 * @warning Los índices deben estar dentro de rangos válidos para evitar errores
 * 
 * @par Ejemplo de uso para lectura:
 * @code
 *     long int int_data[10];
 *     double dbl_data[20];
 *     long int length = 0;
 *     
 *     // Obtener información de un nodo específico
 *     long int result = db(NODE_DOF, 100, int_data, dbl_data, length, 
 *                         VERSION_NORMAL, GET);
 *     if (result == 1) {
 *         // Datos obtenidos exitosamente, length contiene la longitud real
 *     }
 * @endcode
 * 
 * @par Ejemplo de uso para escritura:
 * @code  
 *     long int new_int_data[] = {1, 2, 3};
 *     double new_dbl_data[] = {1.0, 2.0};
 *     long int length = 3;  // longitud de los datos a escribir
 *     
 *     // Almacenar información para un elemento
 *     long int result = db(ELEMENT_GROUP, 50, new_int_data, new_dbl_data, 
 *                         length, VERSION_NORMAL, PUT);
 * @endcode
 * 
 * @sa db_int - Acceso directo a datos enteros
 * @sa db_dbl - Acceso directo a datos dobles  
 * @sa db_active_index - Verificar si un índice está activo
 * @sa db_max_index - Obtener/establecer el índice máximo
 */
// La firma real de la función se mantiene igual
// long int db( long int idat, long int index, long int int_arr[], 
//              double dbl_arr[], long int &length, long int version, 
//              long int action );

/**
 * @brief Obtiene un puntero directo a datos enteros en la base de datos
 * 
 * Función auxiliar que proporciona acceso directo a los datos enteros
 * almacenados en la base de datos, evitando copias innecesarias.
 * 
 * @param idat Identificador del tipo de dato (siempre positivo)
 * @param index Índice dentro del tipo de dato
 * @param version Versión de los datos a acceder
 * 
 * @return long int* Puntero al inicio de los datos enteros para el registro
 *                  NULL si el índice no es válido o el tipo de dato no es entero
 * 
 * @note El puntero es válido solo mientras la base de datos no cambie
 * @note No se debe liberar manualmente el puntero devuelto
 * @warning Acceder a memoria fuera del rango puede causar errores
 * 
 * @par Ejemplo de uso:
 * @code
 *     // Obtener acceso directo a los grados de libertad de un nodo
 *     long int* node_dof = db_int(NODE_DOF, 100, VERSION_NORMAL);
 *     if (node_dof != NULL) {
 *         // Acceder a los datos directamente
 *         long int first_dof = node_dof[0];
 *         // Modificar datos directamente
 *         node_dof[1] = 42;
 *     }
 * @endcode
 */
// long int* db_int( long int idat, long int index, long int version );

#endif // DOCUMENTATION_EXAMPLE_H