

# OpenLB

**OpenLB** es un paquete de software de código abierto (open-source) y una librería en **C++** de alto rendimiento, diseñada específicamente para implementar simulaciones de fluidos y fenómenos de transporte utilizando el **Método Lattice Boltzmann (LBM)**.

OpenLB es una herramienta reconocida en la comunidad de investigación y fue una de las primeras plataformas genéricas de LBM disponibles para la comunidad (licencia GPLv2).

## Características Clave de OpenLB

### 1. Enfoque y Diseño

- **Orientado a Objetos (C++):** El código está escrito en C++ de manera modular y bien legible, lo que facilita tanto a los programadores de aplicaciones como a los desarrolladores avanzados implementar nuevos modelos LBM personalizados o extender las funcionalidades existentes.
- **Kernel LBM Completo:** El núcleo de OpenLB se basa en una amplia variedad de modelos Lattice Boltzmann, lo que permite abordar una gran diversidad de problemas físicos.

### 2. Alto Rendimiento y Paralelización

Una de las fortalezas más grandes de OpenLB es su arquitectura optimizada para la computación de alto rendimiento:

- **Paralelismo Híbrido:** Es eficiente en plataformas de memoria compartida y distribuida. Soporta:
  - **MPI (Message Passing Interface):** Para paralelismo distribuido (clusters de CPU).
  - **OpenMP:** Para paralelismo de memoria compartida (multihilo en una sola CPU).
  - **CUDA:** Para el uso de **GPUs (Unidades de Procesamiento Gráfico)**, permitiendo simulaciones extremadamente rápidas.
- **Vectorización (SIMD):** Incluye optimizaciones para el procesamiento de datos a nivel de instrucción, mejorando la eficiencia en CPUs.

### 3. Pre-procesamiento y Geometría

A diferencia de muchos otros *solvers* de CFD, OpenLB tiene capacidades robustas de pre-procesamiento integradas:

- **Generación de Malla Automatizada:** Permite utilizar archivos de geometría en formato **STL** (muy común en CAD) o formas geométricas primitivas (cilindros, esferas, etc.). OpenLB genera automáticamente la malla de volumen (voxelización) adaptada a esa geometría.
- **Condiciones de Contorno:** Soporta la configuración automática de condiciones de contorno complejas.

### 4. Áreas de Aplicación

Debido a su diseño y flexibilidad, OpenLB aborda una vasta gama de problemas de transporte y fluidos, similar a Palabos, pero con una posible mayor especialización en ciertos dominios:

- **Flujos con Geometría Compleja:** Ideal para simular flujos en componentes de ingeniería, medicina o medios porosos.
- **Flujos Multifásicos y Multicomponentes:** Manejo de la interacción entre diferentes fluidos.
- **Flujos Térmicos y Radiación:** Modelado de transferencia de calor y radiación luminosa.
- **Flujos Turbulentes:** Incluye modelos de turbulencia validados.
- **Optimización Topológica y Flujo de Partículas:** Soporta métodos Euler-Euler y Euler-Lagrange.

En resumen, **OpenLB** es un *framework* de LBM avanzado, modular y altamente optimizado, especialmente enfocado en ofrecer soluciones eficientes para problemas complejos de CFD y multi-física, aprovechando las últimas capacidades de cómputo paralelo (CPU y GPU).

# Instalar, configurar y testear OpenLB 

Instalar, configurar y testear **OpenLB** en Linux es un proceso similar a la instalación de cualquier librería científica de C++: se descargan las dependencias, se obtiene el código fuente y se compila usando un sistema de construcción (típicamente CMake o Make).

Aquí se tiene una guía paso a paso, priorizando el uso de **CMake**, que es el método recomendado y más moderno.

## 1. Requisitos y Dependencias

Asegúrarse de tener instalados los siguientes paquetes en tu sistema Linux.

| **Requisito**           | **Propósito**                                          | **Comando de Instalación (Ej. Ubuntu/Debian)**        |
| ----------------------- | ------------------------------------------------------ | ----------------------------------------------------- |
| **Compilador C++**      | Necesario para compilar el código.                     | `sudo apt update && sudo apt install build-essential` |
| **CMake**               | Sistema de construcción recomendado por OpenLB.        | `sudo apt install cmake`                              |
| **MPI (OpenMPI/MPICH)** | Esencial para el paralelismo (distribución de tareas). | `sudo apt install openmpi-bin libopenmpi-dev`         |
| **Git**                 | Para clonar el repositorio de OpenLB.                  | `sudo apt install git`                                |

------

## 2. Instalación de OpenLB

### A. Obtener el Código Fuente

Clonar el repositorio oficial de OpenLB usando Git:

Bash

```
git clone https://github.com/OpenLB/openlb.git
```

Esto creará una carpeta llamada `openlb`.

### B. Creación del Directorio de Compilación (Out-of-Source Build)

Es una práctica recomendada realizar la compilación en un directorio separado del código fuente (`out-of-source build`). Esto mantiene el código fuente limpio.

Bash

```
cd openlb
mkdir build
cd build
```

### C. Configuración con CMake

Ejecutar `cmake` para configurar el proceso de construcción, decirle dónde encontrar las dependencias y definir las opciones de compilación.

- **Configuración Básica:**

  Bash

  ```
  cmake ..
  ```

- Configuración con Paralelismo (MPI):

  Si se quiere asegurar que la librería se construya con soporte para MPI (lo cual es crucial para el rendimiento), puedes especificarlo explícitamente (aunque a menudo CMake lo detecta automáticamente):

  Bash

  ```
  cmake -DWITH_MPI=ON ..
  ```

  > **Nota:** Puede añadir otras opciones aquí, como `-DWITH_CUDA=ON` si quieres soporte para GPU (requiere tener el kit de desarrollo CUDA instalado).

### D. Compilación

Una vez que CMake ha generado los `Makefiles`, puedes compilar la librería:

Bash

```
make -j <N>
```

Donde `<N>` es el número de núcleos (threads) que quieres usar para acelerar la compilación (ej., `make -j 8`).

## 3. Testeo y Verificación

Para verificar que OpenLB se instaló correctamente y que las dependencias como MPI funcionan, puedes ejecutar los tests de la librería.

### A. Ejecutar los Tests de Unidad

La forma más directa de probar la funcionalidad de la librería es mediante el comando `ctest`, que fue configurado por CMake:

Bash

```
ctest
```

Este comando ejecutará una serie de pruebas automatizadas incluidas en OpenLB. Si todos los tests pasan, verás un resultado de éxito al final.

### B. Ejecutar un Ejemplo Simple

Para un test funcional en un problema real, puedes compilar y ejecutar uno de los ejemplos que vienen con la distribución.

1. **Encuentra un Ejemplo:** Navegar a la carpeta de ejemplos (ej. `cd ../examples/poiseuille/`).

2. **Compila el Ejemplo:** Dentro de la carpeta del ejemplo, puede usar un `Makefile` específico (si está presente) o usar el sistema de construcción de CMake si el ejemplo está integrado.

3. Ejecución Paralela (Test Principal):

   Una vez compilado el ejecutable (ej. poiseuille), probar con mpirun:

   Bash

   ```
   mpirun -np 4 ./poiseuille
   ```

   - **`-np 4`**: Ejecutar la simulación usando 4 procesos (núcleos de CPU).

Si la simulación se ejecuta hasta el final sin errores de MPI o de segmentación, y genera archivos de salida (típicamente VTK para post-procesamiento), significa que OpenLB está instalado y configurado correctamente.