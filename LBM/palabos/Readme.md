# Palabos

Palabos es una **librería de código abierto (open-source)** escrita en **C++** que se utiliza para desarrollar y ejecutar simulaciones de fluidos basadas en el **Método de Lattice Boltzmann (LBM)**, del cual hablamos anteriormente.

El nombre "Palabos" es un acrónimo de **PArLAs-BOlS-Solver** (Parallel Lattice Boltzmann Solver).

### Características Principales

- **Método Base:** Se basa exclusivamente en el Método Lattice Boltzmann (LBM).
- **Lenguaje:** Es una librería escrita en C++, lo que permite a los usuarios escribir sus propias aplicaciones de CFD aprovechando las estructuras y optimizaciones que ofrece Palabos.
- **Paralelización:** Está diseñado desde su base para ser **altamente escalable y eficiente** en sistemas paralelos, utilizando la Interfaz de Paso de Mensajes (MPI) para ejecutarse en clústeres y supercomputadoras. Esto es clave dada la naturaleza local del LBM.
- **Plataforma:** Es un software multiplataforma, compatible con Linux, Windows y macOS.

###  Desarrollo y Soporte

El desarrollo y mantenimiento de Palabos están supervisados principalmente por la empresa **FlowKit Ltd.** y tienen una fuerte colaboración con la **Universidad de Ginebra** en Suiza, que es pionera en la investigación del método Boltzmann.

## Palabos Aplicaciones Típicas de Palabos

Palabos es valorado por su capacidad para manejar simulaciones complejas, donde el LBM demuestra superioridad sobre los métodos CFD tradicionales. Sus principales áreas de aplicación incluyen:

1. **Flujos Multifásicos:** Simulación de la interacción entre diferentes fases de fluidos, como la mezcla de aceite y agua o la dinámica de burbujas y gotas.
2. **Medios Porosos:** Modelado del flujo de fluidos a través de estructuras porosas complejas (p. ej., filtración, pilas de combustible).
3. **Flujos Térmicos:** Simulación de convección natural y forzada, incluyendo la transferencia de calor.
4. **Flujos con Reacciones Químicas:** Permite integrar modelos de reacción dentro de la dinámica del fluido.

En resumen, Palabos es una herramienta esencial para investigadores y profesionales que buscan implementar de forma eficiente simulaciones de fluidos de vanguardia utilizando el enfoque mesoscópico del Lattice Boltzmann.

# Instalar, configurar y testear Palabos

Instalar, configurar y testear **Palabos** en Linux generalmente se realiza compilando el código fuente. Palabos depende de un conjunto de bibliotecas de terceros, siendo **MPI (Message Passing Interface)** crucial para la ejecución paralela.

Aquí se tiene una guía paso a paso para la instalación y configuración básica, seguida de cómo ejecutar el test de la librería.

## 1. Requisitos Previos

Antes de instalar Palabos, asegúrense de tener las siguientes herramientas y librerías instaladas en su distribución Linux (p. ej., Ubuntu, Debian, Fedora):

| **Requisito**           | **Propósito**                                      | **Comando de Instalación (Ej. Ubuntu/Debian)**        |
| ----------------------- | -------------------------------------------------- | ----------------------------------------------------- |
| **Compilador C++**      | Necesario para compilar el código fuente.          | `sudo apt update && sudo apt install build-essential` |
| **MPI (OpenMPI/MPICH)** | Esencial para la ejecución paralela (obligatorio). | `sudo apt install openmpi-bin libopenmpi-dev`         |
| **Git**                 | Para clonar el repositorio de Palabos.             | `sudo apt install git`                                |

------

## 2. Instalación de Palabos

### A. Obtener el código fuente

Utilizar Git para clonar el repositorio de Palabos (o puede descargar el archivo ZIP desde su sitio oficial):

Bash

```
git clone https://gitlab.com/unigespc/palabos.git
```

Esto creará una carpeta llamada `palabos` en tu directorio actual.

### B. Instalación de prerequisitos

Para que palabos funcione correctamente se requieren los siguientes paquetes.

#### Distribuciones derivadas Debian

Bash
```
$ sudo apt install gcc clang clang-format cmake make libtbb-dev
```
### C. Compilación de la Librería

Palabos no requiere un paso de "instalación" tradicional (`make install`); en su lugar, se compila las aplicaciones vinculándolas a la librería Palabos. La configuración se realiza principalmente a través del **Makefile** incluido.

1. **Navega a la carpeta:**

   Bash

   ```
   cd palabos
   ```

2. **Configurar el Makefile:**

   - Palabos utiliza una estructura de directorios modular. Para el test inicial, se puede usar uno de los ejemplos. Navegar a un directorio de prueba, por ejemplo:

     Bash

     ```
     cd examples/showCases/laminarChannel/build
     ```

   - Abrir el `Makefile` con un editor de texto (p. ej. `nano Makefile` o `vim Makefile`).
   - **Verificar la configuración de MPI:** Asegurarse de que las variables de compilación (como `CXX` y `MPICXX`) apunten a los compiladores de MPI correctos (por defecto suelen ser `mpicxx` o `mpic++`, que ya deberían estar configurados si se instalo `libopenmpi-dev`).
   - **Opcional - Optimización:** Puede ajustar las banderas de optimización si lo desea (p. ej. `-O3`).


3. **Compilar el ejemplo**:

   Ejecutar el comando make para compilar el código fuente del ejemplo:

   Bash

   ```
   make
   ```

   Si la compilación es exitosa, se creará un archivo ejecutable (con el mismo nombre que el directorio, o el nombre definido en el Makefile, en este caso probablemente `laminarChannel`).

## 3. Testeo y Ejecución

Una vez construido el ejecutable, se puede probar así:

### A. Ejecución en un Solo Núcleo

Para ejecutar la simulación en un solo núcleo de CPU (secuencial):

Bash

```
./laminarChannel
```

### B. Ejecución Paralela (Usando MPI)

Para probar la funcionalidad paralela (el corazón de Palabos), utiliza el comando `mpirun` o `mpiexec`. Aquí se ejecutará el test en 4 núcleos (cambiar el número según tu CPU):

Bash

```
mpirun -np 4 ./laminarChannel
```

- **`-np 4`**: Indicar a MPI que lance el programa y distribuya la carga de trabajo entre **4** procesos (núcleos).

### C. Verificación de Resultados

Una ejecución exitosa producirá mensajes de progreso en la consola y, generará archivos de salida para la visualización.

- **Archivos de Salida:** El ejemplo `laminarChannel` generará archivos en formato **VTK** (o similar) que contienen la distribución de velocidad y presión. Estos archivos son legibles por software de posprocesamiento como **ParaView**.

Si la simulación se ejecuta hasta el final sin errores de segmentación o MPI y los archivos de salida son generados, la instalación de Palabos es exitosa.

## 4. Configuración Avanzada (Para Proyectos Propios)

Para empezar un propio proyecto basado en Palabos:

1. **Crear un Directorio de Proyecto:** Crear una nueva carpeta fuera de la estructura de `palabos/examples`.

2. **Copiar el Makefile:** Copiar el `Makefile` de alguno de los ejemplos de Palabos al nuevo directorio de proyecto.

3. **Modificar la fuente:** En el `Makefile` del proyecto, ajustar la variable `SRC_FILES` para apuntar al código fuente C++ personalizado (p. ej., `mySimulation.cpp`).

4. **Inclusión:** Asegurarse que el código C++ incluya la cabecera principal de Palabos:

   C++

   ```
   #include "palabos.h"
   ```

  # Adiciones al README de Palabos

## 1. Requisitos Previos

Palabos en sí mismo depende de muy poco: en su forma más básica solo necesita un compilador C++ y, opcionalmente, MPI. La confusión suele venir de mezclar el flujo de compilación **clásico por Makefile** (el que se usa históricamente en cada carpeta de `examples/showCases/...`) con el flujo **CMake**, que es el que usa el repositorio principal (`CMakeLists.txt` en la raíz) y el que recomienda hoy el proyecto para compilar la librería completa o integrarla en un proyecto propio. Ambos flujos existen y son válidos, pero no piden exactamente lo mismo.

### A. Necesario en cualquiera de los dos flujos

| Requisito | Propósito | Comando (Ubuntu/Debian) |
|---|---|---|
| **Compilador C++ (GCC o Clang)** | Compila el código fuente de Palabos y de tus aplicaciones. | `sudo apt install build-essential` (incluye GCC) |
| **Git** | Clonar el repositorio. | `sudo apt install git` |
| **Make** | Backend de compilación, invocado directamente en el flujo de Makefile o generado por CMake en el otro. | `sudo apt install make` |

### B. Solo si compilas con el flujo Makefile clásico (por ejemplo)

Este es el flujo que describe el resto de esta guía: entrar a `examples/showCases/<ejemplo>/build` (o carpeta equivalente) y ejecutar `make` directamente contra el Makefile de ese ejemplo. No necesitas CMake para esto — el Makefile ya trae las banderas de compilación resueltas.

- Compilador C++ y MPI (ver tabla A y C).

### C. Solo si compilas con el flujo CMake (librería completa o proyecto propio vía `CMakeLists.txt`)

| Requisito | Propósito | Comando (Ubuntu/Debian) |
|---|---|---|
| **CMake** | Genera los archivos de compilación a partir de `CMakeLists.txt`. | `sudo apt install cmake` |
| **Clang** | Alternativa/complemento a GCC; algunos `showCases` y el pipeline de CI del proyecto compilan con ambos para detectar problemas de portabilidad. | `sudo apt install clang` |
| **clang-format** | Solo si vas a contribuir código y quieres respetar el estilo del proyecto; no afecta la compilación. | `sudo apt install clang-format` |

### D. Opcionales / recomendados (afectan funcionalidad, no la compilación mínima)

| Paquete | Para qué sirve | Comando |
|---|---|---|
| **OpenMPI** (`openmpi-bin`, `libopenmpi-dev`) | Ejecución paralela. La mayoría de los `showCases` lo usan por defecto, así que en la práctica es casi obligatorio si vas a correr ejemplos tal cual vienen, pero la librería puede compilarse y correr en un solo núcleo sin él. | `sudo apt install openmpi-bin libopenmpi-dev` |
| **HDF5** (`libhdf5-dev`, `libhdf5-mpi-dev`) | Habilita la salida en formato HDF5 en los ejemplos que lo soportan. | `sudo apt install libhdf5-dev libhdf5-mpi-dev` |
| **ImageMagick** | Algunos ejemplos generan `.gif` invocando `convert`. | `sudo apt install imagemagick` |
| **ccache** | Acelera recompilaciones repetidas; útil en desarrollo iterativo. | `sudo apt install ccache` |

> **Nota sobre `libtbb-dev`:** no aparece como requisito en la documentación oficial del proyecto (ni en el flujo Makefile ni en el CMake estándar). Es posible que tu versión o el ejemplo específico que estás siguiendo (por ejemplo, algún acoplamiento GPU o con otra librería) sí lo necesite, pero antes de dejarlo como requisito general en el README conviene confirmar contra el `CMakeLists.txt` o el Makefile del ejemplo puntual que lo pide — de lo contrario alguien podría instalarlo sin necesitarlo, o al revés, omitirlo pensando que es solo para formato de código (fácil de confundir con `clang-format`).

Con esta separación, la tabla original de la sección 1 del README puede quedar solo con lo que es verdaderamente necesario para *cualquier* instalación (tabla A), y remitir a las tablas C y D según el flujo que la persona vaya a seguir.

---
   
