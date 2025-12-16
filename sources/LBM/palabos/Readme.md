# Palabos

**Palabos** es un nombre que se refiere a un popular y poderoso software en el campo de la Mecánica de Fluidos Computacional (CFD).

## ¿Qué es Palabos?

Palabos es una **librería de código abierto (open-source)** escrita en **C++** que se utiliza para desarrollar y ejecutar simulaciones de fluidos basadas en el **Método Lattice Boltzmann (LBM)**, del cual hablamos anteriormente.

El nombre "Palabos" es un acrónimo de **PArLAs-BOlS-Solver** (Parallel Lattice Boltzmann Solver).

### Características Principales

- **Método Base:** Se basa exclusivamente en el Método Lattice Boltzmann (LBM).
- **Lenguaje:** Es una librería escrita en C++, lo que permite a los usuarios escribir sus propias aplicaciones de CFD aprovechando las estructuras y optimizaciones que ofrece Palabos.
- **Paralelización:** Está diseñado desde su base para ser **altamente escalable y eficiente** en sistemas paralelos, utilizando la Interfaz de Paso de Mensajes (MPI) para ejecutarse en clusters y supercomputadoras. Esto es clave dada la naturaleza local del LBM.
- **Plataforma:** Es un software multiplataforma, compatible con Linux, Windows y macOS.

###  Desarrollo y Soporte

El desarrollo y mantenimiento de Palabos están supervisados principalmente por la empresa **FlowKit Ltd.** y tiene una fuerte colaboración con la **Universidad de Ginebra** en Suiza, que es pionera en la investigación del método Boltzmann.

## Palabos Aplicaciones Típicas de Palabos

Palabos es valorado por su capacidad para manejar simulaciones complejas, donde el LBM demuestra superioridad sobre los métodos CFD tradicionales. Sus principales áreas de aplicación incluyen:

1. **Flujos Multifásicos:** Simulación de la interacción entre diferentes fases de fluidos, como la mezcla de aceite y agua, o la dinámica de burbujas y gotas.
2. **Medios Porosos:** Modelado del flujo de fluidos a través de estructuras porosas complejas (ej. filtración, pilas de combustible).
3. **Flujos Térmicos:** Simulación de convección natural y forzada, incluyendo la transferencia de calor.
4. **Flujos con Reacciones Químicas:** Permite integrar modelos de reacción dentro de la dinámica del fluido.

En resumen, Palabos es una herramienta esencial para investigadores y profesionales que buscan implementar de forma eficiente simulaciones de fluidos de vanguardia utilizando el enfoque mesoscópico del Lattice Boltzmann.

# Instalar, configurar y testear Palabos

Instalar, configurar y testear **Palabos** en Linux generalmente se realiza compilando el código fuente. Palabos depende de un conjunto de bibliotecas de terceros, siendo **MPI (Message Passing Interface)** crucial para la ejecución paralela.

Aquí se tiene una guía paso a paso para la instalación y configuración básica, seguida de cómo ejecutar el test de la librería.

## 1. Requisitos Previos

Antes de instalar Palabos, asegúrarse de tener las siguientes herramientas y librerías instaladas en tu distribución Linux (ej. Ubuntu, Debian, Fedora):

| **Requisito**           | **Propósito**                                      | **Comando de Instalación (Ej. Ubuntu/Debian)**        |
| ----------------------- | -------------------------------------------------- | ----------------------------------------------------- |
| **Compilador C++**      | Necesario para compilar el código fuente.          | `sudo apt update && sudo apt install build-essential` |
| **MPI (OpenMPI/MPICH)** | Esencial para la ejecución paralela (obligatorio). | `sudo apt install openmpi-bin libopenmpi-dev`         |
| **Git**                 | Para clonar el repositorio de Palabos.             | `sudo apt install git`                                |

------

## 2. Instalación de Palabos

### A. Obtener el Código Fuente

Utilizar Git para clonar el repositorio de Palabos (o puede descargar el archivo ZIP desde su sitio oficial):

Bash

```
git clone https://gitlab.com/unige.ch-cfd/palabos.git
```

Esto creará una carpeta llamada `palabos` en tu directorio actual.

### B. Compilación de la Librería

Palabos no requiere un paso de "instalación" tradicional (`make install`); en su lugar, se compila la aplicación del usuario vinculándola a la librería Palabos. La configuración se realiza principalmente a través del **Makefile** incluido.

1. **Navega a la Carpeta:**

   Bash

   ```
   cd palabos
   ```

2. **Configurar el Makefile:**

   - Palabos utiliza una estructura de directorios modular. Para el test inicial, vamos a usar uno de los ejemplos. Navega a un directorio de prueba, por ejemplo:

     Bash

     ```
     cd examples/showCases/laminarChannel
     ```

   - Abrir el `Makefile` con un editor de texto (ej. `nano Makefile` o `vim Makefile`).

   - **Verificar la configuración de MPI:** Asegúrarse que las variables de compilación (como `CXX` y `MPICXX`) apunten a los compiladores de MPI correctos (por defecto suelen ser `mpicxx` o `mpic++`, que ya deberían estar configurados si se insto  `libopenmpi-dev`).

   - **Opcional - Optimización:** Puede ajustar las banderas de optimización si lo desea (ej. `-O3`).

3. Compilar el Ejemplo:

   Ejecutar el comando make para compilar el código fuente del ejemplo:

   Bash

   ```
   make
   ```

   Si la compilación es exitosa, se creará un archivo ejecutable (con el mismo nombre que el directorio, o el nombre definido en el Makefile, en este caso probablemente `laminarChannel`).

## 3. Testeo y Ejecución

Una vez compilado el ejecutable, puedes probarlo para verificar que Palabos y MPI están funcionando correctamente.

### A. Ejecución en un Solo Núcleo

Para ejecutar la simulación en un solo núcleo de CPU (secuencial):

Bash

```
./laminarChannel
```

### B. Ejecución Paralela (Usando MPI)

Para probar la funcionalidad paralela (el corazón de Palabos), utiliza el comando `mpirun` o `mpiexec`. Aquí ejecutaremos el test en 4 núcleos (cambia el número según tu CPU):

Bash

```
mpirun -np 4 ./laminarChannel
```

- **`-np 4`**: Indicar a MPI que lance el programa y distribuya la carga de trabajo entre **4** procesos (núcleos).

### C. Verificación de Resultados

Una ejecución exitosa producirá mensajes de progreso en la consola y, lo más importante, generará archivos de salida para la visualización.

- **Archivos de Salida:** El ejemplo `laminarChannel` generará archivos en formato **VTK** (o similar) que contienen la distribución de velocidad y presión. Estos archivos son legibles por software de post-procesamiento como **ParaView**.

Si la simulación se ejecuta hasta el final sin errores de segmentación o MPI, y los archivos de salida son generados, la instalación de Palabos es exitosa.

## 4. Configuración Avanzada (Para Proyectos Propios)

Para empezar su propio proyecto basado en Palabos:

1. **Crea un Directorio de Proyecto:** Crear una nueva carpeta fuera de la estructura de `palabos/examples`.

2. **Copia el Makefile:** Copiar el `Makefile` de alguno de los ejemplos de Palabos a tu nuevo directorio de proyecto.

3. **Modifica la Fuente:** En el `Makefile` del proyecto, ajustar la variable `SRC_FILES` para apuntar a tu código fuente C++ personalizado (ej. `mySimulation.cpp`).

4. **Inclusión:** Asegúrarse que el código C++ incluya la cabecera principal de Palabos:

   C++

   ```
   #include "palabos.h"
   ```
   
   