# Palabos

Palabos es una **librería de código abierto (open-source)** escrita en **C++** que se utiliza para desarrollar y ejecutar simulaciones de fluidos basadas en el **Método de Lattice Boltzmann (LBM)**, del cual hablamos anteriormente.

El nombre "Palabos" es un acrónimo de **PArLAs-BOlS-Solver** (Parallel Lattice Boltzmann Solver).

## Características principales

- **Método Base:** Se basa exclusivamente en el Método Lattice Boltzmann (LBM).
- **Lenguaje:** Es una librería escrita en C++, lo que permite a los usuarios escribir sus propias aplicaciones de CFD aprovechando las estructuras y optimizaciones que ofrece Palabos.
- **Paralelización:** Está diseñado desde su base para ser **altamente escalable y eficiente** en sistemas paralelos, utilizando la Interfaz de Paso de Mensajes (MPI) para ejecutarse en clústeres y supercomputadoras. Esto es clave dada la naturaleza local del LBM.
- **Plataforma:** Es un software multiplataforma, compatible con Linux, Windows y macOS.

### Desarrollo y soporte

El desarrollo y mantenimiento de Palabos están supervisados principalmente por la empresa **FlowKit Ltd.** y tienen una fuerte colaboración con la **Universidad de Ginebra** en Suiza, que es pionera en la investigación del método Boltzmann.

## Aplicaciones típicas

Palabos es valorado por su capacidad para manejar simulaciones complejas, donde el LBM demuestra superioridad sobre los métodos CFD tradicionales. Sus principales áreas de aplicación incluyen:

1. **Flujos Multifásicos:** Simulación de la interacción entre diferentes fases de fluidos, como la mezcla de aceite y agua o la dinámica de burbujas y gotas.
2. **Medios Porosos:** Modelado del flujo de fluidos a través de estructuras porosas complejas (p. ej., filtración, pilas de combustible).
3. **Flujos Térmicos:** Simulación de convección natural y forzada, incluyendo la transferencia de calor.
4. **Flujos con Reacciones Químicas:** Permite integrar modelos de reacción dentro de la dinámica del fluido.

En resumen, Palabos es una herramienta esencial para investigadores y profesionales que buscan implementar de forma eficiente simulaciones de fluidos de vanguardia utilizando el enfoque mesoscópico del Lattice Boltzmann.

## Instalación, configuración y testeo

Instalar, configurar y testear **Palabos** en Linux generalmente se realiza compilando el código fuente. Palabos depende de un conjunto de bibliotecas de terceros, siendo **MPI (Message Passing Interface)** crucial para la ejecución paralela.

Aquí se tiene una guía paso a paso para la instalación y configuración básica, seguida de cómo ejecutar el test de la librería.

### 1. Requisitos previos

Antes de instalar Palabos, asegúrense de tener las siguientes herramientas y librerías instaladas en su distribución Linux (p. ej., Ubuntu, Debian, Fedora):

| **Requisito**           | **Propósito**                                      | **Comando de Instalación (Ej. Ubuntu/Debian)**        |
| ----------------------- | -------------------------------------------------- | ----------------------------------------------------- |
| **Compilador C++**      | Necesario para compilar el código fuente.          | `sudo apt update && sudo apt install build-essential` |
| **MPI (OpenMPI/MPICH)** | Esencial para la ejecución paralela (obligatorio). | `sudo apt install openmpi-bin libopenmpi-dev`         |
| **Git**                 | Para clonar el repositorio de Palabos.             | `sudo apt install git`                                |

------

### 2. Instalación de Palabos

#### 2.1 Obtener el código fuente

Utilizar Git para clonar el repositorio de Palabos (o puede descargar el archivo ZIP desde su sitio oficial):

Bash

```
git clone https://gitlab.com/unigespc/palabos.git
```

Esto creará una carpeta llamada `palabos` en tu directorio actual.

#### 2.2 Instalar prerrequisitos

Para que palabos funcione correctamente se requieren los siguientes paquetes.

##### Distribuciones derivadas de Debian

Bash
```
$ sudo apt install gcc clang clang-format cmake make libtbb-dev
```
#### 2.3 Compilar la biblioteca

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

### 3. Testeo y ejecución

Una vez construido el ejecutable, se puede probar así:

#### 3.1 Ejecución en un solo núcleo

Para ejecutar la simulación en un solo núcleo de CPU (secuencial):

Bash

```
./laminarChannel
```

#### 3.2 Ejecución paralela con MPI

Para probar la funcionalidad paralela (el corazón de Palabos), utiliza el comando `mpirun` o `mpiexec`. Aquí se ejecutará el test en 4 núcleos (cambiar el número según tu CPU):

Bash

```
mpirun -np 4 ./laminarChannel
```

- **`-np 4`**: Indicar a MPI que lance el programa y distribuya la carga de trabajo entre **4** procesos (núcleos).

#### 3.3 Verificación de resultados

Una ejecución exitosa debe terminar sin errores de segmentación ni errores de MPI y debe generar mensajes de progreso o un resumen final. La generación de archivos de salida debe verificarse por separado, porque depende del ejemplo y de la versión de Palabos.

Desde el directorio del ejemplo, localizar los resultados:

```bash
find . -type f \( -name '*.vtk' -o -name '*.vti' -o -name '*.vtu' \
   -o -name '*.pvtu' -o -name '*.pvd' -o -name '*.h5' \) -print | sort
```

Los archivos VTK, VTU, PVTU o PVD pueden abrirse con ParaView. En una ejecución MPI se debe abrir preferentemente el archivo maestro (`.pvtu` o `.pvd`) y no una pieza individual.

Si no se generan resultados, revisar el código fuente del ejemplo, el `Makefile` y las llamadas de escritura. No todos los casos escriben archivos automáticamente.

## 4. Caso de estudio: flujo laminar en un canal

El ejemplo `laminarChannel` se utiliza aquí como caso de estudio para verificar el ciclo completo de Palabos: compilar una aplicación, ejecutarla, inspeccionar sus archivos y extraer magnitudes físicas. Es un caso de validación de flujo, no un modelo de propiedades mecánicas de un sólido.

### 4.1 Objetivo

El objetivo es comprobar que el campo de velocidad y las variables hidrodinámicas evolucionan de forma razonable en un canal laminar. La comparación debe hacerse con la documentación y los parámetros de la versión descargada, ya que la geometría, el número de pasos, el caudal y los nombres de salida pueden cambiar entre releases.

### 4.2 Ejecutar el caso

Después de compilar desde `examples/showCases/laminarChannel/build`, ejecutar primero en serie:

```bash
./laminarChannel
```

Guardar la salida de consola para revisar errores y tiempos:

```bash
./laminarChannel 2>&1 | tee laminarChannel.serial.log
```

Cuando el caso secuencial funcione, probar MPI:

```bash
mpirun -np 4 ./laminarChannel 2>&1 | tee laminarChannel.mpi.log
```

No comparar dos ejecuciones como si fueran idénticas si utilizan distinto número de procesos sin comprobar que el caso y las condiciones iniciales son reproducibles.

### 4.3 Identificar los archivos de salida

Al finalizar cada ejecución, revisar los archivos generados:

```bash
find . -type f \( -name '*.vtk' -o -name '*.vti' -o -name '*.vtu' \
   -o -name '*.pvtu' -o -name '*.pvd' \) -printf '%TY-%Tm-%Td %TH:%TM %10s %p\n' | sort
```

Interpretar los archivos según su función:

- `.vtk`, `.vti` o `.vtu`: datos de un instante o una malla concreta.
- `.pvtu`: colección de piezas producidas por una ejecución paralela.
- `.pvd`: colección temporal que permite recorrer varios instantes.

Los nombres exactos deben confirmarse en la carpeta del ejemplo. Si solo hay archivos de geometría o no hay resultados, la aplicación puede no tener habilitada la escritura o puede requerir una configuración adicional.

### 4.4 Abrir y visualizar en ParaView

Instalar ParaView si es necesario:

```bash
sudo apt install -y paraview
paraview
```

Seguir este flujo:

1. Seleccionar **File > Open**.
2. Abrir el `.pvd` si existe; así se carga la secuencia temporal completa.
3. Si el caso fue paralelo, abrir el `.pvtu` maestro.
4. Si solo existe una salida individual, abrir `.vtk`, `.vti` o `.vtu`.
5. Pulsar **Apply**.
6. Usar **Color By** para seleccionar velocidad, presión, densidad u otra variable disponible.
7. Pulsar **Rescale to Data Range** y activar la leyenda desde **View > Color Map Editor**.
8. Usar **Play** para comprobar la evolución temporal.

No seleccionar una variable solo por su nombre: comprobar si está almacenada como `Point Data` o `Cell Data` y documentar qué campo representa.

### 4.5 Extraer un perfil de velocidad

Para analizar cuantitativamente el flujo del canal:

1. Seleccionar el conjunto de datos cargado.
2. Aplicar **Filters > Data Analysis > Plot Over Line**.
3. Definir una línea transversal al canal en una posición donde el flujo esté desarrollado.
4. Pulsar **Apply**.
5. Seleccionar la componente de velocidad disponible.
6. Exportar la tabla mediante **File > Save Data** en CSV.

Guardar junto al CSV las coordenadas de los extremos de la línea, el instante temporal, la variable utilizada y las unidades lattice o físicas. Si el campo de velocidad tiene componentes separadas, indicar cuál se analizó.

### 4.6 Mediciones adicionales

Según las variables exportadas por el ejemplo, pueden utilizarse estos filtros:

- **Calculator:** calcular la magnitud de la velocidad o una expresión derivada.
- **Slice:** inspeccionar un plano longitudinal o transversal.
- **Contour:** localizar superficies de presión o velocidad constante.
- **Stream Tracer:** observar líneas de corriente.
- **Integrate Variables:** obtener integrales sobre una superficie o volumen.
- **Plot Selection Over Time:** seguir una variable en una posición seleccionada.

Si se necesita una magnitud que no aparece en los archivos, debe añadirse su escritura en el código del ejemplo y recompilar. Para repetir el análisis, guardar el estado de ParaView con **File > Save State** en un archivo `.pvsm`.

### 4.7 Validar el caso

Antes de considerar validado el caso `laminarChannel`, comprobar:

- que la ejecución alcanza el número de pasos previsto;
- que los archivos temporales aparecen en el orden esperado;
- que el perfil de velocidad cambia de forma coherente entre el inicio y el régimen estacionario;
- que el perfil y el caudal no dependen de forma inesperada del número de procesos MPI;
- que el refinamiento espacial y el paso temporal son suficientes;
- que las unidades lattice están documentadas antes de convertirlas a unidades físicas;
- que, cuando sea posible, el perfil se compara con una solución analítica o una referencia del ejemplo.

La captura de pantalla demuestra que el resultado puede visualizarse, pero el perfil exportado, la convergencia temporal y la comparación cuantitativa son los elementos que convierten la ejecución en un caso de estudio reproducible.

## 5. Caso de estudio: C-S-H con inclusiones de PET

El documento [case_silicato_calcio_hidratado.md](case_silicato_calcio_hidratado.md) propone un caso de estudio LBM para analizar agua de poro y una segunda fase no mojante dentro de una matriz porosa de silicato de calcio hidratado (C-S-H), usada como proxy mesoscópico de inclusiones de PET.

### 5.1 Qué estudia

El caso combina:

- una geometría porosa sintética o derivada de micro-CT;
- un modelo Shan-Chen multicomponente;
- calibración de la interacción fluido-fluido y fluido-sólido;
- comparación de una matriz C-S-H con y sin inclusiones;
- análisis de saturación de agua, mojabilidad y permeabilidad efectiva.

### 5.2 Estado de implementación

Este caso es una especificación metodológica y un esqueleto de adaptación, no un ejecutable incluido en este repositorio. Los fragmentos de C++ deben integrarse en una versión concreta de `examples/showCases/multiComponent2d` y ajustarse a la API de Palabos descargada. En particular, deben implementarse y verificarse la inicialización de los lattices, el acoplamiento Shan-Chen, las condiciones de frontera, la escritura VTK/PVD y la definición de las unidades.

Por tanto, la ejecución no debe comenzar con una interpretación física de los valores propuestos. Primero debe compilar un caso mínimo de dos fases y comprobarse la estabilidad de la interfaz.

### 5.3 Secuencia recomendada

1. Compilar y ejecutar un ejemplo multicomponente oficial de Palabos.
2. Reproducir una gota sobre una superficie plana para calibrar el ángulo de contacto.
3. Sustituir la superficie plana por una geometría porosa simple.
4. Añadir la máscara de C-S-H y comprobar conservación de masa.
5. Añadir las inclusiones de PET como segunda componente y comparar saturación.
6. Ejecutar los casos con y sin inclusiones bajo las mismas condiciones.
7. Extraer $S_w(t)$, ángulo de contacto, caudal, caída de presión y permeabilidad.
8. Repetir con otra resolución lattice y documentar la sensibilidad a `G_fluidFluid` y `G_fluidSolid`.

### 5.4 Postprocesamiento del caso

El postprocesamiento debe seguir este orden:

1. Abrir el archivo `.pvd` o `.pvtu` en ParaView.
2. Verificar que se exportaron la máscara sólida, el indicador de agua/PET, velocidad y presión.
3. Aplicar **Threshold**, **Slice** y **Contour** para separar matriz, poros e interfaz.
4. Calcular la saturación de agua sobre el volumen poroso, excluyendo el sólido.
5. Ajustar el ángulo de contacto en el caso de calibración y reportar después el ángulo aparente en la geometría porosa.
6. Usar **Calculator** e **Integrate Variables** para obtener caudal y caída de presión.
7. Calcular la permeabilidad mediante la ley de Darcy, manteniendo explícita la conversión de unidades lattice a físicas.
8. Exportar las series temporales y guardar un estado `.pvsm` reproducible.

Las definiciones, ecuaciones, criterios de validación y comandos de inspección están documentados en [case_silicato_calcio_hidratado.md](case_silicato_calcio_hidratado.md).

## 6. Configuración avanzada para proyectos propios

Para empezar un propio proyecto basado en Palabos:

1. **Crear un Directorio de Proyecto:** Crear una nueva carpeta fuera de la estructura de `palabos/examples`.

2. **Copiar el Makefile:** Copiar el `Makefile` de alguno de los ejemplos de Palabos al nuevo directorio de proyecto.

3. **Modificar la fuente:** En el `Makefile` del proyecto, ajustar la variable `SRC_FILES` para apuntar al código fuente C++ personalizado (p. ej., `mySimulation.cpp`).

4. **Inclusión:** Asegurarse que el código C++ incluya la cabecera principal de Palabos:

   C++

   ```
   #include "palabos.h"
   ```

## 7. Notas sobre los flujos de compilación

La guía anterior describe el flujo clásico mediante `Makefile`, que es el más directo para ejecutar los ejemplos de `examples/showCases`. En este flujo no es necesario instalar CMake: el `Makefile` del ejemplo contiene las reglas y rutas de compilación.

Para compilar la biblioteca completa o desarrollar un proyecto mediante `CMakeLists.txt`, instalar además:

```bash
sudo apt install -y cmake clang clang-format
```

Los requisitos opcionales dependen del caso:

| Paquete | Uso | Instalación |
|---|---|---|
| HDF5 | Salida HDF5 en ejemplos compatibles | `sudo apt install libhdf5-dev libhdf5-mpi-dev` |
| ImageMagick | Generación de GIF en algunos ejemplos | `sudo apt install imagemagick` |
| ccache | Acelerar recompilaciones | `sudo apt install ccache` |

`libtbb-dev` no debe tratarse como requisito general. Instalarlo únicamente si el `Makefile`, el `CMakeLists.txt` o el ejemplo concreto lo solicita.

