

# OpenLB

**OpenLB** es un paquete de software de código abierto (open-source) y una librería en **C++** de alto rendimiento, diseñada específicamente para implementar simulaciones de fluidos y fenómenos de transporte utilizando el **Método Lattice Boltzmann (LBM)**.

OpenLB es una herramienta reconocida en la comunidad de investigación y fue una de las primeras plataformas genéricas de LBM disponibles para la comunidad (licencia GPLv2).

## Características Clave de OpenLB

### 1. Enfoque y Diseño

- **Orientado a Objetos (C++):** El código está escrito en C++ de manera modular y bien legible, lo que facilita tanto a los programadores de aplicaciones como a los desarrolladores avanzados implementar nuevos modelos LBM personalizados o extender las funcionalidades existentes.
- **Kernel LBM Completo:** El núcleo de OpenLB se basa en una amplia variedad de modelos Lattice Boltzmann, lo que permite abordar una gran diversidad de problemas físicos.

### 2. Alto rendimiento y paralelización

Una de las fortalezas más grandes de OpenLB es su arquitectura optimizada para la computación de alto rendimiento:

- **Paralelismo Híbrido:** Es eficiente en plataformas de memoria compartida y distribuida. Soporta:
  - **MPI (Message Passing Interface):** Para paralelismo distribuido (clusters de CPU).
  - **OpenMP:** Para paralelismo de memoria compartida (multihilo en una sola CPU).
  - **CUDA:** Para el uso de **GPUs (unidades de procesamiento gráfico)**, permitiendo simulaciones extremadamente rápidas.
- **Vectorización (SIMD):** Incluye optimizaciones para el procesamiento de datos a nivel de instrucción, mejorando la eficiencia en CPUs.

### 3. Preprocesamiento y Geometría

A diferencia de muchos otros *solvers* de CFD, OpenLB tiene capacidades robustas de preprocesamiento integradas:

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

Asegúrate de tener instalados los siguientes paquetes en tu sistema Linux.

| **Requisito**           | **Propósito**                                          | **Comando de Instalación (Ej. Ubuntu/Debian)**        |
| ----------------------- | ------------------------------------------------------ | ----------------------------------------------------- |
| **Compilador C++**      | Necesario para compilar el código.                     | `sudo apt update && sudo apt install build-essential` |
| **CMake**               | Sistema de construcción recomendado por OpenLB.        | `sudo apt install cmake`                              |
| **MPI (OpenMPI/MPICH)** | Esencial para el paralelismo (distribución de tareas). | `sudo apt install openmpi-bin libopenmpi-dev`         |
| **Git**                 | Para clonar el repositorio de OpenLB.                  | `sudo apt install git`                                |

------

## 2. Empleo de OpenLB

### Obtener el código fuente

Clonar el repositorio oficial de OpenLB usando Git:

Bash

```
git clone https://gitlab.com/openlb/release.git olb-release
```

Esto creará una carpeta llamada `olb-release`.

## 3. Testeo y verificación

Para verificar OpenLB se pueden ejecutar los ejemplos

### Ejecutar un Ejemplo Simple

Para un test funcional en un problema, se puede compilar y ejecutar uno de los ejemplos que vienen en la fuente.

1. **Encuentre un Ejemplo:** Navegar a la carpeta de ejemplos (p. ej. `cd ../examples/laminar/cavity2d/`).

2. **Compilar el Ejemplo:** Dentro de la carpeta del ejemplo, puede usar el comando 

Bash

```
make
```

Si la simulación se ejecuta hasta el final sin errores de MPI o de segmentación y genera archivos de salida (típicamente VTK para post-procesamiento), significa que el ejemplo se ha construido correctamente.

## 4. Postprocessing

El postprocesamiento comienza después de que el ejemplo termina correctamente y genera sus archivos de salida. OpenLB suele escribir resultados en formatos compatibles con ParaView, como VTK, VTU, PVTU o PVD, aunque el formato exacto y el nombre de la carpeta dependen del ejemplo y de la versión descargada.

### 4.1 Verificar las salidas

Desde el directorio del ejemplo, comprobar qué archivos se generaron:

```bash
find . -type f \( -name '*.vtk' -o -name '*.vti' -o -name '*.vtu' \
  -o -name '*.pvtu' -o -name '*.pvd' \) | sort | head -50
```

También conviene revisar el tamaño y la fecha de modificación:

```bash
find . -type f \( -name '*.vtk' -o -name '*.vti' -o -name '*.vtu' \
  -o -name '*.pvtu' -o -name '*.pvd' \) -printf '%TY-%Tm-%Td %TH:%TM %10s %p\n' | sort
```

Si no se genera ningún archivo, revisar la salida de la simulación, el `Makefile` y las llamadas de escritura del código del ejemplo. No todos los casos de OpenLB escriben resultados automáticamente.

### 4.2 Instalar y abrir ParaView

Para instalar ParaView en Ubuntu/Debian:

```bash
sudo apt update
sudo apt install -y paraview
```

Abrir ParaView desde la terminal:

```bash
paraview
```

En la interfaz:

1. Seleccionar **File > Open**.
2. Navegar hasta la carpeta de salida del ejemplo.
3. Abrir el archivo `.pvd` si existe. Este archivo suele agrupar varios instantes temporales.
4. Si el caso se ejecutó con MPI y produjo piezas paralelas, abrir el archivo `.pvtu` o el archivo maestro equivalente, no una pieza individual.
5. Pulsar **Apply** en el panel de propiedades.
6. Presionar **Play** en la barra de animación para recorrer los pasos guardados.

Para una ejecución secuencial, abrir directamente el archivo `.vtk`, `.vti` o `.vtu` correspondiente al instante que se quiera inspeccionar.

### 4.3 Visualizar campos físicos

Una vez cargados los datos, utilizar el panel **Coloring** para seleccionar las variables disponibles, por ejemplo velocidad, presión, densidad o número de Mach. El procedimiento típico es:

1. Seleccionar el campo en **Color By**.
2. Usar **Rescale to Data Range** para ajustar la escala de colores.
3. Cambiar la representación a **Surface**, **Wireframe** o **Volume** según la dimensión y el tipo de dato.
4. Activar una barra de escala desde **View > Color Map Editor**.
5. Añadir una leyenda, título y unidades antes de guardar una figura.

Para obtener magnitudes derivadas, seleccionar **Filters > Alphabetical** y utilizar, según corresponda:

- **Calculator** para crear una expresión a partir de componentes del campo.
- **Gradient** para gradientes o vorticidad cuando el campo y la malla lo permitan.
- **Contour** para superficies de nivel, por ejemplo una presión o fracción de fase constante.
- **Slice** para cortar el dominio y observar perfiles internos.
- **Stream Tracer** para visualizar líneas de corriente a partir del campo de velocidad.
- **Warp By Vector** para deformar visualmente una superficie usando un vector.

Después de añadir cada filtro, pulsar **Apply** y comprobar que la variable seleccionada corresponde a **Point Data** o **Cell Data**. Si una variable no aparece, puede estar almacenada en la otra categoría.

### 4.4 Extraer perfiles y valores numéricos

Para estudiar un perfil en una línea:

1. Seleccionar el conjunto de datos.
2. Aplicar **Filters > Data Analysis > Plot Over Line**.
3. Definir los puntos inicial y final de la línea.
4. Pulsar **Apply**.
5. Exportar la tabla con **File > Save Data** en CSV.

Para obtener valores en puntos, líneas o superficies concretas se pueden utilizar **Probe Location**, **Plot Selection Over Time** y **Integrate Variables**. Estos resultados deben guardarse junto con las coordenadas, el instante temporal, las unidades y la versión del caso.

### 4.5 Exportar figuras y animaciones

Para guardar una imagen de la vista actual:

1. Ajustar cámara, escala de colores, leyenda y unidades.
2. Seleccionar **File > Save Screenshot**.
3. Guardar en PNG o TIFF con una resolución explícita.

Para guardar una evolución temporal:

1. Abrir el archivo temporal `.pvd` o la colección correspondiente.
2. Revisar la animación con **Play**.
3. Seleccionar **File > Save Animation**.
4. Elegir una secuencia de imágenes o un formato de vídeo disponible.

Las figuras científicas deben conservar el nombre del caso, el paso temporal y las unidades. No conviene interpretar una imagen aislada como resultado convergido sin comprobar la evolución temporal.

### 4.6 Automatizar el postprocesamiento

Para repetir el mismo análisis, ParaView permite guardar la configuración mediante **File > Save State** en un archivo `.pvsm`. Después se puede abrir el estado desde la interfaz o ejecutarlo con `pvpython` si la instalación incluye ese comando:

```bash
pvpython --version
pvpython analizar_openlb.py
```

El script `analizar_openlb.py` debe ser creado para el caso concreto y debe indicar explícitamente el archivo de entrada, las variables, los filtros, las unidades y los nombres de salida. Si no existe ese script, realizar el análisis manualmente y guardar el estado `.pvsm` para reproducirlo.

### 4.7 Comprobaciones antes de interpretar resultados

Antes de usar los datos de OpenLB para analizar transporte o propiedades efectivas, comprobar:

- que la simulación terminó sin errores y alcanzó el número de pasos previsto;
- que los archivos contienen todos los instantes temporales esperados;
- que la malla, el campo y las condiciones de contorno corresponden al caso ejecutado;
- que las unidades físicas y la conversión de unidades lattice a físicas están documentadas;
- que el paso temporal y el refinamiento espacial son suficientes;
- que los perfiles y magnitudes integradas son estables durante el intervalo considerado;
- que los resultados se comparan con una solución analítica, una referencia publicada o un estudio de convergencia cuando sea posible.

El postprocesamiento visual sirve para inspeccionar y comunicar los resultados, pero no sustituye la validación numérica ni el análisis de convergencia del modelo LBM.