
# Flujo Multifásico — Silicato de Calcio Hidratado (C-S-H) con inclusiones de PET

## Contexto y motivación

Una línea de investigación activa en materiales de construcción sostenibles es reemplazar parte del agregado mineral del concreto por partículas de **PET reciclado** (tereftalato de polietileno), reduciendo residuo plástico sin sacrificar demasiado desempeño mecánico. Un punto crítico para que esto funcione es entender cómo se comporta el **agua de poro** (la solución intersticial que satura la matriz de **silicato de calcio hidratado, C-S-H**, el principal producto de hidratación del cemento) cuando encuentra inclusiones de un material no mojante como el PET: la mojabilidad deficiente del PET frente al agua puede dejar vacíos de aire atrapados en la interfaz pasta-plástico, que luego son puntos de arranque de fisuras.

El método Lattice Boltzmann, y en particular el **modelo Shan-Chen multicomponente**, es una herramienta natural para esto: permite representar dos fases inmiscibles (agua de poro y una fase "no mojante" que hace de proxy del PET) e imponer un ángulo de contacto específico ajustando la fuerza de interacción fluido-fluido y fluido-sólido, sin necesidad de rastrear explícitamente la interfaz.

Este ejemplo reutiliza la misma familia de ejemplos que ya trae Palabos en `examples/showCases/multiComponent2d` (basados en el descriptor `ForcedShanChenD2Q9Descriptor`), adaptándola al caso de estudio.

### Paso 1 — Definir el dominio y el modelo

Se usa una malla D2Q9 con el descriptor de Shan-Chen forzado, que añade al lattice estándar los campos necesarios para el término de fuerza de interacción entre componentes:

```cpp
#include "palabos2D.h"
#include "palabos2D.hh"
#include <cstdlib>
#include <iostream>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::ForcedShanChenD2Q9Descriptor

// Dimensiones del dominio (en unidades de lattice)
const plint nx = 400;   // dirección x
const plint ny = 200;   // dirección y
```

### Paso 2 — Generar la geometría porosa de C-S-H

La matriz de C-S-H no es un medio homogéneo: es un sólido con porosidad interconectada. Para un primer ejemplo pedagógico, esa porosidad se puede generar de forma sintética (por ejemplo, con un campo aleatorio umbralizado) en vez de importar una micro-CT real, dejando la puerta abierta a reemplazar esta función más adelante por datos experimentales:

```cpp
// Máscara binaria: true = poro (fluido), false = sólido C-S-H
std::vector<std::vector<bool>> generatePorousCSH(plint nx, plint ny, T porosity, plint seed) {
    std::vector<std::vector<bool>> mask(nx, std::vector<bool>(ny, false));
    srand(seed);
    for (plint i = 0; i < nx; ++i) {
        for (plint j = 0; j < ny; ++j) {
            mask[i][j] = (static_cast<T>(rand()) / RAND_MAX) < porosity;
        }
    }
    return mask;
}
```

### Paso 3 — Insertar las inclusiones de PET como segunda fase

En vez de tratar el PET como un obstáculo sólido rígido, se modela como la **segunda componente fluida** del par Shan-Chen (de alta densidad relativa y no mojante), lo que permite que el propio modelo resuelva la forma de la interfaz agua/PET dentro de los poros, en vez de imponerla a mano:

```cpp
// Coloca inclusiones circulares de "fase PET" dentro del espacio poroso
void placePETInclusions(MultiScalarField2D<T>& petPhase,
                         std::vector<Array<T,2>> const& centers,
                         T radius)
{
    for (plint i = 0; i < petPhase.getNx(); ++i) {
        for (plint j = 0; j < petPhase.getNy(); ++j) {
            for (auto const& c : centers) {
                T dx = i - c[0];
                T dy = j - c[1];
                if (dx*dx + dy*dy < radius*radius) {
                    petPhase.get(i, j) = (T) 1.9;  // densidad inicial alta -> fase PET
                }
            }
        }
    }
}
```

### Paso 4 — Configurar la interacción Shan-Chen (mojabilidad)

El parámetro `G` de interacción fluido-fluido controla la tensión interfacial, y un segundo término de interacción fluido-sólido controla el ángulo de contacto contra las paredes de C-S-H. Valores típicos de partida (a calibrar contra el ángulo de contacto experimental agua-PET, ~70-80°):

```cpp
T G_fluidFluid  = -1.2;   // controla tensión superficial agua/PET
T G_fluidSolid  = -0.4;   // controla mojabilidad contra la matriz C-S-H (ajustar para no-mojante)
```

> Estos valores son un punto de partida razonable, no una calibración final: antes de sacar conclusiones físicas hay que ajustar `G_fluidSolid` iterativamente hasta reproducir el ángulo de contacto medido experimentalmente para el par agua-PET, siguiendo el procedimiento estándar de calibración del modelo Shan-Chen (por ejemplo, gota sobre superficie plana como caso de referencia antes de correr la geometría porosa completa).

### Paso 5 — Condiciones de frontera y ejecución

Se imponen condiciones periódicas o de flujo impuesto en los bordes exteriores del dominio, según si el objetivo es estudiar imbibición espontánea (periódico) o flujo forzado (gradiente de presión), y luego se corre la simulación con el mismo patrón de los demás ejemplos del README:

```bash
cd examples/showCases/multiComponent2d/build
make
mpirun -np 4 ./cshPetMultiphase
```

### Paso 6 — Postprocesamiento

El postprocesamiento debe comenzar comprobando qué variables escribe realmente la aplicación. No basta con abrir una imagen de la interfaz: para este caso hay que conservar los campos de fase, la geometría sólida, la velocidad y la presión, además del instante temporal y la conversión entre unidades lattice y físicas.

#### 6.1 Verificar los archivos de salida

Después de ejecutar el caso, localizar las salidas:

```bash
find . -type f \( -name '*.vtk' -o -name '*.vti' -o -name '*.vtu' \
    -o -name '*.pvtu' -o -name '*.pvd' \) -print | sort
```

Abrir en ParaView el archivo `.pvd` si existe, porque normalmente agrupa los instantes temporales. Para una ejecución MPI, abrir el archivo maestro `.pvtu`; no abrir una pieza individual. Si no se generan archivos, revisar las llamadas de escritura del código y habilitar explícitamente las variables necesarias.

#### 6.2 Visualizar la geometría y las fases

En ParaView:

1. Seleccionar **File > Open** y abrir la colección temporal o el archivo de salida.
2. Pulsar **Apply**.
3. Mostrar la máscara sólida de C-S-H con **Threshold** o **Contour** si fue exportada.
4. Seleccionar el campo indicador de agua o PET en **Color By**.
5. Aplicar **Slice** para observar el interior de la matriz porosa.
6. Usar **Contour** para aproximar la interfaz entre agua y la segunda fase.
7. Activar la leyenda y documentar si el campo es `Point Data` o `Cell Data`.

La variable de fase debe estar definida de forma inequívoca. Por ejemplo, debe documentarse qué intervalo representa agua, PET y sólido, y si `petPhase = 1.9` es una densidad inicial o una fracción de volumen. No debe llamarse directamente "PET" a una segunda componente fluida sin aclarar que es un proxy mesoscópico de la inclusión.

#### 6.3 Calcular la saturación de agua

La saturación de agua se calcula sobre el volumen poroso, excluyendo la matriz sólida y las celdas que correspondan a la segunda fase según la definición del modelo:

$$
S_w = \frac{V_{agua}}{V_{poros}}
$$

En ParaView se puede obtener mediante **Threshold** sobre el indicador de agua, seguido de **Integrate Variables**, si el campo y el criterio de fase lo permiten. Si el campo es una fracción continua, se debe integrar esa fracción en lugar de contar celdas. Exportar el resultado a CSV para cada instante y graficar $S_w(t)$.

Documentar el umbral utilizado, el volumen de la matriz excluido y la definición de poro. La saturación debe compararse entre el sistema sin inclusiones y el sistema con PET, usando la misma geometría y el mismo criterio.

#### 6.4 Medir el ángulo de contacto

El ángulo de contacto no debe leerse directamente de una imagen sin calibración. El procedimiento recomendado es:

1. Ejecutar primero una gota de agua sobre una superficie plana con los mismos parámetros de interacción.
2. Ajustar `G_fluidSolid` hasta reproducir el ángulo experimental agua-PET dentro de la incertidumbre elegida.
3. Guardar la relación entre `G_fluidSolid` y el ángulo medido.
4. Ejecutar la geometría porosa de C-S-H/PET con el parámetro calibrado.
5. Extraer un corte local de la interfaz y ajustar la geometría de la interfase lejos de la región de contacto.
6. Reportar el ángulo, la posición, el instante y el método de ajuste.

La presencia de porosidad, resolución lattice y curvatura local puede hacer que el ángulo aparente difiera del ángulo de la superficie plana. Por ello debe informarse como ángulo aparente local y no como una propiedad universal del material.

#### 6.5 Calcular la permeabilidad efectiva

Para estudiar transporte de humedad se necesitan dos simulaciones comparables: C-S-H sin inclusiones y C-S-H con inclusiones de PET. En régimen estacionario, medir el caudal volumétrico $Q$, la longitud del dominio $L$, el área transversal $A$ y la caída de presión $\Delta p$.

La permeabilidad puede estimarse con la ley de Darcy:

$$
K = \frac{Q\,\mu\,L}{A\,\Delta p}
$$

donde $\mu$ es la viscosidad dinámica en unidades físicas. Antes de convertir el resultado, registrar la relación entre unidades lattice y unidades físicas para longitud, tiempo, velocidad, presión y viscosidad. Verificar que el caudal y la caída de presión son aproximadamente constantes durante el intervalo usado.

En ParaView, utilizar **Calculator** para obtener la componente de velocidad relevante, **Slice** o una superficie de integración para calcular el flujo y **Integrate Variables** cuando la malla lo permita. Exportar $Q$, $\Delta p$ y $K$ junto con el tiempo y la resolución de la malla.

#### 6.6 Validación y reproducibilidad

Antes de interpretar los resultados:

- comprobar que la simulación alcanza un régimen estacionario o justificar el análisis transitorio;
- repetir el caso con al menos dos resoluciones lattice;
- comprobar la sensibilidad a `G_fluidFluid` y `G_fluidSolid`;
- comparar la calibración del ángulo de contacto con el caso de gota sobre superficie plana;
- verificar conservación de masa de cada componente;
- comparar permeabilidad y saturación con y sin inclusiones bajo las mismas condiciones;
- guardar parámetros, semilla, geometría, versión de Palabos y archivos de salida.

El resultado mínimo del caso debe ser una tabla con $S_w(t)$, ángulo de contacto calibrado y aparente, $Q$, $\Delta p$ y $K$, acompañada de las conversiones de unidades y de la incertidumbre numérica.

---

*Nota general:* los fragmentos de código anteriores son un esqueleto ilustrativo del flujo de trabajo (definición de dominio → geometría → segunda fase → interacción → condiciones de frontera → postproceso), inspirado en la estructura real de `examples/showCases/multiComponent2d` de Palabos. Antes de compilarlos tal cual, conviene contrastarlos contra el ejemplo `rayleighTaylor2D.cpp` de esa carpeta (o el que corresponda en tu versión instalada), ya que la API exacta de inicialización de campos y acoplamientos Shan-Chen puede variar entre versiones de Palabos.
