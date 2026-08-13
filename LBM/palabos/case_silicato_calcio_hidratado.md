
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

La salida en formato VTK se abre en ParaView, igual que en los demás ejemplos del README. Las magnitudes de interés para este caso de estudio son:

- **Saturación de agua** dentro del espacio poroso, para cuantificar aire atrapado cerca de las inclusiones de PET.
- **Ángulo de contacto local** en la interfaz agua-PET-poro, comparado contra el valor de calibración del Paso 4.
- **Permeabilidad efectiva** de la matriz C-S-H con y sin inclusiones, para relacionar la mojabilidad con el desempeño de transporte de humedad del material compuesto.

---

*Nota general:* los fragmentos de código anteriores son un esqueleto ilustrativo del flujo de trabajo (definición de dominio → geometría → segunda fase → interacción → condiciones de frontera → postproceso), inspirado en la estructura real de `examples/showCases/multiComponent2d` de Palabos. Antes de compilarlos tal cual, conviene contrastarlos contra el ejemplo `rayleighTaylor2D.cpp` de esa carpeta (o el que corresponda en tu versión instalada), ya que la API exacta de inicialización de campos y acoplamientos Shan-Chen puede variar entre versiones de Palabos.
