## Simulación en Lammps usando el potencial AIREBO

Esta simulación utiliza el potencial **AIREBO** para modelar **polietileno (polyethylene)**. 

### Objetivo de la Simulación

Esta simulación tiene como objetivo evaluar el comportamiento de un sistema de **polietileno** modelado con el potencial **AIREBO**, que es adecuado para hidrocarburos (átomos de carbono e hidrógeno). Es una simulación de referencia o *benchmark* para validar rendimiento o precisión.

### Detalle del Input Script

```lammps
# AIREBO polyethelene benchmark

units               metal
```

- Establece el sistema de unidades. `"metal"` implica unidades típicas en física de materiales:
    - Distancia en Ångstroms (Å)
    - Tiempo en picosegundos (ps)
    - Energía en eV
    - Temperatura en Kelvin

```lammps
atom_style          atomic
```

- Se utiliza el estilo `"atomic"`: adecuado para modelos sin enlaces explícitos, cargas parciales u otros atributos atómicos complejos.

```lammps
read_data           data.airebo
```

- Carga una estructura inicial del archivo `data.airebo`, que contiene:
    - Número y tipo de átomos
    - Posiciones atómicas
    - Tipos atómicos (probablemente C y H)
    - Información de la celda de simulación

```lammps
replicate           17 16 2
```

- Replica la celda leída en las direcciones **x (17 veces), y (16 veces), z (2 veces)**, incrementando significativamente el tamaño del sistema. Proporciona un sistema más representativo para simulación a granel.

```lammps
neighbor            0.5 bin
neigh_modify        delay 5 every 1
```

- Establece parámetros para la lista de vecinos, usada en el cálculo de interacciones interatómicas. Un skin de `0.5 Å` y actualizaciones cada `5` pasos, verificando cada `1` paso.

```lammps
pair_style          airebo 3.0 1 1
pair_coeff          * * CH.airebo C H
```

- Se utiliza el **potencial AIREBO**, apropiado para hidrocarburos. Sus parámetros:
    - Corte de 3.0 Å,
    - Primer `1` activa la parte REBO (enlaces covalentes),
    - Segundo `1` activa las interacciones de van der Waals tipo LJ.
- El archivo de potencial `CH.airebo` contiene los parámetros del modelo.
- `pair_coeff * * CH.airebo C H`: asigna elementos a los tipos atómicos.

```lammps
velocity            all create 300.0 761341
```

- Inicializa velocidades (distribuidas mediante Maxwell-Boltzmann) para una temperatura de **300 K** con semilla aleatoria `761341`.

```lammps
fix                 1 all nve
```

- Integra las ecuaciones de movimiento con el método **NVE** (Número de partículas, Volumen, y Energía constantes). No hay control de temperatura o presión directamente.

```lammps
timestep            0.00025
```

- El paso de tiempo es de **0.25 femtosegundos (fs)**.

```lammps
thermo              10
```

- Imprime el resumen termodinámico (energía, temperatura, etc.) cada **10 pasos**.

```lammps
run                 100
```

- Ejecuta la simulación por **100 pasos**. Aunque inicialmente es corto (50 fs totales), se ajusta a una prueba de rendimiento (benchmark) o prueba inicial del sistema.


### Interpretación

La simulación inicializa un sistema grande de polietileno y simula sus primeros instantes con un potencial que permite modelar enlaces, interacciones débiles, y posibles reacciones. Este benchmark ideal para evaluar el desempeño del potencial AIREBO y verificar la estabilidad del sistema.

### Perspectivas

- Para simulaciones físicas significativas, se recomienda extender considerablemente el número de pasos (`run 100000` o más).
- Se puede usar un *thermostat* (como `fix nvt` o `fix langevin`) para controlar la temperatura.
- Se podrían agregar análisis estructurales si se desean evaluar propiedades del polímero.

## Requisitos para estudiar propiedades mecánicas

La ejecución actual es un benchmark corto: 100 pasos en NVE no permiten obtener propiedades de equilibrio ni una curva mecánica estadísticamente confiable. Para convertir este caso en una simulación completa de polietileno se requiere el siguiente flujo.

### 1. Minimización y equilibrio

Antes de medir propiedades, hay que eliminar contactos desfavorables y equilibrar el sistema:

```lammps
minimize            1.0e-8 1.0e-10 10000 100000
velocity            all create 300.0 761341 mom yes rot yes dist gaussian
fix                 1 all nvt temp 300.0 300.0 0.1
thermo              1000
run                 200000
unfix               1
```

Con `timestep 0.00025`, esos 200 000 pasos corresponden a 50 ps. Para propiedades de equilibrio más robustas conviene usar una etapa de producción de al menos 0.5-2 ns, comprobando que temperatura, energía, densidad y presión hayan alcanzado una meseta. Si se busca la densidad experimental, puede añadirse una etapa NPT antes de la producción; para sólidos o cajas anisotrópicas debe revisarse cuidadosamente qué dimensiones se permiten variar.

El paso de tiempo se cambió de `0.0005` a `0.00025 ps` (0.25 fs). El valor anterior no es necesariamente inválido, pero este paso más conservador representa mejor las vibraciones rápidas de los enlaces C-H y reduce el riesgo de inestabilidad durante la minimización, el equilibrio y la deformación. Debe verificarse la estabilidad observando la energía total y la ausencia de advertencias de LAMMPS.

### 2. Medición de propiedades de equilibrio

Durante una corrida NVT o NVE de producción se deben registrar, como mínimo, temperatura, energía potencial, energía total, presión, volumen y densidad. Las medias y desviaciones estándar deben calcularse descartando el transitorio inicial. Para obtener propiedades elásticas también se recomienda evaluar la matriz de constantes elásticas mediante fluctuaciones o pequeñas deformaciones.

### 3. Ensayos mecánicos

Se necesitan simulaciones separadas y réplicas independientes para:

- tracción uniaxial en las direcciones `x`, `y` y `z`;
- compresión, si el objetivo incluye respuesta compresiva;
- cizallamiento, para obtener módulos de corte;
- relajación lateral durante la tracción, para estimar el coeficiente de Poisson.

Cada ensayo debe partir de una configuración equilibrada nueva. La deformación debe aplicarse lentamente, por ejemplo hasta 5-10 %, usando deformaciones pequeñas para la región elástica y una velocidad de deformación reportada explícitamente. En LAMMPS, una configuración típica debe combinar `fix deform` con un termostato compatible y guardar la tensión, la deformación y el volumen:

```lammps
variable            Lx0 equal ${lx}
variable            strain equal (lx-v_Lx0)/v_Lx0
fix                 2 all nvt temp 300.0 300.0 0.1
fix                 3 all deform 1 x erate 1.0e-6 units box remap x
thermo              1000
thermo_style        custom step temp pe ke etotal press pxx pyy pzz lx ly lz v_strain
dump                1 all custom 1000 tensile_x.dump id type x y z
run                 200000
```

La velocidad y la duración deben ajustarse al tamaño del sistema y al coste computacional. Para un ensayo de tracción, el módulo de Young se obtiene de la pendiente inicial de `sigma_xx` frente a `epsilon_x`; no debe confundirse la presión hidrostática `Press` con una tensión normal. Las tensiones deben convertirse a las unidades deseadas y corregirse por el área o volumen instantáneo según la convención utilizada.

### 4. Análisis y validación

Para que el resultado sea una propiedad mecánica y no una sola trayectoria, se requiere:

- varias semillas aleatorias y, preferiblemente, varias configuraciones iniciales;
- estimación de incertidumbre mediante promedios por bloques o réplicas;
- comprobación de convergencia respecto al tamaño de la caja, duración de equilibrio y velocidad de deformación;
- comparación con densidad, módulo de Young, coeficiente de Poisson y resistencia experimental o bibliográfica;
- inspección de la estructura, orientación de cadenas, distribución de enlaces y temperatura durante cada ensayo;
- guardar versiones de LAMMPS, archivo `CH.airebo`, parámetros, semillas y condiciones de frontera para reproducibilidad.

Los scripts `01-temperature.py` y `2-more-variables.py` sirven para una inspección inicial del log. Para un estudio completo deben ampliarse para calcular promedios, fluctuaciones, pendientes de las curvas tensión-deformación y barras de error. Ambos deben ejecutarse desde este directorio, donde el nombre correcto del archivo de salida es `log.airebo`.
