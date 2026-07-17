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
timestep            0.0005
```

- El paso de tiempo es de **0.5 femtosegundos (fs)**.

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
