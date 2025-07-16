## Simulación simple de metano en LAMMPS con salida de trayectoria

Esta simulación realiza una **dinámica molecular** sencilla de **moléculas de metano (CH₄)** utilizando LAMMPS, enfocándose en describir su evolución a temperatura ambiente (300 K) y generando una **salida de trayectoria en formato XYZ** para análisis visual o posterior procesamiento. A continuación, se explica las componentes clave y el objetivo de cada instrucción:

### Propósito general

- Simular el comportamiento molecular de metano a temperatura ambiente, usando el modelo Lennard-Jones para describir las interacciones entre átomos de carbono e hidrógeno.
- Obtener trayectorias de las posiciones atómicas durante la simulación para visualización o análisis estructural.


### Desglose y explicación del script

- **units real:** Define unidades físicas tipo “real”, como las usadas en química: distancias en *Å*, energías en *kcal/mol*, tiempos en *fs*, útil para moléculas orgánicas.
- **atom_style full:** Indica que se modelarán moléculas completas, permitiendo potenciales tipo Coulomb, ángulos y enlaces (aunque aquí sólo se usa Lennard-Jones, el formato permite versatilidad).
- **read_data methane.data:** Carga la estructura inicial, tipos atómicos y sus posiciones desde un archivo externo, que describe la geometría y constitución de una o más moléculas de metano.
- **pair_style lj/cut 10.0:** Selecciona el potencial de Lennard-Jones con un *cutoff* de 10 Å para describir la interacción entre todos los átomos.
- **pair_coeff:** Especifica los parámetros Lennard-Jones para cada combinación de tipos atómicos:
    - 1-1 (C-C): ε=0.1094, σ=3.40
    - 1-2 (C-H): ε=0.0300, σ=2.50
    - 2-2 (H-H): ε=0.0157, σ=2.42
- **neighbor / neigh_modify:** Define el rango y la frecuencia de actualización de la lista de vecinos, optimizando el cálculo de interacciones, importante en sistemas con potencial truncado.
- **velocity all create 300.0 12345:** Inicializa las velocidades siguiendo una distribución Maxwell-Boltzmann a 300 K, con una semilla aleatoria para reproducibilidad.
- **fix 1 all nve:** Usa el integrador NVE, es decir, simula un sistema aislado con Energía, Volumen y Número de partículas constantes.
- **thermo / thermo_style:** Imprime cada 10 pasos propiedades termodinámicas clave, como:
    - paso (step), temperatura (temp), energía potencial (pe), energía total (etotal) y presión (press).
- **dump / dump_modify:** Genera un archivo de **trayectoria** llamado `traj_methane.xyz` cada 10 pasos, en formato XYZ, con asignación explícita de elementos “C” (tipo 1) y “H” (tipo 2).
- **minimize:** Realiza una **minimización de energía** antes de iniciar la dinámica molecular para evitar fuerzas iniciales muy altas debidas a posiciones iniciales no optimizadas.
- **timestep 1.0:** Establece un paso de tiempo de 1 fs.
- **run 100:** Simula la evolución del sistema durante 100 pasos (100 fs), un periodo muy corto, pero suficiente para probar la operación y obtener una pequeña muestra de trayectoria.


### Características y aplicación

- **Simplicidad:** Este es un ejemplo básico didáctico, útil para aprender el flujo de trabajo en LAMMPS y generar una primera *trayectoria molecular*.
- **Visualización:** El archivo XYZ producido puede visualizarse en programas como VMD, Ovito o PyMOL para analizar la evolución de la molécula.
- **Flexibilidad:** Aunque aquí se usa sólo Lennard-Jones, la preparación del formato “full” permite en el futuro agregar otros potenciales (enlaces, ángulos, cargas, etc.).
- **Mínimo computacional:** Simulación diseñada para ser ejecutada rápidamente en cualquier PC estándar como prueba o punto de partida pedagógico.

Se puede adaptar este ejemplo para un sistema más grande, incorporar enlaces explícitos, o extender el tiempo de simulación, solo es necesario que modificar parámetros como el número de pasos, la selección del potencial, o el formato de salida.
