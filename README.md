# Study of Mechanical Properties

Repositorio de trabajo para estudiar propiedades mecánicas y de transporte en materiales cementicios, polímeros y sus interfaces mediante simulación computacional.

El proyecto reúne ejemplos independientes y una propuesta de integración multiescala. LAMMPS y NWChem contienen casos ejecutables o generadores de estructuras; OpenLB y Palabos contienen guías y casos metodológicos que requieren descargar sus fuentes; el acoplamiento MD-DPD-DEM funciona como hoja de ruta para una implementación posterior.

## Índice de documentación

- [Guía de ejecución y entornos](docs/README.md): Conda/Mamba, Distrobox, Visual Studio Code, dependencias y comandos generales.
- [Ejemplos de NWChem y ASE](nwchem/README.md): moléculas, agua, hielo, clusters de calcio-silicato y tobermorita.
- [Guía de LBM](LBM/Readme.md): introducción a Lattice Boltzmann, OpenLB y Palabos.
- [Guía de OpenLB](LBM/openLB/Readme.md): instalación, compilación, visualización y postprocesamiento con ParaView.
- [Guía de Palabos](LBM/palabos/Readme.md): compilación con Makefile, MPI, caso `laminarChannel` y caso C-S-H/PET.
- [Caso C-S-H/PET](LBM/palabos/case_silicato_calcio_hidratado.md): propuesta multicomponente Shan-Chen y protocolo de postprocesamiento.
- [Plan MD-DPD-DEM](coupling/MD-DPD-DEM/Readme.md): flujo multiescala para composites C-S-H/PE.

## Mapa de tecnologías y proyectos

| Escala o método | Proyecto | Sistema o fenómeno | Estado |
|---|---|---|---|
| Dinámica molecular | [LAMMPS/AIREBO](lammps/airebo/Readme.md) | Polietileno y potencial AIREBO | Benchmark ejecutable; requiere `CH.airebo` |
| Dinámica molecular | [LAMMPS/ClayFF](lammps/clayff/Readme.md) | C-S-H y C-S-H/PE | Inputs de prueba; la parametrización ClayFF completa está pendiente |
| Dinámica molecular | [LAMMPS/calcium](lammps/calcium/Readme.md) | Sistemas con calcio y C-S-H | Input y análisis RDF documentados |
| Dinámica molecular | [LAMMPS/water](lammps/water/Readme.md) | Agua con TIP4P | Input de agua molecular |
| Dinámica molecular | `lammps/methane` | Metano, enlaces y análisis geométrico | Inputs y scripts de análisis |
| Química cuántica | [NWChem/ASE](nwchem/README.md) | H₂, H₂O, hielo, calcio-silicato y tobermorita | Ejemplos Python; NWChem se instala aparte |
| Lattice Boltzmann | [OpenLB](LBM/openLB/Readme.md) | Fluidos, transporte y medios complejos | Guía; las fuentes se descargan externamente |
| Lattice Boltzmann | [Palabos](LBM/palabos/Readme.md) | Flujos, multifase y C-S-H/PET | Guía y casos metodológicos |
| Acoplamiento multiescala | [MD-DPD-DEM](coupling/MD-DPD-DEM/Readme.md) | C-S-H/PE desde átomos hasta sólido | Hoja de ruta de implementación |
| Simulación estadística | `molsim/Exercise_1` | Distribución de partículas en compartimentos | Ejemplo Octave/Monte Carlo |

## Relación entre las escalas

La propuesta científica del repositorio sigue esta lógica:

1. **NWChem/ASE** puede aportar geometrías, energías y referencias atomísticas para estructuras pequeñas o modelos de C-S-H.

2. **LAMMPS** permite estudiar polímeros, agua, calcio-silicatos, interfaces y propiedades estructurales o mecánicas a escala molecular.

3. **DPD** representa la evolución mesoscópica del gel y la formación de redes o percolación.

4. **OpenLB/Palabos** permiten estudiar transporte, multifase, mojabilidad y flujo en medios porosos; el caso C-S-H/PET usa Palabos como propuesta de modelo Shan-Chen.

5. **DEM** puede representar el estado sólido y la fractura usando propiedades constitutivas derivadas de las escalas anteriores.

La trazabilidad propuesta es MD → DPD → DEM, con LBM como herramienta complementaria para transporte y fenómenos multifásicos. Los parámetros no deben transferirse automáticamente sin calibración, conversión de unidades y validación independiente.

## Primeros pasos

### Preparar el entorno

Seguir [docs/README.md](docs/README.md) para instalar Conda/Mamba, dependencias Python, LAMMPS, NWChem, MPI, Distrobox y herramientas de Visual Studio Code.

### Ejecutar un primer caso LAMMPS

El benchmark AIREBO es el punto de entrada más corto:

```bash

cd lammps/airebo

test -s CH.airebo

lmp -in in.airebo -log log.airebo

python 01-temperature.py

python 02-more-variables.py

```

El sistema usa 32 640 átomos, un paso de 0.25 fs y 100 pasos. Sirve para comprobar la instalación, pero no constituye todavía una medición completa de propiedades mecánicas.

### Ejecutar un ejemplo NWChem

Después de activar el entorno Conda/Mamba:

```bash

cd nwchem/test-h2-ASE

python h2-ase.py

```

Este caso verifica que ASE puede comunicarse con NWChem y optimizar una molécula pequeña.

## Estado y límites del repositorio

- Los ejemplos de LAMMPS son puntos de partida y no todos tienen una parametrización validada para publicar resultados.

- El caso ClayFF actual utiliza un modelo Lennard-Jones genérico de prueba; no debe confundirse con ClayFF completo.

- Los clusters de NWChem son aproximaciones y requieren estudios de convergencia de base, funcional, carga, multiplicidad y periodicidad.

- Los casos C-S-H/PET de Palabos son esqueletos metodológicos que deben adaptarse a una versión concreta de la API y validarse antes de ejecutar estudios físicos.

- OpenLB y Palabos no tienen sus árboles fuente completos dentro de este repositorio.

- El plan MD-DPD-DEM referencia scripts futuros que aún deben implementarse.


## Convenciones de reproducibilidad

Para cada simulación se deben registrar:

- versión del software y del entorno Conda/Mamba;

- archivo de potencial, funcional, base y parámetros;

- semillas aleatorias y condiciones iniciales;

- unidades y conversiones entre escalas;

- tamaño de sistema, malla, paso temporal y número de pasos;

- archivos de entrada, logs, trayectorias y resultados procesados;

- criterios de convergencia y comparación con referencias.


El objetivo es que cada ejemplo pueda pasar de una prueba de instalación a un caso científicamente trazable sin ocultar qué partes son demostrativas y cuáles están listas para una simulación de producción.
## Guia de ejecucion

La guia para preparar los entornos y ejecutar los ejemplos de LAMMPS, NWChem, OpenLB y Palabos desde Linux, Distrobox, Conda/Mamba y Visual Studio Code esta en [docs/README.md](docs/README.md).

## Investigation

### Molecular Dynamics Method - LAMMPS

El método de dinámica molecular (MD) es una técnica computacional que simula el comportamiento de sistemas a nivel atómico o molecular resolviendo las ecuaciones de movimiento de Newton para cada partícula en el sistema. Permite investigar propiedades físicas, químicas y estructurales de materiales, biomoléculas y líquidos bajo diferentes condiciones. LAMMPS es uno de los paquetes más populares para MD, altamente escalable y versátil, soportando una gran variedad de potenciales de interacción, modelos atómicos y condiciones de frontera, ideal para simular materiales, polímeros y biomoléculas a gran escala.

### Lattice Boltzmann Method - OpenLB

El método de Lattice Boltzmann (LBM) es una herramienta numérica alternativa para la dinámica de fluidos computacional (CFD), basada en una representación mesoscópica de los fluidos. En vez de resolver directamente las ecuaciones de Navier-Stokes, LBM modela la evolución espacio-temporal de una función de distribución de partículas en una red discreta o "lattice". Así, captura fenómenos multiescala y complejos (como flujos multifásicos, interacción fluido-estructura y transporte con reacciones). OpenLB es un paquete de software abierto diseñado para ejecutar simulaciones LBM en arquitecturas paralelas, ideal para modelado de flujos complejos y geometrías irregulares.[6][10]

### Finite Volume Method - OpenFOAM

El método de volúmenes finitos (FVM) es una técnica numérica ampliamente utilizada en la simulación CFD, especialmente en OpenFOAM. Consiste en dividir el dominio espacial en volúmenes de control, sobre los cuales se integran las ecuaciones de conservación (masa, momento, energía). Los valores de campo (como velocidad y presión) se almacenan en los centros de cada celda, y los flujos a través de las caras se calculan explícita o implícitamente usando esquemas de interpolación y discretización flexibles. OpenFOAM implementa FVM sobre mallas no estructuradas y permite personalizar los esquemas de cada término de la ecuación, proporcionando poder y flexibilidad para aplicaciones científicas e ingenieriles avanzadas.[2][4][11]

***

Cada método y framework se especializa en dominios físicos y escalas diferentes: LAMMPS para fenómenos moleculares y materiales, OpenLB para microfluídica y flujos complejos, y OpenFOAM para simulaciones hidrodinámicas robustas en geometrías complejas a escala continua.

[1](https://www.wolfdynamics.com/training/OF_WS2020/traning_session2020.pdf)
[2](https://openfoamwiki.net/index.php/OpenFOAM_guide/Finite_volume_method_(OpenFOAM))
[3](https://gidropraktikum.narod.ru/Moukalled-et-al-FVM-OpenFOAM-Matlab.pdf)
[4](https://pearl.plymouth.ac.uk/cgi/viewcontent.cgi?article=2945&context=secam-research)
[5](https://repositum.tuwien.at/bitstream/20.500.12708/152237/1/Florian%20Tobias%20-%202023%20-%20Space-time%20finite%20volume%20method%20in%20OpenFOAM.pdf)
[6](https://www.openlb.net/lattice-boltzmann-methods/)
[7](https://www.youtube.com/watch?v=4v7xJulFCjM)
[8](https://www.youtube.com/watch?v=jfk4feD7rFQ)
[9](https://openfoamwiki.net/index.php/OpenFOAM_guide/Finite_volume_method)
[10](https://www.youtube.com/watch?v=oxaxoeDAiuo)
[11](https://www.cfd-online.com/Forums/openfoam/88437-openfoam-fem-fvm-fdm.html)
[12](https://www.cfd-online.com/Forums/openfoam/60907-lattice-boltzmann-approach.html)
[13](https://www.nas.nasa.gov/assets/nas/pdf/ams/2020/AMS_20201201_Krause.pdf)
