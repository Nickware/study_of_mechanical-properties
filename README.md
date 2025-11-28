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
