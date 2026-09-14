# Guia de ejecucion del repositorio

Esta guia describe como preparar un entorno reproducible y ejecutar los ejemplos de este repositorio desde Linux y Visual Studio Code. La ruta recomendada es trabajar dentro de un entorno Conda/Mamba y usar Distrobox solo cuando se necesite aislar el sistema operativo o fijar una distribucion concreta.

## 1. Mapa de tecnologias

| Tecnologia | Ejemplos locales | Lenguaje | Proposito |
|---|---|---|---|
| LAMMPS + AIREBO | `lammps/airebo` | LAMMPS y Python | Dinamica molecular de polietileno |
| NWChem + ASE | `nwchem/` | Python y NWChem | Geometrias, energia y optimizacion DFT |
| OpenLB | `LBM/openLB/Readme.md` | C++ | Lattice Boltzmann para fluidos y transporte |
| Palabos | `LBM/palabos/Readme.md` | C++ | Lattice Boltzmann paralelo |

Los directorios `LBM/openLB` y `LBM/palabos` contienen principalmente documentacion de instalacion y comparacion. No incluyen en este repositorio un arbol fuente completo ni un caso C++ listo para compilar. Para ejecutarlos hay que descargar sus fuentes oficiales y usar sus ejemplos.

## 2. Requisitos del sistema

En Ubuntu o Debian se recomienda instalar las herramientas basicas:

```bash
sudo apt update
sudo apt install -y git curl wget build-essential cmake make g++ \
  openmpi-bin libopenmpi-dev
```

Para Visual Studio Code, instalar tambien las extensiones **Python**, **C/C++** y, si se trabaja dentro de contenedores, **Dev Containers** o la extension de Distrobox disponible en la instalacion local de VS Code.

## 3. Instalacion y configuracion de entornos

### Distrobox

Distrobox permite ejecutar una distribucion Linux aislada usando Podman o Docker, compartiendo el directorio de trabajo y el usuario con el sistema anfitrion. Es util cuando LAMMPS, NWChem o los compiladores de LBM requieren versiones distintas de las del sistema.

#### Instalacion base

En Ubuntu/Debian, instalar primero Podman y herramientas auxiliares:

```bash
sudo apt update
sudo apt install -y podman distrobox
```

Si el paquete `distrobox` no existe en la version de la distribucion, seguir la instalacion oficial de Distrobox. Comprobar:

```bash
distrobox --version
podman --version
```

#### Crear un entorno de trabajo

```bash
distrobox create --name mechanical-dev --image ubuntu:24.04
```

Entrar y preparar el contenedor:

```bash
distrobox enter mechanical-dev
sudo apt update
sudo apt install -y git curl wget build-essential cmake make g++ \
  python3 python3-venv python3-pip openmpi-bin libopenmpi-dev
```

Dentro del contenedor se puede instalar Conda o Mamba y crear el mismo entorno Python descrito abajo. El repositorio se comparte desde el directorio anfitrion; por eso los cambios permanecen fuera del contenedor.

Salir con `exit` y volver a entrar con:

```bash
distrobox enter mechanical-dev
```

No instalar simultaneamente varias implementaciones MPI sin revisar el `PATH`; escoger OpenMPI o MPICH por entorno para evitar mezclas de bibliotecas.

### Entorno Conda/Mamba

Con Conda o Mamba se crea y activa el entorno de trabajo del repositorio:

```bash
mamba create -n mechanical-properties -c conda-forge \
  python=3.11 pip numpy scipy pandas matplotlib ase
eval "$(mamba shell hook --shell bash)"
mamba activate mechanical-properties
python -m pip install lammps-logfile
```

Los comandos `mamba` pueden sustituirse por `conda` cuando se utilice Conda. No se recomienda mezclar paquetes de distintos canales sin necesidad; `conda-forge` debe ser el canal principal.

Comprobar el entorno:

```bash
python -c "import ase, numpy, pandas, matplotlib, lammps_logfile; print('Python OK')"
which python
python --version
```

### Visual Studio Code

Para Visual Studio Code, instalar las extensiones **Python**, **C/C++** y, si se trabaja dentro de contenedores, **Dev Containers** o la extension de Distrobox disponible en la instalacion local de VS Code.

#### Seleccionar el entorno Python

1. Abrir la raiz `/tmp/study_of_mechanical-properties` en VS Code.
2. Activar el entorno con `mamba activate mechanical-properties` en una terminal integrada.
3. Ejecutar `Python: Select Interpreter` y seleccionar el Python de `mechanical-properties`.
4. Abrir cualquier script Python y usar el boton **Run Python File**.

La terminal integrada debe iniciar en la raiz del repositorio o en el directorio del ejemplo. Los scripts que leen `log.airebo` deben ejecutarse desde `lammps/airebo`.

#### Usar VS Code dentro de Distrobox

Desde una terminal del contenedor, comprobar si el comando `code` esta disponible:

```bash
code --version
code /tmp/study_of_mechanical-properties
```

Si `code` no esta disponible dentro del contenedor, abrir VS Code en el anfitrion y usar una terminal integrada con:

```bash
distrobox enter mechanical-dev
cd /tmp/study_of_mechanical-properties
```

Tambien puede usarse la extension **Dev Containers** con una imagen equivalente, pero eso es un flujo diferente de Distrobox. No se deben mezclar automaticamente sus configuraciones: primero escoger si VS Code se conectara al contenedor o si solo se usara el contenedor desde la terminal integrada.

## 4. Instalacion de tecnologias y ejecucion de ejemplos

### LAMMPS + AIREBO

La forma mas directa para este repositorio es instalar el paquete de Conda/Mamba:

```bash
mamba install -n mechanical-properties -c conda-forge lammps
mamba activate mechanical-properties
which lmp
lmp -help
```

El nombre del ejecutable puede variar entre instalaciones. Comprobar tambien:

```bash
which lammps
which lmp_serial
```

El ejemplo AIREBO necesita el archivo de potencial `CH.airebo`. LAMMPS debe poder encontrarlo en su directorio de potenciales o mediante la ruta configurada en el sistema. Antes de ejecutar, localizarlo:

```bash
find "$CONDA_PREFIX" -name CH.airebo -print
```

Si no aparece, descargar o copiar el archivo desde una distribucion de LAMMPS compatible y colocar una copia en `lammps/airebo/`, o definir una ruta de potenciales apropiada. No se debe ejecutar el caso hasta comprobar que `pair_coeff * * CH.airebo C H` puede abrir ese archivo.

Desde la raiz del repositorio:

```bash
cd /tmp/study_of_mechanical-properties
mamba activate mechanical-properties
cd lammps/airebo
lmp -in in.airebo -log log.airebo
```

El input actual replica la estructura hasta 32 640 atomos y ejecuta un benchmark corto de 100 pasos. Usa `timestep 0.00025`, es decir, 0.25 fs. Esto es mas conservador para las vibraciones C-H, pero el benchmark sigue siendo demasiado corto para medir propiedades mecanicas.

Comprobar el resultado:

```bash
tail -30 log.airebo
grep -E "Dangerous builds|Loop time|Total wall time" log.airebo
```

Para analizar temperatura y energia:

```bash
python 01-temperature.py
python 2-more-variables.py
```

Los scripts leen `log.airebo`. Necesitan una sesion grafica para mostrar las figuras; en un servidor sin interfaz se puede configurar Matplotlib con un backend no interactivo y guardar las figuras en archivos.

### NWChem + ASE

Los scripts Python de `nwchem/` generan o usan estructuras y llaman al ejecutable externo `nwchem`. ASE es la interfaz Python; no reemplaza al programa NWChem.

Instalar ASE y localizar NWChem:

```bash
mamba activate mechanical-properties
mamba install -c conda-forge nwchem
python -m pip install ase
which nwchem
nwchem --help
```

Si NWChem no esta disponible en Conda para la plataforma utilizada, instalarlo desde los binarios o fuentes oficiales y añadir su directorio `bin` a `PATH`. Verificarlo con `which nwchem` antes de ejecutar cualquier script.

### Agua con ASE y NWChem

```bash
cd /tmp/study_of_mechanical-properties/nwchem/test-h2o-ASE
python h2o.py
```

El script optimiza una molecula de agua con DFT y guarda `H2O.xyz`.

### Agua con entrada NWChem directa

```bash
cd /tmp/study_of_mechanical-properties/nwchem/test-h2o-nwchem
python h2o-nwchem.py
```

Este caso genera `h2o.nwi` y llama a NWChem. La entrada puede ejecutarse directamente para conservar el output:

```bash
nwchem h2o.nwi > h2o.nwo
```

### Cluster de calcio-silicato

```bash
cd /tmp/study_of_mechanical-properties/nwchem/test-calcium
python calcium-silicate.py
nwchem csh_opt.nwi > csh_opt.nwo
```

El cluster es una geometria inicial aproximada, no un modelo cristalino validado de C-S-H. Para el ejemplo ASE de tobermorita se necesita el archivo `tobermorite_11A.cif` que ya se encuentra en `test-calcium-ASE`:

```bash
cd /tmp/study_of_mechanical-properties/nwchem/test-calcium-ASE
python calcium-ase.py
```

Los calculos DFT pueden consumir bastante memoria y tiempo. Primero ejecutar moleculas pequeñas y revisar convergencia, carga, multiplicidad, base y funcional antes de interpretar las energias como propiedades del material.

### OpenLB

OpenLB es una biblioteca C++ de Lattice Boltzmann. El README local propone instalar compilador, CMake, MPI y Git, clonar las fuentes y compilar un ejemplo.

```bash
cd /tmp
git clone https://gitlab.com/openlb/release.git olb-release
cd olb-release
```

La estructura exacta de ejemplos depende de la version descargada. Localizar un caso:

```bash
find . -maxdepth 5 -type f \( -name Makefile -o -name CMakeLists.txt \) | head -30
```

Seguir el README oficial de la version descargada. En un ejemplo basado en Makefile, el flujo suele ser:

```bash
cd ruta/al/ejemplo
make
./ejecutable
```

Para un caso MPI:

```bash
mpirun -np 4 ./ejecutable
```

Los resultados suelen escribirse en VTK y pueden inspeccionarse con ParaView. El contenido local de `LBM/openLB` es una guia, no un ejecutable del proyecto.

### Palabos

Palabos tambien es una biblioteca C++ de Lattice Boltzmann y normalmente se compila junto con sus aplicaciones de ejemplo:

```bash
cd /tmp
git clone https://gitlab.com/unigespc/palabos.git
cd palabos
```

La ruta y el nombre de los ejemplos pueden cambiar entre versiones. Para un caso con Makefile:

```bash
cd examples/showCases/laminarChannel/build
make
./laminarChannel
```

Para ejecutar en paralelo:

```bash
mpirun -np 4 ./laminarChannel
```

Revisar el `Makefile` del ejemplo para confirmar el nombre real del ejecutable. Los resultados VTK pueden abrirse con ParaView. El directorio local `LBM/palabos` documenta la tecnologia, pero no incluye las fuentes completas de Palabos.

### Ejecutar desde tareas de VS Code

Las tareas pueden invocar el entorno activado y dejar los comandos reproducibles. Como ejemplo, desde la terminal integrada:

```bash
mamba activate mechanical-properties
cd lammps/airebo
lmp -in in.airebo -log log.airebo
python 01-temperature.py
```

Para NWChem:

```bash
cd nwchem/test-h2o-ASE
python h2o.py
```

Para OpenLB o Palabos, abrir la carpeta de la fuente descargada y ejecutar `make` desde el directorio del ejemplo; para MPI, ejecutar `mpirun` con el numero de procesos que el equipo pueda soportar.

## 5. Notas sobre las simulaciones

### LAMMPS y propiedades mecanicas del polietileno

El benchmark no equivale a un estudio mecanico. Para obtener propiedades del polietileno se necesita, como minimo:

1. minimizacion de energia;
2. equilibrio NVT y, si procede, NPT;
3. produccion suficientemente larga y varias semillas;
4. ensayos separados de traccion en x, y y z, compresion y cizallamiento;
5. registro de `pxx`, `pyy`, `pzz`, deformacion, volumen y temperatura;
6. ajuste de la region elastica, calculo de modulo de Young y Poisson;
7. estudio de convergencia respecto al tamaño, velocidad de deformacion y tiempo;
8. incertidumbre estadistica y comparacion con datos experimentales.

El flujo detallado y los bloques LAMMPS de referencia estan en [lammps/airebo/Readme.md](../lammps/airebo/Readme.md).

### NWChem

El cluster de calcio-silicato es una geometria inicial aproximada, no un modelo cristalino validado de C-S-H. Los calculos DFT pueden consumir bastante memoria y tiempo. Primero ejecutar moleculas pequeñas y revisar convergencia, carga, multiplicidad, base y funcional antes de interpretar las energias como propiedades del material.

### OpenLB y Palabos

Los resultados suelen escribirse en VTK y pueden inspeccionarse con ParaView. Los directorios locales documentan las tecnologias, pero no incluyen las fuentes completas ni un caso C++ listo para compilar; los resultados dependen de la version externa descargada y del ejemplo seleccionado.

## 6. Lista de comprobacion

Antes de considerar un ejemplo reproducible, verificar:

- `python`, `lmp` o `nwchem` se encuentran con `which`;
- el entorno activo es el esperado con `which python`;
- `CH.airebo` esta disponible para el caso AIREBO;
- las fuentes de OpenLB o Palabos corresponden a la version documentada;
- MPI funciona con `mpirun --version`;
- los archivos de salida se generan en el directorio esperado;
- la simulacion termina sin errores, segmentaciones ni `Dangerous builds`;
- la version de cada software, semilla y parametros quedan registrados.
