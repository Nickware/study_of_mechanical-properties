# Ejemplos de NWChem y ASE

Esta carpeta reúne ejemplos pequeños para probar NWChem desde Python mediante ASE y para generar estructuras que después pueden utilizarse como entradas de NWChem. Los casos sirven como pruebas de instalación, construcción de geometrías, optimización molecular y preparación de modelos de agua o silicato de calcio.

## 1. Relación entre ASE y NWChem

- **ASE** (`Atomic Simulation Environment`) crea y manipula estructuras atómicas, escribe archivos y puede controlar optimizaciones.
- **NWChem** es el programa de química cuántica que calcula energías, gradientes y optimizaciones DFT.
- `ase.calculators.nwchem.NWChem` conecta ambos programas cuando un script asigna un calculador NWChem a un objeto `Atoms`.
- Algunos scripts solo generan `.xyz` o `.nwi`; esos scripts no ejecutan NWChem por sí mismos.

Los resultados de estos ejemplos no deben interpretarse todavía como propiedades convergidas de C-S-H, hielo o agua. Son casos de prueba que requieren revisar funcional, base, carga, multiplicidad, convergencia y condiciones de frontera antes de usarlos científicamente.

## 2. Mapa de ejemplos

| Carpeta | Archivo principal | Qué hace | Ejecuta NWChem |
|---|---|---|---|
| `test-h2-ASE` | `h2-ase.py` | Optimiza una molécula H2 con ASE | Sí, mediante ASE |
| `test-h2o-ASE` | `h2o.py` | Optimiza una molécula de agua con DFT | Sí, mediante ASE |
| `test-h2o-nwchem` | `h2o-nwchem.py` | Genera un input y llama al ejecutable NWChem | Sí |
| `test-h2o-cube` | `h2o-cube.py` | Genera una caja periódica de ocho moléculas de agua | No |
| `test-h2o-cube` | `h2o-cube-v2.py` | Genera una supercelda simplificada de hielo Ih | No |
| `test-calcium-tiny` | `calcium-tiny.py` | Genera un cluster pequeño Ca-Si-O-H y su input | No, se ejecuta aparte |
| `test-calcium` | `calcium-silicate.py` | Genera un cluster C-S-H aproximado y un input de optimización | No, se ejecuta aparte |
| `test-calcium-ASE` | `calcium-ase.py` | Optimiza una estructura de tobermorita desde CIF | Sí, mediante ASE |

## 3. Requisitos

Se necesita Python 3.11 o una versión compatible con el entorno del proyecto, ASE y NumPy:

```bash
mamba activate mechanical-properties
mamba install -n mechanical-properties -c conda-forge python=3.11 ase numpy
```

Comprobar la instalación de Python y ASE:

```bash
python -c "import ase, numpy; print('ASE y NumPy OK')"
which python
python --version
```

Para los casos que realizan cálculos cuánticos también se necesita el ejecutable NWChem:

```bash
mamba install -n mechanical-properties -c conda-forge nwchem
which nwchem
nwchem --help
```

Si NWChem no está disponible para la plataforma mediante Conda/Mamba, instalarlo desde una distribución oficial y añadir su directorio `bin` al `PATH`. ASE debe poder encontrar el ejecutable mediante `which nwchem`.

## 4. Flujo general de ejecución

Desde la raíz del repositorio:

```bash
cd /tmp/study_of_mechanical-properties
mamba activate mechanical-properties
```

Cada script debe ejecutarse desde su propia carpeta, porque utiliza nombres de archivos relativos y escribe sus resultados en el directorio de trabajo actual.

Un flujo típico es:

1. Crear o leer una geometría.
2. Verificar número de átomos, fórmula, distancias y celda.
3. Asignar un calculador ASE cuando corresponda.
4. Ejecutar energía, gradientes u optimización.
5. Guardar la geometría final y el log.
6. Revisar convergencia antes de interpretar la energía.

## 5. Ejemplos moleculares

### 5.1 Molécula de hidrógeno: `test-h2-ASE`

El script `h2-ase.py` crea una molécula H2 con una distancia inicial de 0.7 Å, asigna el calculador NWChem con el funcional PBE y optimiza con BFGS hasta una fuerza máxima de 0.02 eV/Å.

```bash
cd /tmp/study_of_mechanical-properties/nwchem/test-h2-ASE
python h2-ase.py
ls -lh H2.xyz
```

Este es el caso más pequeño para verificar que ASE puede iniciar NWChem y obtener fuerzas. Conviene ejecutarlo antes de casos con agua o estructuras de C-S-H.

### 5.2 Molécula de agua con ASE: `test-h2o-ASE`

`h2o.py` construye H2O, la centra en una caja de vacío, utiliza PBE y la base 3-21g, y optimiza con BFGS usando el criterio `fmax=0.02`.

```bash
cd /tmp/study_of_mechanical-properties/nwchem/test-h2o-ASE
python h2o.py
ls -lh H2O.xyz
```

El archivo `H2O.xyz` contiene la geometría final escrita por ASE. Revisar el archivo de salida de NWChem generado por ASE si se necesita auditar iteraciones, energía y convergencia.

### 5.3 Molécula de agua con input directo: `test-h2o-nwchem`

`h2o-nwchem.py` genera `h2o.nwi` y ejecuta `nwchem h2o.nwi` mediante `subprocess`. Para conservar el output de forma explícita:

```bash
cd /tmp/study_of_mechanical-properties/nwchem/test-h2o-nwchem
python h2o-nwchem.py
```

El script puede producir el input y lanzar NWChem sin guardar el output en un archivo. Para ejecutar el input manualmente y conservar el resultado:

```bash
nwchem h2o.nwi > h2o.nwo
```

El input usa la base 3-21g, el funcional PBE96 y una multiplicidad singlete. Verificar que los parámetros sean adecuados antes de reutilizarlo.

## 6. Ejemplos de agua y hielo

### 6.1 Caja de agua: `h2o-cube.py`

Este script crea ocho moléculas de agua en una caja periódica de 10 Å y escribe `water_box.xyz`. No llama a NWChem; es un generador de geometría.

```bash
cd /tmp/study_of_mechanical-properties/nwchem/test-h2o-cube
python h2o-cube.py
ls -lh water_box.xyz
```

La densidad impresa es una estimación geométrica. Antes de usar esta estructura en un cálculo periódico, comprobar que las moléculas no se solapen, que la celda y las condiciones periódicas sean compatibles con el calculador y que el modelo electrónico incluya las condiciones necesarias.

### 6.2 Hielo Ih: `h2o-cube-v2.py`

Este script crea una celda cristalina simplificada con `ase.spacegroup.crystal`, la replica 2x2x2 y escribe `ice_ih.xyz`:

```bash
cd /tmp/study_of_mechanical-properties/nwchem/test-h2o-cube
python h2o-cube-v2.py
ls -lh ice_ih.xyz
```

La estructura es pedagógica. El propio script advierte que no representa necesariamente la red de hielo Ih completa, por lo que debe verificarse la estequiometría, las posiciones fraccionarias, los enlaces de hidrógeno y la densidad antes de realizar un estudio físico.

## 7. Ejemplos de calcio-silicato hidratado

### 7.1 Cluster pequeño: `test-calcium-tiny`

`calcium-tiny.py` construye una geometría Ca-Si-O-H, escribe `csh_small_inicial.xyz` y genera el input `csh_small.nwi` con `noautoz noautosym`. Este bloque es útil para verificar la generación de entradas con geometría explícita.

```bash
cd /tmp/study_of_mechanical-properties/nwchem/test-calcium-tiny
python calcium-tiny.py
nwchem csh_small.nwi > csh_small.nwo
```

Es un cluster aproximado, no una estructura cristalina de C-S-H. Revisar cargas, multiplicidad, base y geometría antes de extraer conclusiones químicas.

### 7.2 Cluster de silicato de calcio: `test-calcium`

`calcium-silicate.py` genera `csh_inicial.xyz` y `csh_opt.nwi`. El input usa DFT PBE96, base 3-21g, multiplicidad 2 y una optimización con `driver`.

```bash
cd /tmp/study_of_mechanical-properties/nwchem/test-calcium
python calcium-silicate.py
nwchem csh_opt.nwi > csh_opt.nwo
ls -lh csh_inicial.xyz csh_opt.nwi csh_opt.nwo
```

La geometría es una aproximación de cluster con hidrógenos añadidos; no sustituye una estructura periódica de tobermorita, jennita o C-S-H amorfo validada.

### 7.3 Tobermorita con ASE: `test-calcium-ASE`

`calcium-ase.py` lee `tobermorite_11A.cif`, centra la estructura, activa condiciones periódicas y utiliza ASE con NWChem para una optimización BFGS. Es el caso más cercano a una estructura de material del repositorio.

```bash
cd /tmp/study_of_mechanical-properties/nwchem/test-calcium-ASE
python calcium-ase.py
ls -lh csh_inicial.xyz csh_optimizado.xyz csh_optimizado.cif csh_opt.traj csh_opt.log
```

El script requiere que `tobermorite_11A.cif` esté en el mismo directorio. Guarda la geometría inicial, la trayectoria de optimización, un log y las estructuras finalizadas en XYZ y CIF.

La configuración actual utiliza PBE, base 3-21g, corrección de dispersión D3, carga cero y multiplicidad 1. Para un estudio serio se debe comprobar la convergencia con bases mayores, funcionales apropiados para sólidos, tamaño de celda y tratamiento periódico.

## 8. Inspección y postprocesamiento

### 8.1 Archivos que deben conservarse

Según el ejemplo, conservar:

- `.xyz`: geometría inicial o final.
- `.cif`: estructura cristalográfica con celda cuando corresponda.
- `.nwi`: entrada reproducible de NWChem.
- `.nwo`: salida completa de NWChem.
- `.traj`: trayectoria de optimización ASE.
- `.log`: historial del optimizador ASE.

### 8.2 Comprobaciones mínimas

Después de cada cálculo:

```bash
grep -iE "converg|error|failed|total energy|final energy" *.nwo *.log 2>/dev/null
```

Revisar también:

- que el cálculo terminó y no fue interrumpido;
- que las fuerzas finales cumplen el criterio elegido;
- que la energía no cambia por falta de convergencia SCF;
- que no hay átomos excesivamente próximos;
- que la carga y la multiplicidad son las previstas;
- que la geometría final conserva la química esperada.

Para inspección visual se pueden abrir XYZ o CIF con OVITO, VMD o Avogadro. Para trayectorias ASE, utilizar un lector compatible con `.traj` o convertir la trayectoria a XYZ desde Python.

### 8.3 Interpretación

Una energía de un cluster aislado no se puede comparar directamente con la energía de una celda periódica. Tampoco se deben comparar energías calculadas con distintas bases, funcionales, cargas o multiplicidades sin documentar esas diferencias. Para propiedades de materiales se necesitan además convergencia de celda, tamaño de supercelda y, cuando corresponda, cálculos periódicos consistentes.

## 9. Errores frecuentes

- **`No module named ase`:** activar el entorno correcto e instalar `ase`.
- **`nwchem: command not found`:** instalar NWChem o añadir su directorio `bin` al `PATH`.
- **No se encuentra el CIF:** ejecutar `calcium-ase.py` desde `test-calcium-ASE`.
- **El script genera archivos pero no calcula:** comprobar si el script solo es un generador y ejecutar manualmente el `.nwi`.
- **Fallo de convergencia SCF:** revisar funcional, base, multiplicidad, carga, geometría inicial y opciones de damping.
- **Fallo del optimizador BFGS:** revisar que el calculador entregue fuerzas y que la estructura inicial no tenga distancias imposibles.

## 10. Reproducibilidad

Registrar para cada cálculo:

- versión de NWChem y ASE;
- versión de Python y entorno Conda/Mamba;
- funcional, base, carga y multiplicidad;
- criterio de convergencia y número máximo de iteraciones;
- estructura inicial y archivos de entrada;
- energía, fuerzas y geometría final;
- memoria, número de procesos y tiempo de ejecución.

Estos ejemplos son puntos de partida para la cadena multiescala del repositorio. Antes de conectar sus resultados con LAMMPS, DPD o DEM, validar la estructura y las propiedades electrónicas con referencias independientes.
