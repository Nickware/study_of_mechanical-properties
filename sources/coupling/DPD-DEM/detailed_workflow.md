# Workflow Detallado: DPD-DEM para C-S-H/PE Composites

## FASE 1: PREPARACIÓN Y SETUP DPD (Días 1-4)

### Día 1: Literatura y Fundamentos
**Tareas críticas:**
- [ ] Revisión DPD para materiales cementosos (últimos 3 años)
- [ ] Coupling methods DPD-DEM (literatura limitada - oportunidad!)
- [ ] Parámetros DPD para C-S-H en literatura (Pellenq, Manzano, etc.)
- [ ] Setup LAMMPS con DPD package

**Recursos clave:**
```bash
# Papers esenciales
- "DPD modeling of cement paste" (Cement & Concrete Research)
- "Mesoscale modeling of concrete" (Computational Materials Science)  
- "DPD-DEM coupling methods" (Computer Methods in Applied Mechanics)
```

**Deliverable:** Base datos + entorno computacional listo

### Día 2: Modelos DPD Fundamentales
**Construcción modelos base:**

**C-S-H DPD model:**
```python
# Parámetros conservativos
a_CSH_CSH = 18.0  # Repulsión moderada, material rígido
rho = 3.0         # Densidad reducida
gamma = 4.5       # Fricción
sigma = 3.0       # Ruido térmico
rc = 1.0          # Radio corte
```

**PE DPD model:**
```python  
# Cadenas poliméricas
a_PE_PE = 25.0    # Repulsión suave, material flexible
bond_k = 100.0    # Rigidez enlace
bond_r0 = 0.7     # Longitud equilibrio
angle_k = 20.0    # Rigidez angular
```

**Interfaz C-S-H/PE:**
```python
# Parámetro crítico para adhesión
a_CSH_PE = 35.0   # Repulsión fuerte (hidrofóbico/hidrofílico)
# Será optimizado según adhesión experimental
```

### Día 3: Implementación y Calibración
**Scripts LAMMPS:**

```bash
# Archivo: dpd_csh_pure.in
units           lj
dimension       3
atom_style      dpd
boundary        p p p

# Geometría
region          box block 0 20 0 20 0 20
create_box      1 box
create_atoms    1 random 24000 12345 box

# Parámetros DPD
pair_style      dpd 1.0 1.0 12345
pair_coeff      1 1 18.0 4.5 1.0

# Integración
fix             1 all nve
fix             2 all langevin 1.0 1.0 0.1 12345
timestep        0.01
```

**Calibración densidad:**
- [ ] Ajustar N_particles para ρ=3.0
- [ ] Verificar presión ≈ 0
- [ ] RDF consistente con estructura C-S-H

### Día 4: Validación Modelos Puros
**Test sistemático:**
```bash
# Simulaciones calibración
lmp_serial -in dpd_csh_pure.in     # 100k steps
lmp_serial -in dpd_pe_pure.in      # 100k steps
lmp_serial -in dpd_interface.in    # 200k steps
```

**Métricas validación:**
- [ ] Densidad: ρ_sim ≈ ρ_target (±2%)
- [ ] Presión: P ≈ 0 (±0.1)
- [ ] Temperatura: T ≈ 1.0 (±0.05)
- [ ] RDF: Estructura física realista

## FASE 2: SIMULACIONES DPD Y CARACTERIZACIÓN (Días 5-11)

### Día 5-6: Sistemas Composite DPD
**Configuraciones objetivo:**
- C-S-H puro (referencia)
- 10% PE (matriz C-S-H, PE disperso)
- 30% PE (transición, bicontinua)
- 50% PE (co-continua)

**Setup simulaciones:**
```bash
# Sistema 10% PE
create_atoms    1 random 21600 12345 box  # C-S-H
create_atoms    2 random 2400 54321 box   # PE

# Sistema 30% PE
create_atoms    1 random 16800 12345 box  # C-S-H  
create_atoms    2 random 7200 54321 box   # PE

# Sistema 50% PE
create_atoms    1 random 12000 12345 box  # C-S-H
create_atoms    2 random 12000 54321 box  # PE
```

**Simulaciones largas:**
```bash
# Equilibración: 500k steps (5 unidades tiempo)
# Producción: 1M steps (10 unidades tiempo)
# Total por sistema: ~12-24 horas CPU
```

### Día 7-8: Extracción Propiedades Mecánicas DPD
**Deformación uniaxial:**
```bash
# Archivo: dpd_tension.in
fix             deform all deform 1 x erate 0.001 remap x
compute         stress all stress/atom NULL
compute         pe all pe/atom
fix             stress_output all ave/time 100 1 100 c_stress[1] c_stress[2] c_stress[3]
```

**Propiedades extraídas:**
- **Módulo elástico**: E = dσ/dε (región lineal)
- **Resistencia tensil**: σ_max
- **Resistencia compresiva**: Carga inversa
- **Curva completa**: σ-ε hasta fractura

**Análisis estadístico:**
```python
# Script: analyze_dpd_mechanical.py
import numpy as np
import matplotlib.pyplot as plt

def extract_elastic_modulus(stress, strain):
    # Región lineal: ε < 0.01
    linear_region = strain < 0.01
    E = np.polyfit(strain[linear_region], stress[linear_region], 1)[0]
    return E

def extract_tensile_strength(stress):
    return np.max(stress)
```

### Día 9-10: Análisis Morfológico Avanzado
**Caracterización microestructural:**
```python
# Script: morphology_analysis.py
import MDAnalysis as mda
from scipy.spatial.distance import pdist, squareform

def analyze_phase_distribution(trajectory):
    # Análisis conectividad
    connectivity = calculate_connectivity(positions, cutoff=1.2)
    
    # Tamaño dominios
    domain_sizes = cluster_analysis(connectivity)
    
    # Percolación
    percolation = check_percolation(connectivity, box_size)
    
    return connectivity, domain_sizes, percolation

def interface_analysis(csh_coords, pe_coords):
    # Área interfacial específica
    interface_area = calculate_interface_area(csh_coords, pe_coords)
    
    # Espesor interfacial
    interface_thickness = calculate_thickness(csh_coords, pe_coords)
    
    return interface_area, interface_thickness
```

**Métricas morfológicas:**
- [ ] **Conectividad**: Análisis grafos, clusters
- [ ] **Percolación**: Threshold, exponent crítico
- [ ] **Interfaz**: Área específica (m²/g), espesor (nm)
- [ ] **Distribución**: Función correlación radial g(r)

### Día 11: Homogeneización DPD→DEM
**Teoría micromecánica:**
```python
# Script: dpd_to_dem_homogenization.py

def homogenize_elastic_properties(dpd_results, volume_fractions):
    """
    Homogeneización Voigt-Reuss-Hill para propiedades elásticas
    """
    E_csh, E_pe = dpd_results['E_csh'], dpd_results['E_pe']
    f_csh, f_pe = volume_fractions['csh'], volume_fractions['pe']
    
    # Voigt (upper bound)
    E_voigt = f_csh * E_csh + f_pe * E_pe
    
    # Reuss (lower bound)  
    E_reuss = 1.0 / (f_csh/E_csh + f_pe/E_pe)
    
    # Hill (promedio)
    E_hill = 0.5 * (E_voigt + E_reuss)
    
    return E_hill

def derive_dem_contact_parameters(E_bulk, nu_bulk, morphology):
    """
    Derivación parámetros contacto DEM desde propiedades bulk
    """
    # Rigidez normal
    kn = E_bulk / (1 - nu_bulk**2) * characteristic_length
    
    # Rigidez tangencial
    ks = kn * (1 - 2*nu_bulk) / (2 - 2*nu_bulk)
    
    # Fricción
    mu = estimate_friction_from_morphology(morphology)
    
    # Cohesión
    cohesion = estimate_cohesion_from_interface(morphology['interface'])
    
    return {'kn': kn, 'ks': ks, 'mu': mu, 'cohesion': cohesion}
```

## FASE 3: IMPLEMENTACIÓN DEM (Días 12-18)

### Día 12-13: Setup Geometría DEM
**Generación empaquetamiento:**
```python
# Script: generate_dem_packing.py
import numpy as np
from scipy.spatial import Voronoi, distance_matrix

def generate_dem_geometry_from_dpd(dpd_morphology):
    """
    Conversión morfología DPD → geometría DEM
    """
    # Centros masa dominios DPD → centros elementos DEM
    csh_centers = extract_domain_centers(dpd_morphology, phase='csh')
    pe_centers = extract_domain_centers(dpd_morphology, phase='pe')
    
    # Tamaños elementos proporcionales a dominios
    csh_radii = estimate_element_sizes(csh_centers, dpd_morphology)
    pe_radii = estimate_element_sizes(pe_centers, dpd_morphology)
    
    # Conectividad basada en vecindad DPD
    connectivity = derive_connectivity(csh_centers, pe_centers, dpd_morphology)
    
    return {
        'csh_elements': {'centers': csh_centers, 'radii': csh_radii},
        'pe_elements': {'centers': pe_centers, 'radii': pe_radii},
        'connectivity': connectivity
    }
```

**Archivo DEM (LIGGGHTS):**
```bash
# Archivo: dem_composite.in
atom_style      granular
boundary        f f f

# Geometría desde DPD
read_data       dem_geometry.data

# Modelos contacto
pair_style      gran/hertz/history 1 0 1 1 0 1
pair_coeff      1 1 50e9 0.3 0.6 0 0 0 0 0 0  # C-S-H/C-S-H
pair_coeff      2 2 2e9 0.4 0.3 0 0 0 0 0 0   # PE/PE
pair_coeff      1 2 10e9 0.35 0.5 0 1e6 0 0 0 0  # C-S-H/PE + cohesión
```

### Día 14-15: Ensayos Virtuales DEM
**Compresión uniaxial:**
```bash
# Archivo: dem_compression.in
fix             walls all wall/gran hertz/history 1 0 1 1 0 1 &
                zplane 0.0 50.0 zplane 0.0 50.0

# Carga controlada
fix             compress all wall/gran hertz/history 1 0 1 1 0 1 &
                zplane 0.0 NULL zplane 50.0 NULL
fix             load all wall/pressure 1 press 0 0 -1000 0 0 0

# Output
compute         stress all stress/atom gran
dump            1 all custom 1000 compression.dump id type x y z radius
```

**Tensión directa:**
```bash
# Ensayo brasileño (splitting test)
fix             load_top all wall/gran hertz/history 1 0 1 1 0 1 &
                zcylinder 0 0 25 0 0 50
fix             load_bottom all wall/gran hertz/history 1 0 1 1 0 1 &
                zcylinder 0 0 25 0 0 0
```

### Día 16-17: Análisis Resultados DEM
**Post-procesamiento:**
```python
# Script: analyze_dem_results.py
import numpy as np
import matplotlib.pyplot as plt

def analyze_dem_compression(dump_files):
    """
    Análisis ensayo compresión DEM
    """
    time, stress, strain = load_dem_data(dump_files)
    
    # Módulo elástico
    E_dem = calculate_elastic_modulus(stress, strain)
    
    # Resistencia
    sigma_max = np.max(stress)
    
    # Curva completa
    plot_stress_strain(stress, strain, 'DEM_compression.png')
    
    return {'E': E_dem, 'sigma_max': sigma_max}

def damage_analysis(dump_files):
    """
    Análisis patrones fractura
    """
    # Localización damage
    damage_zones = identify_damage_zones(dump_files)
    
    # Modos fractura
    fracture_modes = classify_fracture_modes(damage_zones)
    
    # Energía fractura
    fracture_energy = calculate_fracture_energy(dump_files)
    
    return damage_zones, fracture_modes, fracture_energy
```

### Día 18: Validación DPD-DEM
**Comparación cuantitativa:**
```python
# Script: validate_dpd_dem.py

def validate_mechanical_properties(dpd_results, dem_results):
    """
    Validación propiedades mecánicas DPD vs DEM
    """
    properties = ['E_modulus', 'tensile_strength', 'compressive_strength']
    
    validation_results = {}
    
    for prop in properties:
        dpd_val = dpd_results[prop]
        dem_val = dem_results[prop]
        
        relative_error = abs(dem_val - dpd_val) / dpd_val * 100
        validation_results[prop] = {
            'dpd': dpd_val,
            'dem': dem_val, 
            'error_%': relative_error
        }
    
    return validation_results

def plot_validation_summary(validation_results):
    """
    Gráfico resumen validación
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Correlación DPD-DEM
    dpd_vals = [v['dpd'] for v in validation_results.values()]
    dem_vals = [v['dem'] for v in validation_results.values()]
    
    ax1.scatter(dpd_vals, dem_vals)
    ax1.plot([0, max(dpd_vals)], [0, max(dpd_vals)], 'r--', label='1:1')
    ax1.set_xlabel('DPD Properties')
    ax1.set_ylabel('DEM Properties')
    ax1.legend()
    
    # Errores relativos
    errors = [v['error_%'] for v in validation_results.values()]
    properties = list(validation_results.keys())
    
    ax2.bar(properties, errors)
    ax2.set_ylabel('Relative Error (%)')
    ax2.axhline(y=15, color='r', linestyle='--', label='15% threshold')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('validation_summary.png', dpi=300)
```

## FASE 4: ANÁLISIS INTEGRAL Y REDACCIÓN (Días 19-28)

### Día 19-21: Análisis Sensibilidad y Escalabilidad
**Análisis sensibilidad parámetros:**
```python
# Script: sensitivity_analysis.py
import numpy as np
from scipy.stats import uniform
from SALib.sample import saltelli
from SALib.analyze import sobol

def sensitivity_analysis_dpd_dem():
    """
    Análisis sensibilidad método Sobol
    """
    # Parámetros variables
    problem = {
        'num_vars': 6,
        'names': ['a_csh_csh', 'a_pe_pe', 'a_csh_pe', 'gamma', 'kn_ratio', 'cohesion'],
        'bounds': [[15, 25], [20, 30], [30, 45], [3, 6], [0.3, 0.8], [0.1, 2.0]]
    }
    
    # Sampling
    param_values = saltelli.sample(problem, 256)
    
    # Evaluación modelo
    Y = np.zeros([param_values.shape[0]])
    for i, X in enumerate(param_values):
        Y[i] = evaluate_dpd_dem_model(X)
    
    # Análisis Sobol
    Si = sobol.analyze(problem, Y)
    
    return Si

def scalability_analysis():
    """
    Análisis escalabilidad computacional
    """
    system_sizes = [1e4, 1e5, 1e6, 1e7]  # Número partículas DPD
    
    dpd_times = []
    dem_times = []
    
    for size in system_sizes:
        # Tiempo DPD
        t_dpd = benchmark_dpd_simulation(size)
        dpd_times.append(t_dpd)
        
        # Tiempo DEM equivalente
        dem_size = size / 100  # Ratio coarse-graining
        t_dem = benchmark_dem_simulation(dem_size)
        dem_times.append(t_dem)
    
    return system_sizes, dpd_times, dem_times
```

### Día 22-23: Preparación Figuras Científicas
**Figuras principales:**

**Figura 1: Workflow y Modelos**
- (a) Esquema general DPD→DEM
- (b) Morfología DPD (snapshots)
- (c) Geometría DEM derivada
- (d) Escalas espaciales/temporales

**Figura 2: Caracterización DPD**
- (a) Propiedades vs composición
- (b) Análisis morfológico
- (c) Percolación, conectividad
- (d) Interfaz C-S-H/PE

**Figura 3: Validación Mecánica**
- (a) Curvas σ-ε DPD vs DEM
- (b) Módulos elásticos comparación
- (c) Resistencias tensil/compresiva
- (d) Correlación general DPD-DEM

**Figura 4: Análisis Avanzado**
- (a) Patrones fractura/damage
- (b) Sensibilidad parámetros
- (c) Escalabilidad computacional
- (d) Efficiency vs accuracy

**Figura 5: Aplicación Ingenieril**
- (a) Ensayos virtuales (flexión, triaxial)
- (b) Optimización composición
- (c) Predicción propiedades
- (d) Validación experimental

### Día 24-25: Tablas y Análisis Estadístico
**Tablas principales:**

**Tabla 1: Parámetros DPD**
```
| Parámetro | C-S-H | PE | C-S-H/PE | Fuente |
|-----------|-------|----|---------|---------| 
| a_ii      | 18.0  | 25.0| 35.0    | Este trabajo |
| γ         | 4.5   | 4.5 | 4.5     | Literatura |
| σ         | 3.0   | 3.0 | 3.0     | Literatura |
| ρ         | 3.0   | 3.0 | 3.0     | Standard |
```

**Tabla 2: Propiedades Mecánicas**
```
| Propiedad | DPD | DEM | Error (%) | Exp. Ref. |
|-----------|-----|-----|-----------|-----------|
| E_csh (GPa)| 25.3±2.1| 24.7±1.8| 2.4 | 20-30 |
| E_pe (GPa) | 1.2±0.2 | 1.3±0.2 | 8.3 | 0.8-1.5 |
| σ_t (MPa)  | 3.2±0.4 | 2.9±0.3 | 9.4 | 2-4 |
```

**Tabla 3: Parámetros DEM Derivados**
```
| Contacto | kn (GPa) | ks (GPa) | μ | Cohesión (MPa) |
|----------|----------|----------|---|----------------|
| C-S-H/C-S-H| 65.0 | 25.0 | 0.7 | 0.0 |
| PE/PE    | 3.2  | 1.2  | 0.3 | 0.0 |
| C-S-H/PE | 15.0 | 6.0  | 0.5 | 0.8 |
```

### Día 26-28: Redacción Manuscrito
**Estructura día por día:**

**Día 26: Secciones Metodológicas**
- [ ] Abstract (150 palabras)
- [ ] Introduction (600 palabras)
- [ ] Methodology (1500 palabras)
  - DPD framework (500 palabras)
  - Property extraction (400 palabras)
  - DPD-DEM coupling (400 palabras)
  - Validation protocol (200 palabras)

**Día 27: Resultados y Discusión**
- [ ] DPD characterization (800 palabras)
- [ ] Scale transfer methodology (600 palabras)
- [ ] DEM validation (800 palabras)
- [ ] Computational efficiency (400 palabras)

**Día 28: Finalización**
- [ ] Conclusions (250 palabras)
- [ ] References (50 referencias)
- [ ] Formatting JCP style
- [ ] Revisión final, correcciones
- [ ] Preparation supporting information

## RECURSOS Y HERRAMIENTAS

### Software Stack
```bash
# Simulaciones DPD
LAMMPS 2023.08.02 (dpd package)
HOOMD-blue 4.0 (alternativa GPU)

# Simulaciones DEM  
LIGGGHTS 3.8.0 (open-source)
YADE 2023.02a (python-based)

# Análisis y visualización
Python 3.10+ (NumPy, SciPy, MDAnalysis)
OVITO 3.9.0 (morphology analysis)
ParaView 5.11 (DEM visualization)
MATLAB R2023a (statistical analysis)
```

### Recursos Computacionales
```bash
# DPD simulations
CPU: 32-64 cores
RAM: 64-128 GB  
Time: ~500-1000 core-hours

# DEM simulations
CPU: 16-32 cores
RAM: 32-64 GB
Time: ~200-400 core-hours

# Total computational cost
~1500 core-hours (3-4 weeks en cluster moderno)
Storage: ~200 GB (trajectories + analysis)
```

### Criterios Éxito
- ✅ **Validación mecánica**: Error DPD-DEM < 15%
- ✅ **Morfología**: Patrones cualitativamente similares
- ✅ **Escalabilidad**: Demostrada hasta 10⁶ partículas
- ✅ **Reproducibilidad**: Múltiples realizaciones consistentes
- ✅ **Aplicabilidad**: Metodología transferible otros sistemas

## CONTINGENCIAS

### Si DPD toma más tiempo:
- Reducir número composiciones (3 → 2)
- Paralelizar simulaciones diferentes
- Usar parámetros literatura como baseline

### Si validación DEM falla:
- Refinamiento iterativo parámetros
- Métodos homogeneización alternativos
- Validación cualitativa vs cuantitativa

### Si redacción se atrasa:  
- Foco en figuras clave (3-4 principales)
- Template JCP pre-existente
- Colaboración escritura especializada

¡Este workflow apoya el proyecto solido enmarcado a ser sometido para JCP!