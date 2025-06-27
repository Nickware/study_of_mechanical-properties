# Workflow Detallado: Multiscale MD-DPD para C-S-H/PE Composites

## FASE 1: PREPARACIÓN Y SETUP (Días 1-3)

### Día 1: Literatura y Modelos
**Tareas:**
- [ ] Revisión papers recientes MD de C-S-H (últimos 2 años)
- [ ] Revisión métodos coarse-graining MD→DPD
- [ ] Descarga/preparación force fields (CLAYFF, OPLS-AA)
- [ ] Setup entorno computacional (LAMMPS, scripts)

**Deliverable:** Base de datos literatura + entorno listo

### Día 2-3: Construcción Sistemas Atómicos
**Tareas:**
- [ ] Construir celda C-S-H (Jennite/Tobermorite, ~5000 átomos)
- [ ] Construir cadenas PE (C₂₀₀-C₅₀₀, diferentes pesos moleculares)
- [ ] Crear sistema interfaz C-S-H/PE (geometría laminar)
- [ ] Verificar topologías y parámetros

**Scripts necesarios:**
```bash
# Construcción C-S-H
python build_csh.py --model jennite --size 40x40x20
# Construcción PE
python build_pe.py --chains 50 --length 200
# Interfaz
python build_interface.py --csh csh.data --pe pe.data
```

## FASE 2: SIMULACIONES MD (Días 4-10)

### Día 4-5: Equilibración Sistemas Puros
**C-S-H puro:**
```bash
# Minimización energía
lmp_serial -in minimize_csh.in
# Equilibración NPT
lmp_serial -in equilibrate_csh.in  # 300K, 1atm, 10ns
```

**PE puro:**
```bash
# Minimización
lmp_serial -in minimize_pe.in
# Equilibración NPT
lmp_serial -in equilibrate_pe.in  # 300K, 1atm, 10ns
```

**Análisis convergencia:**
- [ ] Densidad vs tiempo
- [ ] Energía potencial
- [ ] Presión del sistema

### Día 6-8: Interfaz C-S-H/PE
```bash
# Equilibración interfaz (más larga)
lmp_serial -in equilibrate_interface.in  # 300K, 1atm, 20ns
```

**Métricas específicas:**
- [ ] Estabilidad interfaz
- [ ] Perfiles densidad
- [ ] Orientación cadenas PE

### Día 9-10: Propiedades Mecánicas
**Para cada sistema:**
```bash
# Deformación uniaxial (tensión)
lmp_serial -in tension_x.in  # εₓₓ = 0-5%
lmp_serial -in tension_y.in  # εᵧᵧ = 0-5%
lmp_serial -in tension_z.in  # ε_zz = 0-5%

# Deformación hidrostática
lmp_serial -in bulk_compression.in  # P = 0-1 GPa
```

**Extracción datos:**
- [ ] Módulos elásticos (E, G, K)
- [ ] Coeficientes Poisson
- [ ] Energía adhesión interfaz

## FASE 3: ANÁLISIS Y COARSE-GRAINING (Días 11-14)

### Día 11: Análisis Estructural Detallado
**Scripts Python:**
```python
# Análisis propiedades MD
python analyze_mechanical.py --input stress_strain.dat
python analyze_interface.py --input trajectory.dcd
python calculate_adhesion.py --input interface_energy.dat
```

**Outputs esperados:**
- E_CSH ≈ 25 ± 5 GPa
- E_PE ≈ 1.0 ± 0.2 GPa
- Adhesión ≈ 0.1-0.5 J/m²

### Día 12-13: Mapeo Coarse-Graining
**Estrategia mapping:**
```
C-S-H: 1 bead = ~10-20 átomos (Ca-Si-O cluster)
PE: 1 bead = ~5-10 átomos CH₂ (backbone)
```

**Cálculo parámetros DPD:**
```python
# Método Flory-Huggins
python calculate_chi_parameter.py --md_data interface_analysis.dat
# Parámetros conservativos
python derive_dpd_params.py --chi_data chi_values.dat --density rho_md.dat
```

### Día 14: Parámetros DPD Finales
**Tabla parámetros:**
```
a_AA (C-S-H/C-S-H): ~15-25 kBT/rc
a_BB (PE/PE): ~25-35 kBT/rc  
a_AB (C-S-H/PE): ~30-50 kBT/rc
γ: 4.5 kBT·τ/rc²
σ: 3.0 kBT/rc
```

## FASE 4: SIMULACIONES DPD (Días 15-21)

### Día 15-16: Setup y Validación DPD
```bash
# Sistema DPD puro A (C-S-H)
lmp_serial -in dpd_pure_a.in  # ρ=3, 50000 partículas
# Sistema DPD puro B (PE)  
lmp_serial -in dpd_pure_b.in  # ρ=3, 50000 partículas
```

**Validación termodinámica:**
- [ ] Densidad: ρ_DPD ≈ ρ_MD
- [ ] Compresibilidad: κ_DPD ≈ κ_MD
- [ ] Ecuación estado

### Día 17-19: Composites DPD
```bash
# Diferentes fracciones volumétricas
lmp_serial -in dpd_composite_10.in  # 10% PE
lmp_serial -in dpd_composite_30.in  # 30% PE
lmp_serial -in dpd_composite_50.in  # 50% PE
```

**Análisis morfología:**
- [ ] Distribución fases
- [ ] Conectividad/percolación
- [ ] Interfaces DPD vs MD

### Día 20-21: Propiedades Mecánicas DPD
```bash
# Deformación sistemas composite
lmp_serial -in dpd_mechanical.in
```

**Comparación MD-DPD:**
- [ ] Módulos efectivos
- [ ] Respuesta no-lineal
- [ ] Dependencia composición

## FASE 5: ANÁLISIS Y REDACCIÓN (Días 22-28)

### Día 22-23: Análisis Integral
**Comparación sistemática:**
```python
python compare_scales.py --md_data md_results/ --dpd_data dpd_results/
python plot_validation.py --comparison comparison_data.json
python sensitivity_analysis.py --params dpd_params.dat
```

**Métricas clave:**
- Error relativo E_DPD vs E_MD < 20%
- Morfología cualitativa similar
- Escalabilidad demostrada

### Día 24-25: Figuras y Tablas
**Figuras principales:**
1. **Setup sistemas**: MD snapshots, DPD mapping
2. **Propiedades MD**: Stress-strain, interfaz
3. **Validación DPD**: Density, RDF comparison  
4. **Composite properties**: E vs composition
5. **Multiscale comparison**: MD vs DPD morfología

**Tablas:**
1. Parámetros MD (force fields, condiciones)
2. Propiedades mecánicas extraídas
3. Parámetros DPD derivados
4. Comparación cuantitativa MD-DPD

### Día 26-28: Redacción Manuscrito

#### Día 26: Secciones principales
- [ ] Abstract (150 palabras)
- [ ] Introduction (400 palabras)
- [ ] Methodology (1200 palabras)

#### Día 27: Resultados y discusión
- [ ] Results MD (800 palabras)
- [ ] Scale transfer (400 palabras)  
- [ ] DPD validation (600 palabras)

#### Día 28: Finalización
- [ ] Conclusions (200 palabras)
- [ ] References (40-50 refs)
- [ ] Formatting JCP
- [ ] Revisión final

## RECURSOS COMPUTACIONALES

### Hardware necesario:
- **MD simulations**: 32-64 cores, 64 GB RAM
- **DPD simulations**: 16-32 cores, 32 GB RAM
- **Storage**: ~100 GB datos + backups

### Software stack:
```bash
# Simulación
LAMMPS (latest stable)
VMD (visualization)
OVITO (analysis)

# Análisis  
Python 3.8+ (NumPy, SciPy, MDAnalysis)
MATLAB (optional, mechanical analysis)
Origin/Gnuplot (plotting)
```

### Archivos clave generados:
```
md_results/
├── csh_properties.dat
├── pe_properties.dat  
├── interface_analysis.dat
└── mechanical_data/
    ├── stress_strain_csh.dat
    ├── stress_strain_pe.dat
    └── adhesion_energy.dat

dpd_results/
├── validation_thermodynamic.dat
├── composite_morphology/
└── mechanical_comparison.dat

manuscripts/
├── figures/
├── tables/
└── manuscript_v1.tex
```

## CONTINGENCIAS Y ALTERNATIVAS

### Si MD toma más tiempo:
- Reducir tamaño sistemas (trade-off estadística)
- Paralelizar simulaciones diferentes T, P
- Usar datos literatura para validación cruzada

### Si DPD no converge:
- Ajustar parámetros γ, σ iterativamente  
- Probar diferentes mappings (coarser/finer)
- Método alternativo (United Atom → DPD)

### Si análisis es complejo:
- Foco en propiedades principales (E, adhesión)
- Morfología cualitativa vs cuantitativa
- Colaboración análisis especializado

## CRITERIOS DE ÉXITO

### Técnicos:
- ✅ Simulaciones MD convergen (< 5% drift)
- ✅ Parámetros DPD reproducen ρ, κ (< 10% error)
- ✅ Composite DPD muestra física razonable

### Científicos:
- ✅ Metodología transferible otros sistemas
- ✅ Validación cuantitativa MD-DPD
- ✅ Insights nuevos C-S-H/PE interfaces

### Publicación:
- ✅ 8 páginas, formato JCP
- ✅ 5-6 figuras alta calidad
- ✅ Código/datos disponibles (GitHub)
- ✅ Submission antes deadline

¡Este workflow indica una roadmap clara día a día para completar exitosamente el proyecto!