# Workflow Detallado: MD-DPD-DEM Multiescala para C-S-H/PE Composites

## FASE 1: ESCALA ATOMÍSTICA — MOLECULAR DYNAMICS (Días 1-6)

### Día 1: Literatura y Setup
**Tareas:**
- [ ] Revisión force fields C-S-H (CLAYFF, Cygan et al.) y PE (OPLS-AA, PCFF)
- [ ] Revisión literatura reciente MD C-S-H/polímero (interfaz, adhesión)
- [ ] Setup LAMMPS (compilación con paquetes MOLECULE, MC, DPD-BASIC)
- [ ] Preparación estructuras cristalográficas base (tobermorita/jennita)

**Deliverable:** Entorno LAMMPS listo + estructuras base descargadas/construidas

### Día 2-3: Construcción de Sistemas Atómicos
**C-S-H:**
```bash
# Celda tobermorita 14Å, supercelda ~5000 átomos
python build_csh_structure.py --model tobermorite --replicate 4x4x2
```
**PE:**
```bash
# Cadenas C₂₀₀-C₅₀₀, empaquetamiento amorfo
python build_pe_chains.py --n_chains 40 --chain_length 300
```
**Interfaz:**
```bash
python build_interface_md.py --csh csh.data --pe pe.data --separation 3.0
```

### Día 4: Equilibración
```bash
lmp -in minimize_csh.in
lmp -in minimize_pe.in
lmp -in equilibrate_interface.in   # NPT, 300K, 1atm, 5-10 ns
```
**Monitoreo:** densidad, energía potencial, presión — convergencia antes de continuar.

### Día 5: Extracción de Propiedades Mecánicas
```bash
lmp -in tension_x.in    # deformación uniaxial ε=0-5%
lmp -in tension_y.in
lmp -in tension_z.in
lmp -in stress_fluctuation.in   # tensor C_ijkl vía fluctuaciones NVT
```
**Propiedades objetivo:**
- E_CSH ≈ 20-30 GPa, ν_CSH ≈ 0.24-0.30
- E_PE ≈ 0.8-1.2 GPa, ν_PE ≈ 0.40-0.45
- Anisotropía C-S-H (estructura laminar)

### Día 6: Energía de Adhesión Interfacial
```bash
python calculate_work_of_adhesion.py --interface interface_traj.dcd --csh csh_only.dcd --pe pe_only.dcd
```
**Output crítico:** W_adhesión (J/m²) — este valor alimentará directamente el parámetro de cohesión en DEM, cerrando la trazabilidad física de punta a punta.

---

## FASE 2: MAPEO MD → DPD (Días 7-8)

### Día 7: Coarse-Graining Sistemático
**Mapeo:**
```
C-S-H: 1 bead ≈ 10-15 átomos (cluster Ca-Si-O)
PE: 1 bead ≈ 4-6 átomos CH₂ (segmento Kuhn)
```
```python
python derive_dpd_parameters.py --md_energies md_results.json --mapping mapping_scheme.json
# Output: a_CSH-CSH, a_PE-PE, a_CSH-PE derivados de energías de exceso MD
```

### Día 8: Validación Preliminar del Mapeo
- [ ] Comparar densidad reducida DPD objetivo vs densidad atómica MD
- [ ] Verificar compresibilidad isotérmica consistente
- [ ] Ajustar si el error de mapeo > 10%

---

## FASE 3: ESCALA MESOSCÓPICA — DPD, ESTADO GEL (Días 9-15)

### Día 9-10: Setup Simulaciones DPD
```bash
# LAMMPS (pair_style dpd) o alternativamente HOOMD-blue (GPU)
lmp -in dpd_gel_10PE.in
lmp -in dpd_gel_30PE.in
lmp -in dpd_gel_50PE.in
```
**Nota sobre HOOMD-blue:** si el tamaño de sistema crece (>10⁶ partículas) o se quiere aprovechar GPU, HOOMD-blue es preferible a LAMMPS para esta etapa — su motor DPD está optimizado para GPU y facilita correr múltiples réplicas en paralelo, algo útil para el análisis estadístico de percolación que viene después.

### Día 11-13: Evolución del Gel y Formación de Red
**Simulación larga (NVT-DPD, conservando hidrodinámica):**
```bash
# 1-2M steps, monitoreo continuo de conectividad
fix output all ave/time 1000 1 1000 c_cluster_count
```
```python
python monitor_percolation.py --trajectory gel_evolution.dcd --interval 10000
```
**Se registra en cada intervalo:**
- Fracción de clusters C-S-H conectados
- Tamaño del cluster más grande (spanning cluster)
- Viscosidad efectiva del sistema (proxy de fraguado)

### Día 14: Criterio de Transición Gel→Sólido
```python
python detect_percolation_threshold.py --cluster_data cluster_evolution.json
# Criterio: cluster que abarca >90% de la caja en al menos una dimensión
# Validación cruzada: salto en viscosidad efectiva coincide con percolación
```
**Este es el paso metodológicamente más delicado del proyecto** — se recomienda documentar al menos dos definiciones alternativas del umbral (topológica vs reológica) y verificar que ambas convergen a un punto similar en el tiempo de simulación, para blindar el criterio ante revisores.

### Día 15: Extracción de Morfología en el Punto de Transición
```python
python extract_gel_morphology.py --snapshot t_percolation.dcd --output dem_seed_geometry.json
```

---

## FASE 4: HOMOGENEIZACIÓN DPD → DEM (Días 16-17)

### Día 16: Mapeo de Propiedades Constitutivas
```python
python homogenize_dpd_to_dem.py --dpd_morphology dem_seed_geometry.json --md_adhesion W_adhesion.json
```
**Derivación:**
```
kn (C-S-H/C-S-H) desde E_CSH (MD) + geometría local (DPD)
kn (PE/PE) desde E_PE (MD)
cohesión (C-S-H/PE) desde W_adhesión (MD) directamente — no re-derivada en DPD
μ estimado desde morfología de contacto en DPD
```
Este es el punto donde se cierra la cadena de trazabilidad: los parámetros de contacto DEM no se inventan ni se calibran solo contra DPD, sino que heredan valores físicos desde MD.

### Día 17: Generación de Geometría DEM
```python
python generate_dem_packing.py --morphology dem_seed_geometry.json --output dem_geometry.data
```

---

## FASE 5: ESCALA MACROSCÓPICA — DEM, ESTADO SÓLIDO (Días 18-24)

### Día 18-19: Setup DEM
```bash
# LIGGGHTS (recomendado, mayor soporte cohesión/contacto avanzado)
liggghts -in dem_composite_setup.in
```
```
pair_style gran/hertz/history cohesion sjkr
pair_coeff 1 1 kn_csh kt_csh mu_csh 0
pair_coeff 2 2 kn_pe  kt_pe  mu_pe  0
pair_coeff 1 2 kn_int kt_int mu_int cohesion_energy
```
**Nota sobre YADE:** si se prefiere control fino vía scripting Python nativo (útil para automatizar el barrido de sensibilidad del criterio de transición, o si el equipo ya tiene experiencia en Python sobre C++), YADE es una alternativa sólida a LIGGGHTS — con la salvedad de que su comunidad y documentación de modelos cohesivos avanzados es algo menor.

### Día 20-21: Ensayos Virtuales
```bash
liggghts -in dem_compression.in     # uniaxial
liggghts -in dem_tension_brazilian.in
liggghts -in dem_triaxial.in
liggghts -in dem_flexure_3pt.in
```

### Día 22-23: Análisis de Fractura y Damage
```python
python analyze_dem_fracture.py --dump_files compression_*.dump
python classify_fracture_modes.py --damage_data damage_zones.json
```

### Día 24: Validación Trazable de Punta a Punta
```python
python validate_full_chain.py --md_props md_results.json --dpd_morphology dem_seed_geometry.json --dem_results dem_mechanical.json
```
**Se reporta:**
- Consistencia E_CSH(MD) vs rigidez efectiva reconstruida en DEM
- Consistencia W_adhesión(MD) vs comportamiento cohesivo/fractura interfacial en DEM
- Sensibilidad del resultado final ante variaciones del criterio de transición (±10% en umbral de percolación)

---

## FASE 6: ANÁLISIS INTEGRAL Y REDACCIÓN (Días 25-32)

### Día 25-26: Análisis de Sensibilidad y Escalabilidad
```python
python sensitivity_analysis_transition.py   # sensibilidad del umbral gel-sólido
python scalability_three_scales.py          # costo computacional MD/DPD/DEM
```

### Día 27-28: Figuras Científicas
1. Esquema multiescala completo (átomo → gel → sólido) con escalas espaciotemporales
2. Propiedades MD (curvas σ-ε, anisotropía C-S-H)
3. Evolución DPD del gel + curva de percolación con umbral marcado
4. Geometría DEM derivada + ensayos virtuales (σ-ε, fractura)
5. Validación trazable de punta a punta (correlación MD→DPD→DEM)

### Día 29-30: Tablas
- Parámetros MD (force fields, condiciones, propiedades extraídas)
- Parámetros DPD derivados (mapeo, a_ij)
- Criterio de transición (definiciones, valores de umbral, sensibilidad)
- Parámetros DEM homogeneizados y resultados de validación

### Día 31-32: Redacción y Revisión
- [ ] Abstract, Introduction, Methodology (con subsección dedicada al criterio de transición)
- [ ] Results (MD, DPD/percolación, DEM, validación trazable)
- [ ] Conclusions, referencias, formato JCP
- [ ] Revisión final

---

## STACK DE SOFTWARE

### Escala atomística (MD)
- **LAMMPS** — motor principal, force fields CLAYFF + OPLS-AA
- **VMD / OVITO** — visualización y análisis estructural

### Escala mesoscópica (DPD)
- **LAMMPS** (`pair_style dpd`) — continuidad directa con la etapa MD, mismo formato de datos
- **HOOMD-blue** — alternativa recomendada si el sistema crece o se necesitan muchas réplicas en paralelo (GPU-optimizado), útil específicamente para el barrido estadístico del punto de percolación

### Escala macroscópica (DEM)
- **LIGGGHTS** — motor principal, mejor soporte para modelos cohesivos (SJKR) necesarios en la interfaz C-S-H/PE
- **YADE** — alternativa si se prioriza scripting Python nativo y control fino sobre la geometría/condiciones iniciales, con documentación algo más limitada en cohesión avanzada

### Capa de acoplamiento (no estándar, desarrollo propio)
- **Python** (NumPy, SciPy, MDAnalysis, NetworkX) — mapeo MD→DPD, detección de percolación, homogeneización DPD→DEM
- **ParaView** — visualización conjunta de resultados DEM

---

## RECURSOS COMPUTACIONALES

| Etapa | CPU/GPU | RAM | Tiempo estimado |
|---|---|---|---|
| MD (LAMMPS) | 32-64 cores | 64 GB | ~800-1200 core-hours |
| DPD (LAMMPS/HOOMD-GPU) | 32 cores o 1-2 GPU | 64-128 GB | ~500-800 core-hours |
| DEM (LIGGGHTS/YADE) | 16-32 cores | 32-64 GB | ~200-400 core-hours |
| **Total** | | | **~1500-2400 core-hours** |
| Almacenamiento | | | ~250-300 GB (trayectorias + análisis) |

---

## CRONOGRAMA RESUMIDO (5-6 semanas)

| Semana | Actividad principal |
|---|---|
| 1 | MD: construcción, equilibración, propiedades mecánicas base |
| 1.5 | Mapeo MD→DPD |
| 2 | DPD: evolución del gel, monitoreo de percolación |
| 2.5 | Criterio de transición + homogeneización DPD→DEM |
| 3 | DEM: setup, ensayos virtuales |
| 3.5 | Validación trazable MD→DPD→DEM |
| 4 | Análisis de sensibilidad y escalabilidad |
| 4.5-5 | Figuras, tablas, redacción |

**Nota sobre el tiempo:** si el mes original es una restricción estricta, la variable de ajuste más defendible es reducir composiciones PE evaluadas (de 3 a 2) y correr una sola réplica de percolación en vez de un ensemble estadístico completo — documentando esa limitación explícitamente en el paper, en vez de recortar la etapa MD, que es la que sostiene la trazabilidad física de todo el argumento.