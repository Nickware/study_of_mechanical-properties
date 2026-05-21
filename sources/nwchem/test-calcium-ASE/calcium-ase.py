# optimizar_CSH_mejorado.py
# Optimización DFT de tobermorita 11Å con ASE + NWChem
# Este script realiza una optimización geométrica de la estructura de tobermorita 11Å utilizando el funcional PBE con corrección de dispersión D3(BJ) y una base 3-21g. Se han ajustado parámetros para mejorar la convergencia y obtener resultados más realistas.
# Asegúrate de tener ASE y NWChem configurados correctamente en tu entorno para
# ejecutar este script. El archivo 'tobermorite_11A.cif' debe estar en el mismo directorio que este script para cargar la estructura inicial.
# El script guarda la estructura optimizada en formatos XYZ, CIF y JSON, y proporciona información sobre la energía final, fuerzas y recomendaciones para análisis posteriores. La optimización se realiza con un criterio de fuerza máxima de 0.05 eV/Å para asegurar una geometría razonablemente optimizada.
# Nota: La geometría inicial y los parámetros de cálculo pueden necesitar ajustes dependiendo de los resultados obtenidos y la convergencia del cálculo. Es recomendable revisar la estructura optimizada visualmente y considerar cálculos adicionales para validar los resultados.
# Importamos ASE para manejar la geometría y generar el archivo de entrada para NWChem
# ASE es una herramienta útil para crear y manipular estructuras atómicas, así como para escribir
# archivos de entrada en diferentes formatos. En este caso, lo usamos para generar un archivo XYZ con la geometría inicial y para facilitar la creación del archivo de entrada de NWChem.
# Asegúrate de tener ASE instalado en tu entorno de Python para ejecutar este script.
# Puedes instalarlo usando pip:
# pip install ase
# El archivo de entrada para NWChem se configura para realizar una optimización geométrica utilizando DFT con el funcional PBE y corrección de dispersión D3(BJ), junto con una base 3-21g. La multiplicidad se establece en 1 para simular un estado singlete, lo que es común en sistemas cerrados. Se han ajustado parámetros como la malla de integración y el damping para mejorar la convergencia del cálculo.
# El cluster generado es una aproximación y no representa una estructura cristalina real, pero es más complejo que un cluster simple. La geometría inicial puede necesitar ajustes dependiendo de los resultados obtenidos y la convergencia del cálculo.
# El script también incluye recomendaciones para análisis posteriores, como la verificación visual de la geometría optimizada, cálculos single-point con funcionales híbridos o bases mayores, y el cálculo de propiedades adicionales como el band gap o funciones de pair distribution (PDF).
# Nota: Asegúrate de revisar los resultados obtenidos y ajustar los parámetros de cálculo según sea necesario para obtener una optimización exitosa y resultados realistas.
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from ase.units import kB
from ase.calculators.nwchem import NWChem
import numpy as np
import os

# =============================================================================
# Paso 1: Leer estructura desde CIF
# =============================================================================

# Asegurarse de tener el archivo 'tobermorite_11A.cif' en el mismo directorio
if not os.path.exists('tobermorite_11A.cif'):
    raise FileNotFoundError("Archivo 'tobermorite_11A.cif' no encontrado. Asegurarse de guardarlo en la carpeta actual.")

print("Leyendo estructura desde tobermorite_11A.cif...")
atoms = read('tobermorite_11A.cif')

# Opcional: crear una supercelda para reducir efectos de periodicidad
# atoms = atoms * (2, 1, 1)  # Ejemplo: duplicar en a

# Centrar y aplicar condiciones periódicas
atoms.center()
atoms.set_pbc(True)

print(f"Estructura cargada: {atoms.get_chemical_formula()} ({len(atoms)} átomos)")
print(f"Celda: a={atoms.cell.lengths()[0]:.3f}, b={atoms.cell.lengths()[1]:.3f}, c={atoms.cell.lengths()[2]:.3f} Å")
print(f"Ángulos: α={atoms.cell.angles()[0]:.2f}, β={atoms.cell.angles()[1]:.2f}, γ={atoms.cell.angles()[2]:.2f}°")

# Guardar estructura inicial
write('csh_inicial.xyz', atoms)
write('csh_inicial.cif', atoms)

# =============================================================================
# Paso 2: Configurar calculador NWChem con DFT mejorado
# =============================================================================

# --- Análisis de carga total ---
# Fórmula aproximada: Ca₃Si₃O₁₂H₉ → Ca₃²⁺ = +6, Si₃⁴⁺ = +12, O²⁻/OH/H₂O ≈ -24 +3 +0 → carga neta ≈ -3?
# En estructuras cristalinas, la neutralidad se mantiene por ocupación y coordinación.
# Para simulación, mejor partir de **carga 0** y permitir redistribución electrónica.
# O usar carga basada en estequiometría formal. Aquí usamos **charge=0** como punto inicial.

# --- Bases y ECP ---
# Usar bases compatibles con NWChem. Para sistemas grandes, usar bases pequeñas inicialmente
# 3-21g es más manejable que 6-31g* para 208 átomos
basis = '3-21g'  # Base uniforme para todos

# --- Funcional DFT con corrección de dispersión ---
# NWChem soporta DFT-D3(BJ) con 'dispersion d3bj'
# Para sistemas grandes, usar grid estándar en lugar de xfine

calc = NWChem(
    label='csh_opt_d3',
    theory='dft',
    xc='pbe',                     # GGA-PBE
    basis=basis,
    dispersion='d3bj',            # Corrección de dispersión Grimme D3 con Becke-Johnson
    mult=1,                       # Estado singlete (supuesto para sistema cerrado)
    charge=0,                     # Carga total: ajustada a 0 para neutralidad global
    convergence={'energy': 1e-6, 'density': 1e-5},  # Convergencia DFT
    maxiter=150,
    direct=True,
    grid='fine',                  # Malla de integración 'fine' en lugar de 'xfine' (más estable)
    task='optimize',              # Optimización geométrica
    damp=True,                    # Activar damping para mejorar convergencia
    smear=0.001                   # Broadening pequeño para ayudar SCF
)

# Asignar calculador
atoms.calc = calc

# =============================================================================
# Paso 3: Optimización de geometría
# =============================================================================

print("\nIniciando optimización de geometría con NWChem...")
print("Funcional: PBE + D3(BJ), Base: 3-21g, Grid: fine")

optimizer = BFGS(
    atoms,
    trajectory='csh_opt.traj',   # Guarda camino de optimización
    logfile='csh_opt.log'
)

# Critério estricto: fuerza máxima < 0.05 eV/Å
optimizer.run(fmax=0.05)

print(" Optimización finalizada.")

# =============================================================================
# Paso 4: Guardar resultados y análisis
# =============================================================================

# Guardar estructura optimizada
write('csh_optimizado.xyz', atoms)
write('csh_optimizado.cif', atoms)
write('csh_optimizado.json', atoms)  # Formato intermedio para otros scripts

# Energía potencial final
energy = atoms.get_potential_energy()
print(f"\n Resultados finales:")
print(f"Energía potencial final: {energy:.6f} eV")

# Fuerzas máximas y RMS
forces = atoms.get_forces()
max_force = np.max(np.abs(forces))
rms_force = np.sqrt(np.mean(forces ** 2))
print(f"Fuerza máxima: {max_force:.6f} eV/Å")
print(f"Fuerza RMS: {rms_force:.6f} eV/Å")

# Verificación de convergencia
if max_force < 0.05:
    print(" Convergencia alcanzada (fmax < 0.05 eV/Å)")
else:
    print(" Advertencia: convergencia incompleta.")

# Información adicional
print(f"\n Estadísticas:")
print(f"Total de átomos: {len(atoms)}")
print(f"Fórmula química: {atoms.get_chemical_formula()}")

# Opcional: calcular temperatura electrónica (no crítica)
e_kin = 0.5 * np.sum(forces**2) / len(atoms)  # Solo ilustrativo
print(f"Energía ")


# =============================================================================
# Paso 5: Recomendaciones post-optimización
# =============================================================================

print("\n Recomendaciones:")
print("1. Verificar la geometría optimizada visualmente (con VMD, Ovito o ASE-GUI).")
print("2. Considerar un cálculo single-point con funcional híbrido (PBE0) o base mayor.")
print("3. Calcular propiedades: band gap, DOS, funciones de pair distribution (PDF).")
print("4. Probar con superceldas (2x2x1) para reducir interacciones espurias.")
print("5. Comparar energía por fórmula con otras fases (jennite, portlandite).")

# Ejemplo: abrir con ASE-GUI (descomenta si lo deseas)
# from ase.visualize import view
# view(atoms)
