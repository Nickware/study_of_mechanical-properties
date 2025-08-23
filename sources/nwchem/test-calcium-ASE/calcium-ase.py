# optimizar_CSH_mejorado.py
# Optimización DFT de tobermorita 11Å con ASE + NWChem

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

# Asegúrate de tener el archivo 'tobermorite_11A.cif' en el mismo directorio
if not os.path.exists('tobermorite_11A.cif'):
    raise FileNotFoundError("Archivo 'tobermorite_11A.cif' no encontrado. Asegúrate de guardarlo en la carpeta actual.")

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
basis = {
    'Ca': 'lanl2dz',
    'Si': '6-31g*',
    'O': '6-31g*',
    'H': '6-31g*'
}

ecp = {'Ca': 'lanl2dz'}  # Potencial de core efectivo para calcio

# --- Funcional DFT con corrección de dispersión ---
# NWChem soporta DFT-D3(BJ) con 'dispersion d3bj'

calc = NWChem(
    label='csh_opt_d3',
    theory='dft',
    xc='pbe',                     # GGA-PBE
    basis=basis,
    ecp=ecp,
    dispersion='d3bj',            # Corrección de dispersión Grimme D3 con Becke-Johnson
    mult=1,                       # Estado singlete (supuesto para sistema cerrado)
    charge=0,                     # Carga total: ajustada a 0 para neutralidad global
    conv=1e-7,                    # Alta convergencia
    maxiter=200,
    direct=True,
    grid='xfine',                 # Malla de integración muy fina
    task='optimize',              # Optimización geométrica
    # Opcional: aumentar precisión de SCF
    damp=False,                   # Evita problemas de convergencia
    precision='double'
)

# Asignar calculador
atoms.calc = calc

# =============================================================================
# Paso 3: Optimización de geometría
# =============================================================================

print("\nIniciando optimización de geometría con NWChem...")
print("Funcional: PBE + D3(BJ), Base: 6-31g*/LANL2DZ, Grid: xfine")

optimizer = BFGS(
    atoms,
    trajectory='csh_opt.traj',   # Guarda camino de optimización
    logfile='csh_opt.log'
)

# Critério estricto: fuerza máxima < 0.05 eV/Å
optimizer.run(fmax=0.05)

print("✅ Optimización finalizada.")

# =============================================================================
# Paso 4: Guardar resultados y análisis
# =============================================================================

# Guardar estructura optimizada
write('csh_optimizado.xyz', atoms)
write('csh_optimizado.cif', atoms)
write('csh_optimizado.json', atoms)  # Formato intermedio para otros scripts

# Energía potencial final
energy = atoms.get_potential_energy()
print(f"\n🔍 Resultados finales:")
print(f"Energía potencial final: {energy:.6f} eV")

# Fuerzas máximas y RMS
forces = atoms.get_forces()
max_force = np.max(np.abs(forces))
rms_force = np.sqrt(np.mean(forces ** 2))
print(f"Fuerza máxima: {max_force:.6f} eV/Å")
print(f"Fuerza RMS: {rms_force:.6f} eV/Å")

# Verificación de convergencia
if max_force < 0.05:
    print("✅ Convergencia alcanzada (fmax < 0.05 eV/Å)")
else:
    print("⚠️ Advertencia: convergencia incompleta.")

# Información adicional
print(f"\n📊 Estadísticas:")
print(f"Total de átomos: {len(atoms)}")
print(f"Fórmula química: {atoms.get_chemical_formula()}")

# Opcional: calcular temperatura electrónica (no crítica)
e_kin = 0.5 * np.sum(forces**2) / len(atoms)  # Solo ilustrativo
print(f"Energía ")


# =============================================================================
# Paso 5: Recomendaciones post-optimización
# =============================================================================

print("\n📌 Recomendaciones:")
print("1. Verifica la geometría optimizada visualmente (con VMD, Ovito o ASE-GUI).")
print("2. Considera un cálculo single-point con funcional híbrido (PBE0) o base mayor.")
print("3. Calcula propiedades: band gap, DOS, funciones de pair distribution (PDF).")
print("4. Prueba con superceldas (2x2x1) para reducir interacciones espurias.")
print("5. Compara energía por fórmula con otras fases (jennite, portlandite).")

# Ejemplo: abrir con ASE-GUI (descomenta si lo deseas)
# from ase.visualize import view
# view(atoms)