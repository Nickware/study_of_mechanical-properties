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
from ase.calculators.nwchem import NWChem
import numpy as np
import os

# =============================================================================
# Paso 1: Leer estructura desde CIF
# =============================================================================

cif_filename = 'tobermorite_11A.cif'
if not os.path.exists(cif_filename):
    raise FileNotFoundError(f"Archivo '{cif_filename}' no encontrado.")

print(f"Leyendo estructura desde {cif_filename}...")
atoms = read(cif_filename)

# Centrar y aplicar condiciones periódicas
atoms.center()
atoms.set_pbc(True)

print(f"Estructura cargada: {atoms.get_chemical_formula()} ({len(atoms)} átomos)")

# Guardar estructuras previas
write('csh_inicial.xyz', atoms)

# =============================================================================
# Paso 2: Configurar calculador NWChem (Inyección directa de sintaxis)
# =============================================================================

basis = '3-21g'

# Al pasar las directivas como una lista de strings, ASE se limita a imprimirlas
# línea por línea dentro del bloque 'dft ... end', eliminando los errores de parsing.
dft_raw_directives = [
    'xc pbe96',
    'mult 1',
    'iterations 150',
    'convergence energy 1e-6 density 1e-5',
    'iterations damping on',
    'smear 0.001',
    'disp vdw 3',
    'direct',
    'grid fine'
]

calc = NWChem(
    label='csh_opt_d3',
    theory='dft',
    basis=basis,
    charge=0,
    task='energy',         # ASE controlará las fuerzas paso a paso en el BFGS
    dft=dft_raw_directives
)

atoms.calc = calc

# =============================================================================
# Paso 3: Optimización de geometría iónica mediante ASE
# =============================================================================

print("\nIniciando optimización de geometría con NWChem (Vía BFGS de ASE)...")
print("Configuración: DFT-PBE + Grimme D3, Base: 3-21g, Grid: fine")

optimizer = BFGS(
    atoms,
    trajectory='csh_opt.traj',
    logfile='csh_opt.log'
)

# Ejecución del optimizador
optimizer.run(fmax=0.05)

print("\n Optimización finalizada.")

# =============================================================================
# Paso 4: Guardar resultados y análisis de fuerzas
# =============================================================================

write('csh_optimizado.xyz', atoms)
write('csh_optimizado.cif', atoms)

energy = atoms.get_potential_energy()
forces = atoms.get_forces()
max_force = np.max(np.abs(forces))

print(f"\n--- Resultados Finales ---")
print(f"Energía potencial final : {energy:.6f} eV")
print(f"Fuerza máxima registrada : {max_force:.6f} eV/Å")