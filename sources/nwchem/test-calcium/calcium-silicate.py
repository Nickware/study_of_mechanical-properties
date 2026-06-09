# calcium.py - Geometría inicial 
# Este script genera una geometría inicial para un cluster de silicato de calcio (CaSiO₃) con hidrógenos añadidos para simular grupos OH y puentes.
# La geometría es una aproximación y no representa una estructura cristalina real, pero es más compleja que un cluster simple.
# El script también genera un archivo de entrada para NWChem con una configuración básica de DFT utilizando el funcional PBE96 y una base 3-21g.
# Para ejecutar la optimización, guarda este script, ejecútalo para generar el archivo de entrada y luego corre NWChem con el comando indicado al final.
# Nota: Asegúrate de tener NWChem instalado y configurado correctamente para ejecutar el cálculo. La geometría inicial puede necesitar ajustes dependiendo de los resultados obtenidos y la convergencia del cálculo.
# Importamos ASE para manejar la geometría y generar el archivo de entrada para NWChem
# ASE es una herramienta útil para crear y manipular estructuras atómicas, así como para escribir archivos de entrada en diferentes formatos. En este caso, lo usamos para generar un archivo XYZ con la geometría inicial y para facilitar la creación del archivo de entrada de NWChem.
# Asegúrate de tener ASE instalado en tu entorno de Python para ejecutar este script. Puedes instalarlo usando pip:
# pip install ase
# El cluster generado es una aproximación y no representa una estructura cristalina real, pero es más complejo que un cluster simple. La geometría inicial puede necesitar ajustes dependiendo de los resultados obtenidos y la convergencia del cálculo.
# El archivo de entrada para NWChem se configura para realizar una optimización geométrica utilizando DFT con el funcional PBE96 y una base 3-21g. La multiplicidad se establece en 2 para simular un estado con un electrón desapareado, lo que es común en sistemas con metales de transición o en casos donde se espera que haya radicales libres.

from ase import Atoms
from ase.io import write

# 1. Definición del clúster C-S-H
positions = [
    # Ca
    [0.0, 0.0, 0.0], [5.0, 5.0, 5.0], [2.0, 6.0, 2.0],        
    # Si
    [1.5, 1.5, 1.5], [4.0, 4.0, 4.0], [6.0, 6.0, 6.0],        
    # O
    [2.3, 1.5, 1.5], [1.5, 2.3, 1.5], [1.5, 1.5, 2.3], [0.7, 1.5, 1.5], 
    [4.8, 4.0, 4.0], [4.0, 4.8, 4.0], [4.0, 4.0, 4.8], [3.2, 4.0, 4.0],
    [6.8, 6.0, 6.0], [6.0, 6.8, 6.0], [0.5, 0.5, 0.5], [5.5, 5.5, 5.5],
    # H
    [2.8, 1.5, 1.5], [1.5, 2.8, 1.5], [1.5, 1.5, 2.8], [0.7, 2.0, 1.0], 
    [4.8, 4.5, 3.5], [4.5, 4.8, 4.5], [4.0, 4.0, 5.3], [3.2, 4.5, 3.5],
    [6.8, 6.5, 5.5]
]

symbols = ['Ca']*3 + ['Si']*3 + ['O']*12 + ['H']*9

# 2. Inicialización y centrado con ASE
atoms = Atoms(symbols=symbols, positions=positions)
atoms.center(vacuum=5.0)  # Traslada los átomos para centrarlos en la caja
write('csh_inicial.xyz', atoms)

# 3. Generación del input para NWChem usando el objeto atoms corregido
nwi = """start csh_opt
geometry units angstrom
"""

for atom in atoms:
    # Usamos atom.position para exportar las coordenadas trasladadas reales
    nwi += f"  {atom.symbol:<2}  {atom.position[0]:12.6f}  {atom.position[1]:12.6f}  {atom.position[2]:12.6f}\n"

# Bloque de directivas NWChem con tratamiento de convergencia para capa abierta
nwi += """end
basis
  * library 3-21g
end
dft
  xc pbe96
  mult 2
  iterations 150
  smear 0.005          # Resuelve la degeneración HOMO-LUMO
  # --- Reemplazo correcto para controlar oscilaciones en DFT ---
  iterations damping on  # Activa la amortiguación adaptativa de NWChem en el SCF
  convergence energy 1e-6
  convergence density 1e-5
end
driver
  maxiter 50
  tight
end
task dft optimize
"""

# 4. Escritura del archivo de entrada
with open('csh_opt.nwi', 'w') as f:
    f.write(nwi)

print(" -> Archivos 'csh_inicial.xyz' y 'csh_opt.nwi' generados con éxito de forma consistente.")
print(" -> Ejecuta en terminal: nwchem csh_opt.nwi > csh_opt.nwo &")
