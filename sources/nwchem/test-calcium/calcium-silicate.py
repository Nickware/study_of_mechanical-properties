# calcium.py - Geometría inicial más física

from ase import Atoms
from ase.io import write

# Definimos un cluster más realista
# Basado en unidades SiO₄ y Ca²⁺ coordinado a O

positions = [
    # Ca1 - coordinado a 6 O (aprox 2.4 Å)
    [0.0, 0.0, 0.0],  # Ca1
    # Ca2
    [5.0, 5.0, 5.0],  # Ca2
    # Ca3
    [2.0, 6.0, 2.0],  # Ca3

    # Si1 (tetraedro central)
    [1.5, 1.5, 1.5],  # Si1
    # Si2
    [4.0, 4.0, 4.0],  # Si2
    # Si3
    [6.0, 6.0, 6.0],  # Si3

    # O enlazados a Si1
    [2.3, 1.5, 1.5],  # O1  (Si1-O)
    [1.5, 2.3, 1.5],  # O2
    [1.5, 1.5, 2.3],  # O3
    [0.7, 1.5, 1.5],  # O4

    # O enlazados a Si2
    [4.8, 4.0, 4.0],  # O5
    [4.0, 4.8, 4.0],  # O6
    [4.0, 4.0, 4.8],  # O7
    [3.2, 4.0, 4.0],  # O8

    # O enlazados a Si3
    [6.8, 6.0, 6.0],  # O9
    [6.0, 6.8, 6.0],  # O10

    # OH y puentes
    [0.5, 0.5, 0.5],  # O11 (puente Ca1-O-Ca3)
    [5.5, 5.5, 5.5],  # O12 (puente Ca2-O-Si3)

    # H unidos a O (ángulo ~109°, distancia ~0.96 Å)
    [2.8, 1.5, 1.5],  # H1  (en O1)
    [1.5, 2.8, 1.5],  # H2  (en O2)
    [1.5, 1.5, 2.8],  # H3  (en O3)
    [0.7, 2.0, 1.0],  # H4  (en O4, OH)
    [4.8, 4.5, 3.5],  # H5  (en O5)
    [4.5, 4.8, 4.5],  # H6  (en O6)
    [4.0, 4.0, 5.3],  # H7  (en O7)
    [3.2, 4.5, 3.5],  # H8  (en O8)
    [6.8, 6.5, 5.5],  # H9  (en O9)
]

symbols = ['Ca']*3 + ['Si']*3 + ['O']*10 + ['H']*9

atoms = Atoms(symbols=symbols, positions=positions)
atoms.center(vacuum=5.0)  # Espacio alrededor
write('csh_inicial.xyz', atoms)

# Generar input limpio
nwi = """start csh_opt
geometry units angstrom
"""

for s, p in zip(symbols, positions):
    nwi += f"  {s}  {p[0]:.6f}  {p[1]:.6f}  {p[2]:.6f}\n"

nwi += """end
basis
  * library 3-21g
end
dft
  xc pbe96
  mult 2
  iterations 100
  convergence energy 1e-6
  convergence density 1e-5
end
driver
  maxiter 50
  tight
end
task dft optimize
"""

with open('csh_opt.nwi', 'w') as f:
    f.write(nwi)

print("✅ Input generado. Ejecuta: nwchem csh_opt.nwi > csh_opt.nwo")