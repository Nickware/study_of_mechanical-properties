# calcium_small.py - ¡ÚLTIMA VERSIÓN QUE FUNCIONA!

from ase import Atoms
from ase.io import write

# -----------------------------------------------------------------------------
# 1. Definir estructura con distancias seguras
# -----------------------------------------------------------------------------

symbols = [
    'Ca', 'Si', 'O', 'O', 'O', 'O',  # SiO₄
    'O', 'O', 'O', 'O',              # OH-Ca
    'H', 'H', 'H', 'H',              # H-OH Ca
    'H', 'H', 'H'                    # H-OH Si
]

positions = [
    [ 0.000,  0.000,  0.000],  # Ca
    [ 3.500,  3.500,  3.500],  # Si
    [ 4.300,  3.500,  3.500],  # O1
    [ 3.500,  4.300,  3.500],  # O2
    [ 3.500,  3.500,  4.300],  # O3
    [ 2.700,  3.500,  3.500],  # O4
    [ 1.200,  0.000,  0.000],  # O5
    [ 0.000,  1.200,  0.000],  # O6
    [ 0.000,  0.000,  1.200],  # O7
    [-1.200,  0.000,  0.000],  # O8
    [ 1.800,  0.000,  0.000],  # H1
    [ 0.000,  1.800,  0.000],  # H2
    [ 0.000,  0.000,  1.800],  # H3
    [-1.800,  0.000,  0.000],  # H4
    [ 4.700,  3.500,  3.500],  # H5
    [ 3.500,  4.700,  3.500],  # H6
    [ 3.500,  3.500,  4.700],  # H7
]

# Verificación
assert len(symbols) == len(positions)
atoms = Atoms(symbols=symbols, positions=positions)
atoms.center(vacuum=5.0)
write('csh_small_inicial.xyz', atoms)

# -----------------------------------------------------------------------------
# 2. Generar input con NOAUTOZ (¡clave!)
# -----------------------------------------------------------------------------

input_lines = [
    "start csh_small",
    "geometry noautoz noautosym units angstrom",  # ← ¡ESencial!
]

for sym, pos in zip(symbols, positions):
    x, y, z = pos
    input_lines.append(f"  {sym}  {x:8.5f}  {y:8.5f}  {z:8.5f}")

input_lines += [
    "end",
    "basis",
    "  * library 3-21g",
    "end",
    "dft",
    "  xc pbe96",
    "  mult 2",           # Por si hay electrones impares
    "  iterations 100",
    "end",
    "driver",
    "  maxiter 30",
    "  loose",
    "end",
    "task dft optimize"
]

# Escribir archivo
with open('csh_small.nwi', 'w') as f:
    f.write("\n".join(input_lines))

print("✅ Archivo 'csh_small.nwi' generado con 'noautoz'.")
print("🚀 Ejecuta: nwchem csh_small.nwi > csh_small.nwo")