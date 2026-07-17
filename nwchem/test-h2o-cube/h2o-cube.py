from ase import Atoms
from ase.build import molecule
from ase.io import write
import numpy as np

# Crear una molécula de agua optimizada (o usar la que ya tienes)
h2o = molecule('H2O')
h2o.center()

# Parámetros de la caja
n_molecules = 8  # 2x2x2
box_size = 10.0  # Å
spacing = 3.5    # Distancia entre moléculas

# Crear sistema con múltiples moléculas
positions = []
symbols = []

# Grid 2x2x2
for ix in range(2):
    for iy in range(2):
        for iz in range(2):
            offset = np.array([ix, iy, iz]) * spacing
            for atom in h2o:
                positions.append(atom.position + offset)
                symbols.append(atom.symbol)

# Crear el sistema completo
system = Atoms(symbols=symbols, positions=positions)

# Definir la celda periódica (caja)
system.set_cell([box_size, box_size, box_size])
system.set_pbc([True, True, True])  # Periodic boundary conditions

# Centrar en la caja
system.center()

write('water_box.xyz', system)
print(f"Total átomos: {len(system)}")
print(f"Densidad aproximada: {n_molecules * 18.015 / (box_size**3 * 1e-24):.3f} g/cm³")