from ase.build import bulk
from ase.io import write
from ase.spacegroup import crystal
import numpy as np

# Hielo Ih (hexagonal) - estructura cristalina del hielo común
# Parámetros aproximados del hielo Ih
a = 4.497  # Å
c = 7.315  # Å

# Crear celda unitaria de hielo (simplificada)
# Nota: para hielo real se necesitan posiciones fraccionarias específicas
ice = crystal(['O', 'H', 'H'], 
              basis=[(0, 0, 0), (0.37, 0.37, 0), (0.63, 0.63, 0)],
              spacegroup=194,  # P6₃/mmc
              cellpar=[a, a, c, 90, 90, 120])

# Replicar a supercelda 2x2x2
ice_super = ice * (2, 2, 2)
ice_super.center()

write('ice_ih.xyz', ice_super)
print(f"Átomos en supercelda: {len(ice_super)}")