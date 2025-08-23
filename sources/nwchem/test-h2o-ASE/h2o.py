from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.io import write

# Crear molécula de agua
h2o = Atoms('H2O', positions=[(0, 0, 0), (0.7, 0.7, 0), (-0.7, 0.7, 0)])
h2o.center(vacuum=5.0)

# Configurar calculador: usar 'gradient' para obtener fuerzas
h2o.calc = NWChem(xc='pbe', basis='3-21g', task='gradient', label='h2o')

# Optimizar geometría
opt = BFGS(h2o)
opt.run(fmax=0.02)

# Guardar estructura final
write('H2O.xyz', h2o)

# Obtener y mostrar energía
energy = h2o.get_potential_energy()
print(f"Energía potencial final: {energy:.6f} eV")