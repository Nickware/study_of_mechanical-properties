# test_water.py
from ase import Atoms
import subprocess

h2o = Atoms('H2O', positions=[[0, 0, 0], [0.7, 0.7, 0], [-0.7, 0.7, 0]])
h2o.center(vacuum=5.0)

input = """start h2o
geometry
  O  0.0  0.0  0.0
  H  0.7  0.7  0.0
  H -0.7  0.7  0.0
end
basis
  * library 3-21g
end
dft
  xc pbe96
  mult 1
end
task dft energy
"""

with open('h2o.nwi', 'w') as f:
    f.write(input)

subprocess.run(['nwchem', 'h2o.nwi'])
