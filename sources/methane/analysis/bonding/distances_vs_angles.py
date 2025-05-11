import matplotlib.pyplot as plt
import numpy as np

# Datos
distances = [1.1307, 1.0700, 1.1487, 1.1144]
angles = [115.93, 104.92, 103.78, 110.14, 108.21, 113.88]

# Gráfico
plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
plt.bar([f'C-H{i+1}' for i in range(4)], distances, color='skyblue')
plt.axhline(1.09, color='r', linestyle='--')
plt.title('Distancias C-H')
plt.ylabel('Å')

plt.subplot(1,2,2)
plt.hist(angles, bins=6, color='lightgreen', edgecolor='black')
plt.axvline(109.47, color='r', linestyle='--')
plt.title('Distribución de Ángulos H-C-H')
plt.xlabel('Grados')

plt.tight_layout()
plt.show()