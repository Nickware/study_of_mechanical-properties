# Análisis de temperatura y energía en simulación AIREBO
# Este script lee el archivo de log generado por LAMMPS durante una simulación AIREBO, extrae los datos de temperatura y energía, y luego los visualiza y presenta en una tabla.
# Asegúrate de tener el archivo "log.airebo" en el mismo directorio que este script o proporciona la ruta correcta al archivo.
# Requisitos:
# - lammps_logfile: Para leer los archivos de log de LAMMPS
# - matplotlib: Para graficar los datos
# - pandas: Para crear y mostrar la tabla de datos
# Puedes instalar las bibliotecas necesarias usando pip:
# pip install lammps_logfile matplotlib pandas
# Importar las bibliotecas necesarias
# Asegúrate de tener las bibliotecas instaladas antes de ejecutar este script
# pip install lammps_logfile matplotlib pandas
from lammps_logfile import File
import matplotlib.pyplot as plt

log = File("log.airebo")  # Nombre de tu archivo log

# Extraer datos termodinámicos
steps = log.get("Step")
temp = log.get("Temp")
epair = log.get("E_pair")

# Gráfica de temperatura vs tiempo
plt.figure(figsize=(10,5))
plt.plot(steps, temp, 'r-', label='Temperatura (K)')
plt.xlabel('Paso de simulación')
plt.ylabel('Temperatura (K)')
plt.legend()
plt.grid(True)
plt.show()

# Tabla de energía vs paso
import pandas as pd
df = pd.DataFrame({
    'Paso': steps,
    'Energía Potencial (eV)': epair,
    'Energía Total (eV)': log.get("TotEng")
})
print(df.head())
