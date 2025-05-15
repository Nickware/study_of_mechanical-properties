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
