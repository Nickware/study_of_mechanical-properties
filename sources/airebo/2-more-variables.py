from lammps_logfile import File
import matplotlib.pyplot as plt

# Cargar archivo log
log = File("log.lammps")

# Extraer datos
steps = log.get("Step")
temp = log.get("Temp")
epair = log.get("E_pair")
press = log.get("Press")

# Configurar figura con 3 subplots verticales
fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
plt.suptitle("Análisis de simulación AIREBO", y=0.98)

# Gráfica 1: Temperatura vs Paso
axs[0].plot(steps, temp, 'r-', linewidth=1.5)
axs[0].set_ylabel('Temperatura (K)', fontsize=10)
axs[0].grid(True, linestyle='--', alpha=0.7)
axs[0].set_title('Evolución térmica', fontsize=12)

# Gráfica 2: Energía Potencial vs Paso
axs[1].plot(steps, epair, 'b-', linewidth=1.5)
axs[1].set_ylabel('Energía Potencial (eV)', fontsize=10)
axs[1].grid(True, linestyle='--', alpha=0.7)
axs[1].set_title('Dinámica energética', fontsize=12)

# Gráfica 3: Presión vs Paso
axs[2].plot(steps, press, 'g-', linewidth=1.5)
axs[2].set_xlabel('Paso de simulación', fontsize=10)
axs[2].set_ylabel('Presión (bar)', fontsize=10)
axs[2].grid(True, linestyle='--', alpha=0.7)
axs[2].set_title('Comportamiento mecánico', fontsize=12)

plt.tight_layout()
plt.show()
