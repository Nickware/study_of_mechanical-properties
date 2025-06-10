# Simulación mediante Dinámica Molecular de la Estructura del Metano 

Script de LAMMPS para simular moléculas de metano (CH₄) con parámetros optimizados para estabilidad y análisis estructural.

## Tabla de contenidos
- [Descripción general](#descripción-general)
- [Características clave](#características-clave)
- [Requisitos](#requisitos)
- [Ejecución](#ejecución)
- [Archivos de salida](#archivos-de-salida)
- [Personalización](#personalización)

## Descripción general
Simulación que modela moléculas de metano usando:
- **Potenciales armónicos** para enlaces y ángulos
- **Interacciones LJ** con parámetros específicos
- **Protocolo multi-etapa**: minimización de energía → equilibración → producción
- **Salida de trayectorias** en formato XYZ

## Características clave

### Configuración principal
- **Unidades reales** para precisión en sistemas moleculares
- **Condiciones de frontera periódicas** en 3D
- **Algoritmo FIRE** para minimización energética

### Parámetros de interacción
| Tipo          | Coeficiente 1 | Coeficiente 2       | Notas               |
|---------------|---------------|---------------------|---------------------|
| Enlace C-H    | 450 kcal/mol/Å² | 1.09 Å (longitud)  | Rigidez alta        |
| Ángulo H-C-H  | 120 kcal/mol/rad² | 109.5° (ángulo)   | Geometría tetraédrica |
| Par C-C       | ε=0.1094 kcal/mol | σ=3.40 Å          | Interacción débil   |
| Par C-H       | ε=0.0380 kcal/mol | σ=2.42 Å          |                     |

### Protocolo de simulación
1. **Minimización energética** (1M pasos máx)
2. **Termalización** a 300 K con termostato NVT
3. **Fase de producción** de 5 ps (10,000 pasos)

## Requisitos
- LAMMPS compilado con soporte para paquete **FIRE**
- Archivo de datos inicial `methane_final.data`
- 4 GB RAM mínimo para sistemas pequeños

## Ejecución

## Archivos de salida
1. **`methane.log`**: Registro termodinámico con:
   - Temperatura
   - Energías (potencial, cinética, total)
   - Presión
2. **`traj_methane.xyz`**: Trayectoria molecular cada 1 ps
3. **`methane_final_restart.data`**: Estado final para reinicios

## Personalización
### Ajustes comunes
1. **Duración de simulación**:
2. **Frecuencia de salida**:
3. **Condiciones termostáticas**:

### Optimización de rendimiento

## Análisis recomendado
1. **Geometría molecular**:
   - Distribución de ángulos H-C-H
   - Longitudes de enlace C-H
2. **Dinámica**:
   - Coeficiente de difusión (MSD)
   - Función de autocorrelación de velocidades

**Nota importante**: Los parámetros LJ usados son específicos para interacciones CH₄-CH₄. Para sistemas mixtos, se debe consultar reglas de combinación ([Documentación LAMMPS](https://docs.lammps.org/pair_modify.html))[5].

