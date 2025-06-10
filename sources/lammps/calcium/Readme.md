# Simulación Mediante Dinámica Molecular de una Estructura de Silicio de Calcio Hidratado

Script de entrada de LAMMPS para simular sistemas que contienen calcio, con parámetros optimizados para estabilidad y análisis.

## Tabla de contenidos
- [Descripción general de la simulación](#descripción-general-de-la-simulación)
- [Características principales](#características-principales)
- [Estructura del archivo de entrada](#estructura-del-archivo-de-entrada)
- [Cómo ejecutar](#cómo-ejecutar)
- [Archivos de salida](#archivos-de-salida)
- [Consejos para personalizar](#consejos-para-personalizar)

## Descripción general de la simulación

Este script realiza simulaciones de dinámica molecular de sistemas que contienen calcio (probablemente hidratos de silicato de calcio) usando:
- **Potenciales LJ/Coulomb** con mezcla geométrica
- **Parámetros de campo de fuerza modificados** para mayor estabilidad
- **Protocolo en varias etapas**: minimización de energía → equilibración → producción
- **Análisis de función de distribución radial (RDF)**

## Características principales

### Configuración

- **Unidades reales** para simulaciones de materiales
- **Condiciones de frontera periódicas en 3D**
- **Paralelización OpenMP** (4 hilos)

### Interacciones

- **Potenciales de par**: LJ (corte=10Å) + Coulomb
- **Potenciales armónicos** para enlaces y ángulos

### Parámetros clave

| Tipo de interacción | ε (kcal/mol) | σ (Å) | Notas                |
|---------------------|--------------|-------|----------------------|
| Ca-Ca (1-1)         | 0.05         | 3.55  | 30% del ε original   |
| Ca-Ow (1-3)         | 0.1165       | 3.55  |                      |
| Si-Ow (2-3)         | 0.1165       | 3.55  |                      |

**Parámetros de enlaces/ángulos:**
- **Rigidez del enlace:** 100 kcal/mol/Å²
- **Rigidez del ángulo:** 15 kcal/mol/rad²

## Cómo ejecutar

### Requisitos
- LAMMPS compilado con soporte OpenMP
- Archivo de entrada `calcium.data`

### Comando de ejecución


## Archivos de salida

1. **`log.lammps`**: Registro de la simulación con datos termodinámicos
2. **`rdf.dat`**: Datos de la función de distribución radial (RDF)
   - Columnas: Distancia (Å), g(r) para todos los pares de átomos
3. **Archivos de reinicio** (si están configurados)

## Consejos para personalizar

1. **Tamaño del sistema**: Modificar `calcium.data` para diferentes configuraciones
2. **Control de temperatura**:
3. **Duración de la simulación**:
4. **Resolución del RDF**:

## Solución de problemas

- **Fallos en la minimización de energía**: Se debe revisar la estructura inicial en `calcium.data`
- **Problemas de inestabilidad**: Reducir el paso de tiempo (actualmente 0.1 fs)
- **Muestreo insuficiente**: Aumentar los pasos de equilibración/producción

**Nota**: Este script utiliza parámetros de campo de fuerza modificados para optimizar la estabilidad del sistema. La segunda fase implica obtener propiedades materiales precisas, validar los resultados con datos experimentales o cálculos con métodos cuánticos.


