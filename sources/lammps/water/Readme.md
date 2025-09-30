# Simulación en **LAMMPS** para estudiar agua usando el modelo TIP4P 

Este script configura una simulación en **LAMMPS** para estudiar agua usando el modelo de cuatro puntos TIP4P (un modelo realista y ampliamente utilizado en simulaciones moleculares de agua). El objetivo es simular un sistema de moléculas de agua bajo condiciones cercanas a la temperatura ambiente y explorar propiedades estructurales y termodinámicas del líquido.

***

### Conceptos Clave del Script

- **Unidades y Tipos Atómicos:**  
  Usa `"units real"`, es decir, medidas en Ångströms, kcal/mol, fs, etc.  
  Define dos tipos de átomos: oxígeno (tipo 1) e hidrógeno (tipo 2).

- **Caja y Geometría:**  
  Crea una región cúbica de 10 Å de lado centrada en el origen.  
  Reserva espacio para moléculas de agua incluyendo un número fijo de enlaces (bonds) y ángulos, necesarias para la estructura molecular.

- **Modelo de Interacción:**  
  Usa el **TIP4P**, un modelo que requiere un punto adicional “M” para ubicar la carga negativa fuera del átomo de oxígeno, lo que mejora la representación de la polaridad y el comportamiento del agua.  
  - `pair_style lj/cut/tip4p/cut 1 2 1 1 0.15 8.0`: configura la interacción Lennard-Jones y electrostática según TIP4P, indicando los números para identificar O e H, la distancia OM y el cutoff.
  - Define parámetros LJ entre oxígenos, los hidrógenos solo tienen interacciones de repulsión dura.
  - Los enlaces y ángulos se tratan como rígidos (cero energía de deformación), con distancias y ángulos de equilibrio de TIP3P/TIP4P.

- **Construcción de Moléculas:**  
  - Usa una definición de molécula a partir de un archivo externo (`tip3p.mol`), que describe la geometría de la molécula de agua.
  - Inserta 33 moléculas de forma aleatoria y sin solapamientos excesivos.
  - Asigna cargas parciales a O e H: -1.04 (O) y +0.52 (H), características de TIP4P.

- **Restricciones y Dinámica:**  
  - Utiliza el algoritmo SHAKE para mantener las restricciones de geometría intra-molecular.
  - Minimiza la energía inicial para evitar solapamientos.
  - Inicializa velocidades para simular a 300 K y ejecuta la integración con un termostato NVT.

- **Parámetros de Salida:**  
  - Imprime temperatura, presión, energía total y potencial cada 1000 pasos.
  - Escribe el estado final en un archivo de datos.

***

### Para Qué Sirve Este Script

- **Simulación realista de agua líquida:**  
  El TIP4P es un estándar en simulaciones de agua, ya que reproduce con precisión muchas de sus propiedades físicas y termodinámicas.[1][2]
- **Estudios estructurales y de dinámica:**  
  Permite medir difusión, estructura radial, fluctuaciones de energía y termodinámica a escala atómica.
- **Base para simulaciones más complejas:**  
  Se puede ampliar agregando solutos, superficies, proteínas, o cambiando la presión y temperatura para investigar fenómenos como cristalización, nucleación, solvación, interface agua-material, etc.

***

### Adaptabilidad y Oportunidades

- **Sistemas de Solvatación:**  
  Se puede incluir iones, moléculas orgánicas, o superficies para estudiar procesos bioquímicos o materiales.
- **Cambio de Condiciones:**  
  Modifica fácilmente presión, temperatura, tamaño del sistema o concentración.
- **Extensión a Otros Modelos:**  
  Permite reemplazar TIP4P por otros modelos (TIP3P, OPC, SPC/E) para analizar diferencias en propiedades macroscópicas.
- **Fenómenos Avanzados:**  
  Puede ser la base para analizar transporte de agua en canales, interfaces, simulaciones a presión negativa (cavitación), nucleación de hielo, y otras aplicaciones avanzadas.[3][4][1]

***

### Consideraciones Especiales

- El TIP4P es sensible al manejo de cargas y restricciones geométricas; la selección cuidadosa del termostato y el manejo de SHAKE es crucial para simular correctamente el comportamiento del agua.
- La parametrización y archivos externos deben ser coherentes con el emparejamiento de cargas, distancias y geometría definidas.
