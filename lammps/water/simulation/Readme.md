# TIP4P implícito en LAMMPS

Este script configura una simulación de agua usando un modelo **TIP4P implícito** en LAMMPS, partiendo de una geometría TIP3P pero cambiando cargas y potenciales para comportarse como TIP4P y produciendo un archivo de datos equilibrado (`tip4p-implicit.data`). [docs.lammps](https://docs.lammps.org/Howto_tip4p.html)

## Modelo físico

- Se definen **dos tipos de átomos**: tipo 1 = O, tipo 2 = H, con masas de agua estándar (O: 15.9994, H: 1.008).  
- El comando `pair_style lj/cut/tip4p/cut 1 2 1 1 0.15 8.0` indica el uso del modelo TIP4P implícito, donde la carga negativa no está en el oxígeno sino en un sitio ficticio (M) a una distancia OM = 0.15 Å sobre el bisector del ángulo HOH. [afs.enea](https://www.afs.enea.it/software/lammps/doc19/html/Howto_tip4p.html)
- `pair_coeff` asigna los parámetros de Lennard-Jones al oxígeno (ε = 0.1550, σ = 3.1536) y desactiva LJ en el hidrógeno, como es típico en TIP3P/TIP4P. [docs.lammps](https://docs.lammps.org/Howto_tip4p.html)
- Se usa la geometría de **TIP3P** desde `tip3p.mol`, pero luego se ajustan las **cargas** para TIP4P: O = −1.040, H = +0.520, de forma que la distribución de carga corresponda al modelo TIP4P aunque la posición efectiva de la carga negativa se calcule internamente con el pair_style. [blog.csdn](https://blog.csdn.net/2301_81427007/article/details/139868234)

## Construcción del sistema

- Se define una caja cúbica de \([-5,5]\) en cada dirección y se crea un `create_box` con 2 tipos de átomos y espacio para enlaces y ángulos.  
- `molecule water tip3p.mol` lee la topología de una molécula de agua rígida.  
- `create_atoms 0 random 33 34564 NULL mol water 25367 overlap 1.33` inserta 33 moléculas de agua aleatoriamente en la caja, controlando que no se solapen demasiado (parámetro `overlap`).  

## Topología y restricciones

- `bond_style zero` y `angle_style zero` indican que los enlaces y ángulos no aportan energía; solo se usan las distancias/ángulos de equilibrio (0.9574 Å y 104.52°) como geometría objetivo. [gensoft.pasteur](https://gensoft.pasteur.fr/docs/lammps/2020.03.03/Howto_tip4p.html)
- `fix rigid all shake 0.001 10 10000 b 1 a 1` aplica SHAKE para mantener **enlaces y ángulos rígidos**, lo que es lo esperado para un modelo rígido TIP4P/TIP3P. [afs.enea](https://www.afs.enea.it/software/lammps/doc19/html/Howto_tip4p.html)

## Protocolo de simulación

- Primero se hace una **minimización** `minimize 0.0 0.0 1000 10000` para relajar la estructura inicial.  
- Después se resetea el contador de pasos y se usa un **timestep de 1.0 fs**, típico para agua rígida con SHAKE. [gensoft.pasteur](https://gensoft.pasteur.fr/docs/lammps/2020.03.03/Howto_tip4p.html)
- `velocity all create 300.0 5463576` asigna velocidades Maxwellianas a 300 K.  
- `fix integrate all nvt temp 300 300 100.0` corre un ensamble NVT a 300 K con un termostato Nose–Hoover.  

## Salida y propósito

- `thermo_style custom step temp press etotal pe` y `thermo 1000` imprimen cada 1000 pasos temperatura, presión y energías.  
- `run 20000` ejecuta 20 000 pasos (20 ps con 1 fs).  
- `write_data tip4p-implicit.data nocoeff` guarda un archivo de datos ya equilibrado, sin coeficientes de interacción, que puede usarse como estado inicial “limpio” para simulaciones posteriores con el mismo modelo TIP4P.  

Este script es un ejemplo canónico de cómo usar el **pair_style TIP4P implícito** con una geometría de TIP3P, ajustando cargas y aplicando SHAKE para obtener un pequeño “cubo” de agua equilibrada listo para reutilizarse. [blog.csdn](https://blog.csdn.net/2301_81427007/article/details/139868234)
