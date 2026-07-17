# MATLAB/Octave para analizar datos de Lammps

Este script es una herramienta en **MATLAB/Octave** para analizar datos de una simulación de LAMMPS exportados a un archivo CSV, calculando estadísticas básicas y generando varias gráficas de diagnóstico a partir de la evolución temporal de magnitudes termodinámicas. [mathworks](https://www.mathworks.com/help/matlab/ref/csvread.html)

## Flujo general del script

1. Lee un archivo `CSV` con datos de la simulación (por ejemplo, salida de `thermo` exportada).  
2. Extrae columnas: paso, temperatura, presión, energía total y energía potencial.  
3. Calcula estadísticas (media y desviación estándar) de cada magnitud.  
4. Muestra esas estadísticas por pantalla.  
5. Genera gráficas de evolución temporal, correlación temperatura–presión y un histograma de temperaturas.  
6. Guarda las estadísticas en un archivo `.mat` y las gráficas en `.png`.  

Todo esto está encapsulado en la función principal `analyze_lammps(csv_file)` y funciones auxiliares internas.

## Lectura de datos

```matlab
data = csvread(csv_file, 1, 0);  % Skip header row

steps       = data(:,1);
temp        = data(:,2);
press       = data(:,3);
tot_energy  = data(:,4);
pot_energy  = data(:,5);
```

- Usa `csvread` para leer el archivo, saltando la primera fila (cabecera). En MATLAB moderno se recomienda `readmatrix` en lugar de `csvread`, pero en Octave `csvread` sigue siendo habitual. [mathworks](https://www.mathworks.com/help/matlab/ref/csvread.html?nocookie=true)
- Asume un formato de columnas:  
  1. Paso de simulación  
  2. Temperatura  
  3. Presión  
  4. Energía total  
  5. Energía potencial  

## Estadísticas que calcula

```matlab
stats.mean_temp       = mean(temp);
stats.std_temp        = std(temp);
stats.mean_press      = mean(press)/1000;  % kPa
stats.std_press       = std(press)/1000;
stats.mean_tot_energy = mean(tot_energy);
stats.std_tot_energy  = std(tot_energy);
stats.mean_pot_energy = mean(pot_energy);
stats.std_pot_energy  = std(pot_energy);
```

- Calcula medias y desviaciones estándar.  
- Convierte la presión de la unidad original (probablemente atm o bar) a kPa dividiendo por 1000, asumiendo que los datos están en unidades compatibles.  
- Imprime valores tipo:  
  - Temperatura media ± desviación estándar (K).  
  - Presión media ± desviación estándar (kPa).  
  - Energía total y potencial medias ± desviación estándar (kcal/mol).  

Estas métricas sirven para evaluar si la simulación está equilibrada (fluctuaciones razonables, valores cercanos a los esperados).

## Gráficas que genera

La función `create_plots` produce varias figuras:

1. **Temperatura vs pasos**  
2. **Presión vs pasos** (en kPa)  
3. **Energía total vs pasos**  
4. **Energía potencial vs pasos**  

Todas en una figura 2×2, con rejilla y líneas marcadas.

Además:

- Una gráfica adicional de **dispersión temperatura–presión**, útil para ver correlaciones entre ambas magnitudes.  
- Un **histograma de temperaturas** (5 bins), para ver la distribución de T alrededor de la media.  

Guarda:

- `simulation_plots.png`: figura con las curvas de evolución.  
- `temp_press_correlation.png`: figura de correlación T–P.

## Guardado de resultados

```matlab
save_stats(stats, 'simulation_stats.mat');
```

- La función `save_stats` guarda la estructura `stats` en un archivo `simulation_stats.mat` usando `save -struct`, es decir, cada campo pasa a ser una variable independiente en el archivo.  
- Esto permite cargar rápidamente los resultados estadísticos en MATLAB/Octave para análisis posterior.

## Uso desde Octave por línea de comandos

Al final del archivo:

```matlab
if exist('argv', 'var') && length(argv) > 0
    analyze_lammps(argv{1});
else
    disp('Usage en Octave: analyze_lammps(''thermo_data.csv'')');
end
```

- En Octave, `argv()` devuelve los argumentos de la línea de comandos cuando se ejecuta el script como programa. [octave.sourceforge](https://octave.sourceforge.io/octave/function/argv.html)
- Esto permite correr algo como:  
  ```bash
  octave analyze_lammps.m thermo_data.csv
  ```  
  y que el script procese directamente el archivo `thermo_data.csv`.  

## En qué contexto usar este script

- Está pensado para **post-procesar salidas de LAMMPS** (típicamente un CSV generado a partir de `thermo_style custom` exportado).  
- Es útil para:
  - Ver si la temperatura y la energía se han estabilizado.  
  - Analizar fluctuaciones de presión.  
  - Documentar estadísticas de producción en un archivo `.mat` fácil de reutilizar.  
