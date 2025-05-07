% Script Octave para analizar datos de simulación LAMMPS

clear all
close all

% Datos extraídos manualmente del archivo de salida (cada 1000 pasos)
steps = [0:1000:20000]';
temp = [300.00, 343.38, 247.18, 373.22, 311.16, 350.24, 259.69, 285.35, 230.82, 332.09, ...
        315.71, 331.05, 275.33, 282.57, 302.53, 314.47, 270.64, 283.56, 303.27, 245.56, 326.32];
press = [121475.86, -10421.05, -5475.79, -21486.49, -13357.03, -29107.38, -19729.65, ...
         -13047.31, -19795.28, -10259.26, -21736.63, -14740.69, -30370.45, -32254.08, ...
         -12904.56, -18290.65, -11790.66, -7299.60, -32052.83, -11687.59, -2838.80];
tot_energy = [-581.08, -355.88, -143.86, -796.61, -519.39, -1237.41, -796.27, -535.18, ...
              -673.75, -323.69, -741.63, -504.16, -1053.61, -1198.29, -352.61, -613.21, ...
              -427.64, -181.11, -1203.80, -371.39, -47.61];
pot_energy = [-639.20, -422.41, -191.75, -868.92, -579.68, -1305.27, -846.58, -590.47, ...
              -718.47, -388.04, -802.80, -568.30, -1106.96, -1253.04, -411.23, -674.14, ...
              -480.08, -236.05, -1262.56, -418.97, -110.84];

% Análisis de datos
mean_temp = mean(temp)
std_temp = std(temp)
mean_press = mean(press)/1000  % en kPa
std_press = std(press)/1000
mean_energy = mean(tot_energy)
std_energy = std(tot_energy)

% Gráficas
figure;
subplot(2,2,1);
plot(steps, temp, 'b-o');
xlabel('Pasos de simulación');
ylabel('Temperatura (K)');
title('Evolución de la temperatura');
grid on;

subplot(2,2,2);
plot(steps, press/1000, 'r-o');
xlabel('Pasos de simulación');
ylabel('Presión (kPa)');
title('Evolución de la presión');
grid on;

subplot(2,2,3);
plot(steps, tot_energy, 'g-o');
xlabel('Pasos de simulación');
ylabel('Energía total (kcal/mol)');
title('Evolución de la energía total');
grid on;

subplot(2,2,4);
plot(steps, pot_energy, 'm-o');
xlabel('Pasos de simulación');
ylabel('Energía potencial (kcal/mol)');
title('Evolución de la energía potencial');
grid on;

% Correlación entre temperatura y presión
figure;
scatter(temp, press/1000, 'filled');
xlabel('Temperatura (K)');
ylabel('Presión (kPa)');
title('Correlación temperatura-presión');
% lsline;
grid on;

% Histograma de temperaturas
figure;
hist(temp, 5);
xlabel('Temperatura (K)');
ylabel('Frecuencia');
title('Distribución de temperaturas');
grid on;
