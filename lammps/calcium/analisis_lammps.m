%% Script compatible con Octave para visualizar variables en ventanas separadas
clear all; close all; clc;

% Leer el archivo log.lammps
filename = 'log.lammps';
file_content = fileread(filename);

% Extraer todos los datos numéricos
data_pattern = 'Step\s+Temp\s+E_pair\s+E_mol\s+TotEng\s+Press.*?\n([\d\s\.\-eE+]+)';
data_matches = regexp(file_content, data_pattern, 'tokens');

% Combinar todos los bloques de datos
all_data = [];
for i = 1:length(data_matches)
    block = strtrim(data_matches{i}{1});
    lines = strsplit(block, '\n');
    
    for j = 1:length(lines)
        if ~isempty(lines{j})
            values = sscanf(lines{j}, '%f %f %f %f %f %f')';
            all_data = [all_data; values];
        end
    end
end

% Extraer columnas
steps = all_data(:,1);
variables = {
    {'Temp', 'Temperatura (K)', [0, max(all_data(:,2))*1.1], 'b'},
    {'E_pair', 'Energía Potencial (kcal/mol)', [min(all_data(:,3))*1.1, max(all_data(:,3))*1.1], 'r'},
    {'E_mol', 'Energía Molecular (kcal/mol)', [min(all_data(:,4))*1.1, max(all_data(:,4))*1.1], 'g'},
    {'TotEng', 'Energía Total (kcal/mol)', [min(all_data(:,5))*1.1, max(all_data(:,5))*1.1], 'm'},
    {'Press', 'Presión', [min(all_data(:,6))*1.1, max(all_data(:,6))*1.1], 'k'}
};

% Función alternativa a smooth() para Octave
movmean_octave = @(x, k) filter(ones(1, k)/k, 1, x);

% Crear una figura para cada variable
for i = 1:size(variables, 1)
    fig = figure('Name', variables{i}{1}, 'Position', [100+(i-1)*50, 100+(i-1)*50, 800, 400]);
    plot(steps, all_data(:,i+1), variables{i}{4});
    title([variables{i}{1} ' vs Paso de Simulación']);
    xlabel('Step');
    ylabel(variables{i}{2});
    ylim(variables{i}{3});
    grid on;
    set(gca, 'FontSize', 12);
    
    % Añadir línea de tendencia (versión Octave)
    hold on;
    smoothed = movmean_octave(all_data(:,i+1), 50); % Promedio móvil de 50 puntos
    plot(steps, smoothed, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    legend('Datos', 'Tendencia (promedio móvil)');
end