% save this as analyze_lammps.m
function analyze_lammps(csv_file)
    % Cargar datos desde el archivo CSV
    data = csvread(csv_file, 1, 0);  % Skip header row
    
    steps = data(:,1);
    temp = data(:,2);
    press = data(:,3);
    tot_energy = data(:,4);
    pot_energy = data(:,5);
    
    % Análisis estadístico básico
    stats = struct();
    stats.mean_temp = mean(temp);
    stats.std_temp = std(temp);
    stats.mean_press = mean(press)/1000;  % Convert to kPa
    stats.std_press = std(press)/1000;
    stats.mean_tot_energy = mean(tot_energy);
    stats.std_tot_energy = std(tot_energy);
    stats.mean_pot_energy = mean(pot_energy);
    stats.std_pot_energy = std(pot_energy);
    
    % Mostrar estadísticas
    disp('Estadísticas de la simulación:');
    disp(['Temperatura media: ', num2str(stats.mean_temp), ' ± ', num2str(stats.std_temp), ' K']);
    disp(['Presión media: ', num2str(stats.mean_press), ' ± ', num2str(stats.std_press), ' kPa']);
    disp(['Energía total media: ', num2str(stats.mean_tot_energy), ' ± ', num2str(stats.std_tot_energy), ' kcal/mol']);
    disp(['Energía potencial media: ', num2str(stats.mean_pot_energy), ' ± ', num2str(stats.std_pot_energy), ' kcal/mol']);
    
    % Crear gráficas
    create_plots(steps, temp, press, tot_energy, pot_energy);
    
    % Guardar estadísticas
    save_stats(stats, 'simulation_stats.mat');
end

function create_plots(steps, temp, press, tot_energy, pot_energy)
    % Configurar tamaño de figura
    figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
    
    % Gráfica de temperatura
    subplot(2,2,1);
    plot(steps, temp, 'b-o', 'LineWidth', 1.5);
    xlabel('Pasos de simulación');
    ylabel('Temperatura (K)');
    title('Evolución de la temperatura');
    grid on;
    
    % Gráfica de presión
    subplot(2,2,2);
    plot(steps, press/1000, 'r-o', 'LineWidth', 1.5);
    xlabel('Pasos de simulación');
    ylabel('Presión (kPa)');
    title('Evolución de la presión');
    grid on;
    
    % Gráfica de energía total
    subplot(2,2,3);
    plot(steps, tot_energy, 'g-o', 'LineWidth', 1.5);
    xlabel('Pasos de simulación');
    ylabel('Energía total (kcal/mol)');
    title('Evolución de la energía total');
    grid on;
    
    % Gráfica de energía potencial
    subplot(2,2,4);
    plot(steps, pot_energy, 'm-o', 'LineWidth', 1.5);
    xlabel('Pasos de simulación');
    ylabel('Energía potencial (kcal/mol)');
    title('Evolución de la energía potencial');
    grid on;
    
    % Gráfica adicional: correlación temperatura-presión
    figure;
    scatter(temp, press/1000, 100, 'filled');
    xlabel('Temperatura (K)');
    ylabel('Presión (kPa)');
    title('Correlación temperatura-presión');
    %lsline;
    grid on;
    print('temp_press_correlation.png', '-dpng');
    
    % Guardar todas las gráficas
    print('simulation_plots.png', '-dpng');
end

function save_stats(stats, filename)
    save(filename, '-struct', 'stats');
    disp(['Estadísticas guardadas en: ', filename]);
end

% Ejecutar análisis si se llama directamente
if exist('argv', 'var') && length(argv) > 0
    analyze_lammps(argv{1});
else
    disp('Usage en Octave: analyze_lammps(''thermo_data.csv'')');
end
