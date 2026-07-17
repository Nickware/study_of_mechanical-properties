% filepath: /tmp/test/extract_info.m
% Leer el archivo log.lammps
filename = 'log.lammps';
fid = fopen(filename, 'r');
if fid == -1
    error('No se pudo abrir el archivo log.lammps');
end

% Variables para almacenar datos del proceso de minimización
min_steps = [];
min_temp = [];
min_poteng = [];
min_toteng = [];
min_press = [];

% Variables para almacenar datos del proceso de dinámica molecular
dynamics_steps = [];
dynamics_temp = [];
dynamics_poteng = [];
dynamics_toteng = [];
dynamics_press = [];

% Leer línea por línea
is_minimization = true; % Bandera para diferenciar entre minimización y dinámica molecular
while ~feof(fid)
    line = fgetl(fid);

    % Detectar el inicio de la tabla de datos
    if ~isempty(strfind(line, 'Step          Temp          PotEng         TotEng         Press'))
        fgetl(fid); % Saltar línea
        while true
            data_line = fgetl(fid);
            if isempty(data_line) || ~ischar(data_line) || ~isempty(strfind(data_line, 'Loop time'))
                break;
            end
            data = sscanf(data_line, '%f');
            if length(data) == 5 % Validar que la línea contiene exactamente 5 valores
                if is_minimization
                    % Guardar datos de minimización
                    min_steps(end+1) = data(1);
                    min_temp(end+1) = data(2);
                    min_poteng(end+1) = data(3);
                    min_toteng(end+1) = data(4);
                    min_press(end+1) = data(5);
                else
                    % Guardar datos de dinámica molecular
                    dynamics_steps(end+1) = data(1);
                    dynamics_temp(end+1) = data(2);
                    dynamics_poteng(end+1) = data(3);
                    dynamics_toteng(end+1) = data(4);
                    dynamics_press(end+1) = data(5);
                end
            end
        end
        % Cambiar la bandera después de procesar la primera tabla
        is_minimization = false;
    end
end

fclose(fid);

% Graficar el proceso de minimización de energía
figure;
subplot(3, 1, 1);
plot(min_steps, min_poteng, '-o', 'LineWidth', 1.5);
xlabel('Paso');
ylabel('Energía Potencial (kcal/mol)');
title('Minimización de Energía: Energía Potencial');
grid on;

subplot(3, 1, 2);
plot(min_steps, min_toteng, '-o', 'LineWidth', 1.5);
xlabel('Paso');
ylabel('Energía Total (kcal/mol)');
title('Minimización de Energía: Energía Total');
grid on;

subplot(3, 1, 3);
plot(min_steps, min_press, '-o', 'LineWidth', 1.5);
xlabel('Paso');
ylabel('Presión (atm)');
title('Minimización de Energía: Presión');
grid on;

% Graficar el proceso de dinámica molecular
figure;
subplot(3, 1, 1);
plot(dynamics_steps, dynamics_temp, '-o', 'LineWidth', 1.5);
xlabel('Paso');
ylabel('Temperatura (K)');
title('Dinámica Molecular: Temperatura');
grid on;

subplot(3, 1, 2);
plot(dynamics_steps, dynamics_poteng, '-o', 'LineWidth', 1.5);
xlabel('Paso');
ylabel('Energía Potencial (kcal/mol)');
title('Dinámica Molecular: Energía Potencial');
grid on;

subplot(3, 1, 3);
plot(dynamics_steps, dynamics_press, '-o', 'LineWidth', 1.5);
xlabel('Paso');
ylabel('Presión (atm)');
title('Dinámica Molecular: Presión');
grid on;
