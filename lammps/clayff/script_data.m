% Generar 200 posiciones aleatorias para LAMMPS data file
clc; clear;

% Número de átomos
num_atoms = 200;

% Rango de posiciones en cada eje
x_min = 0; x_max = 20;
y_min = 0; y_max = 20;
z_min = 0; z_max = 20;

% Definir tipos de átomos (1: CSH, 2: PE)
atom_types = [1, 2]; % Puedes modificar según necesites

% Abrir archivo de salida
fid = fopen("generated_positions.data", "w");

% Escribir encabezado
fprintf(fid, "LAMMPS data file for generated atoms\n\n");
fprintf(fid, "%d atoms\n\n", num_atoms);
fprintf(fid, "Atoms\n\n");

% Generar posiciones y escribir en el archivo
for i = 1:num_atoms
    atom_id = i;
    atom_type = atom_types(mod(i, length(atom_types)) + 1);
    charge = 0.0; % Modificar si es necesario
    x = x_min + (x_max - x_min) * rand();
    y = y_min + (y_max - y_min) * rand();
    z = z_min + (z_max - z_min) * rand();
    
    % Escribir en formato adecuado para LAMMPS
    fprintf(fid, "%d %d %.1f %.3f %.3f %.3f\n", atom_id, atom_type, charge, x, y, z);
end

% Cerrar archivo
fclose(fid);

disp("Archivo 'generated_positions.data' creado exitosamente!");
