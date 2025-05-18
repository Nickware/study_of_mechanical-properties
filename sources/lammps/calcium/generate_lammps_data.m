function generate_lammps_data(num_atoms, output_filename)
    % Genera un archivo de datos para LAMMPS compatible con atom_style full
    % Uso: generate_lammps_data(20, 'output.data')
    
    % Parámetros del sistema
    box_size = [40.0, 10.0, 10.0];  % x,y,z
    atom_types = {'Ca', 'Si', 'O', 'H'};
    type_masses = [40.078, 28.085, 15.999, 1.008];
    type_charges = [2.0, 4.0, -0.8, 0.4];
    
    % Abrir archivo de salida
    fid = fopen(output_filename, 'w');
    
    % Escribir encabezado (formato exacto para atom_style full)
    fprintf(fid, 'LAMMPS data file for Calcium Silicate Hydrate (C-S-H)\n\n');
    fprintf(fid, '%d atoms\n', num_atoms);
    fprintf(fid, '4 atom types\n');
    fprintf(fid, '0 bonds\n');
    fprintf(fid, '0 angles\n\n');
    
    % Dimensiones de la caja
    fprintf(fid, '%.1f %.1f xlo xhi\n', -box_size(1)/2, box_size(1)/2);
    fprintf(fid, '%.1f %.1f ylo yhi\n', -box_size(2)/2, box_size(2)/2);
    fprintf(fid, '%.1f %.1f zlo zhi\n\n', -box_size(3)/2, box_size(3)/2);
    
    % Masas
    fprintf(fid, 'Masses\n\n');
    for i = 1:4
        fprintf(fid, '%d %.3f\n', i, type_masses(i));
    end
    fprintf(fid, '\n');
    
    % Átomos (formato full: ID molecule-ID type q x y z)
    fprintf(fid, 'Atoms # full\n\n');
    for i = 1:num_atoms
        % Distribución de tipos
        r = rand();
        if r < 0.2       % 20% Ca
            atom_type = 1;
            mol_id = 1;
        elseif r < 0.44  % 24% Si
            atom_type = 2;
            mol_id = 1;
        elseif r < 0.83  % 39% O
            atom_type = 3;
            mol_id = 1;
        else             % 17% H
            atom_type = 4;
            mol_id = 1;
        end
        
        pos = box_size .* rand(1,3) - box_size/2;
        
        % Formato FULL CORRECTO:
        % atom-ID molecule-ID atom-type q x y z
        fprintf(fid, '%d %d %d %.1f %.5f %.5f %.5f\n', ...
                i, mol_id, atom_type, type_charges(atom_type), pos(1), pos(2), pos(3));
    end
    
    % Velocidades (todas cero)
    fprintf(fid, '\nVelocities\n\n');
    for i = 1:num_atoms
        fprintf(fid, '%d 0.0 0.0 0.0\n', i);
    end
    
    fclose(fid);
    fprintf('Archivo %s generado correctamente para atom_style full.\n', output_filename);
end

% Función de conveniencia
function add_atoms(num_atoms)
    generate_lammps_data(num_atoms, 'output.data');
end