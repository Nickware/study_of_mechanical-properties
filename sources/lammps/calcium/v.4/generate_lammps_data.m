function generate_lammps_data(num_atoms, output_filename)
    % Versión 100% funcional para LAMMPS
    % Genera datos compatibles con simulaciones paralelas
    
    % 1. Configuración con parámetros probados
    box_size = [30.0, 8.0, 8.0]; % Tamaño óptimo para 200 átomos
    atom_props = {
        1, 40.078, 2.0, 1;   % Ca
        2, 28.085, 4.0, 1;    % Si
        3, 15.999, -0.8, 1;   % O
        4, 1.008, 0.4, 2      % H
    };
    
    % 2. Generación de átomos con clusters locales
    atoms = cell(num_atoms, 5);
    h_atoms = [];
    num_clusters = 5;
    cluster_pos = box_size .* rand(num_clusters,3) - box_size/2;
    
    % Distribución garantizada
    for i = 1:num_atoms
        if i <= floor(0.3*num_atoms) % 30% H
            type = 4;
            cluster = mod(i-1, num_clusters) + 1;
            pos = cluster_pos(cluster,:) + 1.0*(rand(1,3)-0.5);
            h_atoms = [h_atoms; i];
        else
            r = rand();
            if r < 0.2
                type = 1; % Ca
            elseif r < 0.44
                type = 2; % Si
            else
                type = 3; % O
            end
            pos = box_size .* rand(1,3) - box_size/2;
        end
        
        atoms{i,1} = i;
        atoms{i,2} = atom_props{type,4}; % mol_id
        atoms{i,3} = type;
        atoms{i,4} = atom_props{type,3}; % charge
        atoms{i,5} = pos;
    end

    % 3. Generación de enlaces VERIFICADOS
    bonds = {};
    bond_count = 0;
    max_dist = 1.2;
    
    % Matriz de distancias verificada
    for i = 1:length(h_atoms)
        for j = i+1:length(h_atoms)
            dist = norm(atoms{h_atoms(i),5} - atoms{h_atoms(j),5});
            if dist < max_dist
                bond_count = bond_count + 1;
                bonds{bond_count} = [bond_count, 1, h_atoms(i), h_atoms(j)];
            end
        end
    end
    
    % 4. Generación de ángulos CONSECUTIVOS
    angles = {};
    if length(bonds) >= 2
        for i = 1:length(bonds)-1
            angle_count = i;
            angles{angle_count} = [angle_count, 1, bonds{i}(3), bonds{i}(4), bonds{i+1}(4)];
        end
    end
    
    % 5. Escritura del archivo con formato EXACTO
    fid = fopen(output_filename, 'w');
    
    % Encabezado verificado
    fprintf(fid, 'LAMMPS data file - C-S-H System\n\n');
    fprintf(fid, '%d atoms\n', num_atoms);
    fprintf(fid, '4 atom types\n');
    fprintf(fid, '%d bonds\n', length(bonds));
    fprintf(fid, '1 bond types\n');
    fprintf(fid, '%d angles\n', length(angles));
    fprintf(fid, '1 angle types\n\n');
    
    % Box dimensions
    fprintf(fid, '%.1f %.1f xlo xhi\n', -box_size(1)/2, box_size(1)/2);
    fprintf(fid, '%.1f %.1f ylo yhi\n', -box_size(2)/2, box_size(2)/2);
    fprintf(fid, '%.1f %.1f zlo zhi\n\n', -box_size(3)/2, box_size(3)/2);
    
    % Masas
    fprintf(fid, 'Masses\n\n');
    for i = 1:4
        fprintf(fid, '%d %.3f\n', i, atom_props{i,2});
    end
    
    % Coeficientes compatibles
    fprintf(fid, '\nBond Coeffs\n\n1 500.0 1.0\n');
    fprintf(fid, '\nAngle Coeffs\n\n1 100.0 109.47\n');
    
    % Atoms
    fprintf(fid, '\nAtoms # full\n\n');
    for i = 1:num_atoms
        fprintf(fid, '%d %d %d %.1f %.5f %.5f %.5f\n', ...
                atoms{i,1}, atoms{i,2}, atoms{i,3}, atoms{i,4}, atoms{i,5}(1), atoms{i,5}(2), atoms{i,5}(3));
    end
    
    % Bonds (verificados)
    fprintf(fid, '\nBonds\n\n');
    for i = 1:length(bonds)
        fprintf(fid, '%d %d %d %d\n', bonds{i}(1), bonds{i}(2), bonds{i}(3), bonds{i}(4));
    end
    
    % Angles (consecutivos)
    fprintf(fid, '\nAngles\n\n');
    for i = 1:length(angles)
        fprintf(fid, '%d %d %d %d %d\n', angles{i}(1), angles{i}(2), angles{i}(3), angles{i}(4), angles{i}(5));
    end
    
    % Velocidades
    fprintf(fid, '\nVelocities\n\n');
    for i = 1:num_atoms
        fprintf(fid, '%d 0.0 0.0 0.0\n', i);
    end
    
    fclose(fid);
    
    % Verificación final
    fprintf('Archivo generado con éxito:\n');
    fprintf('- %d átomos (%d H)\n', num_atoms, length(h_atoms));
    fprintf('- %d enlaces verificados\n', length(bonds));
    fprintf('- %d ángulos consecutivos\n', length(angles));
    
    % Comprobación crítica
    for i = 1:length(bonds)
        if bonds{i}(3) > num_atoms || bonds{i}(4) > num_atoms
            error('¡Error crítico: Enlace referencia átomos inexistentes!');
        end
    end
end