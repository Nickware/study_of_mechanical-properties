function generate_lammps_data(num_atoms, output_filename)
    % Versión definitiva para simulaciones MPI en LAMMPS
    % Uso: generate_lammps_data(200, 'calcium_mpi.data')
    
    % 1. Configuración optimizada para MPI
    box_size = [40.0, 15.0, 15.0]; % Caja más grande para mejor distribución
    num_regions = 4; % Igual al número de procesadores que usarás
    
    % Propiedades atómicas
    atom_props = {
        1, 40.078, 2.0, 1;   % Ca
        2, 28.085, 4.0, 1;    % Si
        3, 15.999, -0.8, 1;   % O
        4, 1.008, 0.4, 2      % H (mol_id = 2 para enlaces)
    };
    
    % 2. Generación de átomos con distribución por regiones
    atoms = cell(num_atoms, 5);
    region_H = cell(num_regions, 1); % Hidrógenos por región
    
    % Definir regiones espaciales
    region_edges = linspace(-box_size(1)/2, box_size(1)/2, num_regions+1);
    
    for i = 1:num_atoms
        % Asignar región basada en posición x (para dominio simple en MPI)
        region = mod(i-1, num_regions) + 1;
        x_min = region_edges(region);
        x_max = region_edges(region+1);
        
        % Posición dentro de la región asignada
        pos = [x_min + (x_max-x_min)*rand(), ...
               box_size(2)*(rand()-0.5), ...
               box_size(3)*(rand()-0.5)];
        
        % Asignar tipo atómico
        if rand() < 0.25 % 25% H
            type = 4;
            region_H{region} = [region_H{region}; i];
        else
            r = rand();
            if r < 0.2
                type = 1; % Ca
            elseif r < 0.45
                type = 2; % Si
            else
                type = 3; % O
            end
        end
        
        atoms{i,1} = i;
        atoms{i,2} = atom_props{type,4}; % mol_id
        atoms{i,3} = type;
        atoms{i,4} = atom_props{type,3}; % charge
        atoms{i,5} = pos;
    end
    
    % 3. Generación de enlaces SOLO dentro de cada región
    bonds = {};
    bond_count = 0;
    max_dist = 1.5; % Distancia máxima para enlace
    
    for region = 1:num_regions
        H_in_region = region_H{region};
        for i = 1:length(H_in_region)
            for j = i+1:length(H_in_region)
                dist = norm(atoms{H_in_region(i),5} - atoms{H_in_region(j),5});
                if dist < max_dist
                    bond_count = bond_count + 1;
                    bonds{bond_count} = [bond_count, 1, H_in_region(i), H_in_region(j)];
                end
            end
        end
    end
    
    % 4. Validación CRÍTICA de enlaces
    valid_bonds = {};
    for i = 1:length(bonds)
        id1 = bonds{i}(3);
        id2 = bonds{i}(4);
        % Verificar que ambos átomos existen y son de tipo H
        if id1 <= num_atoms && id2 <= num_atoms && ...
           atoms{id1,3} == 4 && atoms{id2,3} == 4
            valid_bonds{end+1} = bonds{i};
        end
    end
    
    % 5. Generación de ángulos solo si hay suficientes enlaces
    angles = {};
    if length(valid_bonds) >= 2
        for i = 1:length(valid_bonds)-1
            shared = intersect(valid_bonds{i}(3:4), valid_bonds{i+1}(3:4));
            if ~isempty(shared)
                angle_atoms = unique([valid_bonds{i}(3:4), valid_bonds{i+1}(3:4)]);
                if length(angle_atoms) == 3
                    angles{end+1} = [length(angles)+1, 1, angle_atoms(1), shared, angle_atoms(3)];
                end
            end
        end
    end
    
    % 6. Escritura del archivo de datos
    fid = fopen(output_filename, 'w');
    
    % Encabezado
    fprintf(fid, 'LAMMPS data file - MPI Optimized\n\n');
    fprintf(fid, '%d atoms\n', num_atoms);
    fprintf(fid, '4 atom types\n');
    fprintf(fid, '%d bonds\n', length(valid_bonds));
    fprintf(fid, '1 bond types\n');
    fprintf(fid, '%d angles\n', length(angles));
    fprintf(fid, '1 angle types\n\n');
    
    % Box
    fprintf(fid, '%.1f %.1f xlo xhi\n', -box_size(1)/2, box_size(1)/2);
    fprintf(fid, '%.1f %.1f ylo yhi\n', -box_size(2)/2, box_size(2)/2);
    fprintf(fid, '%.1f %.1f zlo zhi\n\n', -box_size(3)/2, box_size(3)/2);
    
    % Masas
    fprintf(fid, 'Masses\n\n');
    for i = 1:4
        fprintf(fid, '%d %.3f\n', i, atom_props{i,2});
    end
    
    % Coeficientes
    fprintf(fid, '\nBond Coeffs\n\n1 500.0 1.0\n');
    fprintf(fid, '\nAngle Coeffs\n\n1 100.0 109.47\n');
    
    % Atoms
    fprintf(fid, '\nAtoms # full\n\n');
    for i = 1:num_atoms
        fprintf(fid, '%d %d %d %.1f %.5f %.5f %.5f\n', ...
                atoms{i,1}, atoms{i,2}, atoms{i,3}, atoms{i,4}, atoms{i,5}(1), atoms{i,5}(2), atoms{i,5}(3));
    end
    
    % Bonds
    fprintf(fid, '\nBonds\n\n');
    for i = 1:length(valid_bonds)
        fprintf(fid, '%d %d %d %d\n', valid_bonds{i}(1), valid_bonds{i}(2), valid_bonds{i}(3), valid_bonds{i}(4));
    end
    
    % Angles
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
    
    % Reporte final
    fprintf('Archivo %s generado exitosamente:\n', output_filename);
    fprintf('- %d átomos totales\n', num_atoms);
    fprintf('- %d enlaces válidos (solo intra-región)\n', length(valid_bonds));
    fprintf('- %d ángulos\n', length(angles));
end