% analyze_water_system.m - Análisis completo de sistemas de agua TIP4P en LAMMPS
function analyze_water_system(data_file)
    % Leer el archivo de datos LAMMPS
    [atoms, bonds, angles, box] = read_lammps_data_octave(data_file);
    
    % Verificar que tenemos datos válidos
    if isempty(atoms)
        error('No se pudieron leer átomos del archivo');
    end
    
    % 1. Estadísticas básicas del sistema
    display_system_stats(atoms, bonds, angles, box);
    
    % 2. Distribución de distancias OH y HH
    plot_bond_distributions(atoms, bonds, box);
    
    % 3. Distribución de ángulos H-O-H
    plot_angle_distributions(angles);
    
    % 4. Función de distribución radial (RDF) O-O
    compute_oo_rdf(atoms, box);
    
    % 5. Visualización 3D del sistema
    plot_3d_structure(atoms, box);
end

%% Funciones auxiliares

function [atoms, bonds, angles, box] = read_lammps_data_octave(filename)
    atoms = [];
    bonds = [];
    angles = [];
    box = zeros(3,2);
    
    fid = fopen(filename, 'r');
    if fid == -1
        error('No se pudo abrir el archivo %s', filename);
    end
    
    % Leer hasta encontrar la sección de la caja
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(strfind(line, 'xlo xhi'))
            box(1,:) = sscanf(line, '%f %f')';
            box(2,:) = sscanf(fgetl(fid), '%f %f')';
            box(3,:) = sscanf(fgetl(fid), '%f %f')';
            break;
        end
    end
    
    % Buscar sección de átomos
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(strfind(line, 'Atoms'))
            while ~feof(fid)
                line = fgetl(fid);
                if isempty(line), continue; end
                if ~isempty(strfind(line, 'Velocities')) || ~isempty(strfind(line, 'Bonds'))
                    break;
                end
                atoms = [atoms; sscanf(line, '%f')'];
            end
            break;
        end
    end
    
    % Buscar sección de bonds
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(strfind(line, 'Bonds'))
            while ~feof(fid)
                line = fgetl(fid);
                if isempty(line), continue; end
                if ~isempty(strfind(line, 'Angles')) || ~isempty(strfind(line, 'Dihedrals'))
                    break;
                end
                bonds = [bonds; sscanf(line, '%f')'];
            end
            break;
        end
    end
    
    % Buscar sección de ángulos
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(strfind(line, 'Angles'))
            while ~feof(fid)
                line = fgetl(fid);
                if isempty(line), continue; end
                if ~isempty(strfind(line, 'Dihedrals')) || ~isempty(strfind(line, 'Impropers'))
                    break;
                end
                angles = [angles; sscanf(line, '%f')'];
            end
            break;
        end
    end
    
    fclose(fid);
end

function display_system_stats(atoms, bonds, angles, box)
    fprintf('\n=== Estadísticas del Sistema ===\n');
    fprintf('Número de átomos: %d\n', size(atoms,1));
    fprintf('Número de moléculas de agua: %d\n', size(atoms,1)/3);
    fprintf('Dimensiones de la caja: %.2f x %.2f x %.2f Å\n', ...
            box(1,2)-box(1,1), box(2,2)-box(2,1), box(3,2)-box(3,1));
    
    % Calcular densidad (masa en uma, volumen en Å³)
    total_mass = sum(atoms(:,4)); % Suma de todas las masas atómicas
    volume = prod(box(:,2)-box(:,1)); % Volumen en Å³
    density = (total_mass * 1.660539) / volume; % Convertir a g/cm³
    
    fprintf('Densidad aproximada: %.3f g/cm³\n', density);
    fprintf('Volumen por molécula: %.3f Å³\n', volume/(size(atoms,1)/3));
end

function plot_bond_distributions(atoms, bonds, box)
    % Calcular distancias OH con PBC
    oh_distances = [];
    box_size = box(:,2) - box(:,1);
    
    for i = 1:size(bonds,1)
        if bonds(i,2) == 1 % Solo bonds tipo OH
            atom1 = atoms(atoms(:,1) == bonds(i,3), 5:7);
            atom2 = atoms(atoms(:,1) == bonds(i,4), 5:7);
            
            % Aplicar condiciones periódicas de contorno
            r = atom1 - atom2;
            r = r - box_size' .* round(r ./ box_size');
            oh_distances = [oh_distances; norm(r)];
        end
    end
    
    % Calcular distancias HH (dentro de la misma molécula)
    hh_distances = [];
    mol_ids = atoms(:,3); % ID de molécula
    
    for mol = unique(mol_ids)'
        h_atoms = atoms(atoms(:,3) == mol & atoms(:,2) == 2, 5:7);
        if size(h_atoms,1) == 2
            r = h_atoms(1,:) - h_atoms(2,:);
            r = r - box_size' .* round(r ./ box_size');
            hh_distances = [hh_distances; norm(r)];
        end
    end
    
    figure;
    subplot(2,1,1);
    hist(oh_distances, 15);
    title('Distribución de distancias O-H');
    xlabel('Distancia (Å)');
    ylabel('Frecuencia');
    grid on;
    
    subplot(2,1,2);
    hist(hh_distances, 15);
    title('Distribución de distancias H-H');
    xlabel('Distancia (Å)');
    ylabel('Frecuencia');
    grid on;
end

function plot_angle_distributions(angles)
    if isempty(angles)
        disp('No se encontraron datos de ángulos');
        return;
    end
    
    angle_values = angles(:,5);
    
    figure;
    hist(angle_values, 15);
    title('Distribución de ángulos H-O-H');
    xlabel('Ángulo (grados)');
    ylabel('Frecuencia');
    grid on;
    
    fprintf('\n=== Estadísticas de ángulos H-O-H ===\n');
    fprintf('Ángulo promedio: %.2f°\n', mean(angle_values));
    fprintf('Desviación estándar: %.2f°\n', std(angle_values));
    fprintf('Ángulo mínimo: %.2f°\n', min(angle_values));
    fprintf('Ángulo máximo: %.2f°\n', max(angle_values));
end

function compute_oo_rdf(atoms, box)
    oxygens = atoms(atoms(:,2) == 1, 5:7);
    if isempty(oxygens)
        disp('No se encontraron átomos de oxígeno');
        return;
    end
    
    num_oxygens = size(oxygens,1);
    box_size = box(:,2) - box(:,1);
    
    max_dist = min(box_size)/2;
    dr = 0.1; % bin size
    edges = 0:dr:max_dist;
    counts = zeros(length(edges)-1,1);
    
    for i = 1:num_oxygens
        for j = i+1:num_oxygens
            r = oxygens(i,:) - oxygens(j,:);
            % Aplicar condiciones periódicas de contorno
            r = r - box_size' .* round(r ./ box_size');
            dist = norm(r);
            
            if dist <= max_dist
                bin = floor(dist/dr) + 1;
                counts(bin) = counts(bin) + 2; % Contar ambos pares
            end
        end
    end
    
    % Cálculo del factor de normalización
    r = edges(1:end-1) + dr/2;
    norm_factor = (4/3)*pi*(edges(2:end).^3 - edges(1:end-1).^3);
    norm_factor = norm_factor * (num_oxygens^2 / prod(box_size));
    g_r = counts ./ norm_factor';
    
    figure;
    plot(r, g_r, 'LineWidth', 2);
    title('Función de Distribución Radial (RDF) O-O');
    xlabel('Distancia r (Å)');
    ylabel('g(r)');
    grid on;
end

function plot_3d_structure(atoms, box)
    oxygens = atoms(atoms(:,2) == 1, :);
    hydrogens = atoms(atoms(:,2) == 2, :);
    
    figure;
    hold on;
    
    % Dibujar caja
    vertices = [box(1,1) box(2,1) box(3,1); box(1,2) box(2,1) box(3,1);
                box(1,2) box(2,2) box(3,1); box(1,1) box(2,2) box(3,1);
                box(1,1) box(2,1) box(3,2); box(1,2) box(2,1) box(3,2);
                box(1,2) box(2,2) box(3,2); box(1,1) box(2,2) box(3,2)];
    
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
    
    patch('Vertices', vertices, 'Faces', faces, ...
          'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1);
    
    % Dibujar átomos
    scatter3(oxygens(:,5), oxygens(:,6), oxygens(:,7), 100, 'r', 'filled');
    scatter3(hydrogens(:,5), hydrogens(:,6), hydrogens(:,7), 50, 'b', 'filled');
    
    % Dibujar bonds (opcional)
    for i = 1:size(oxygens,1)
        mol_id = oxygens(i,3);
        h_atoms = atoms(atoms(:,3) == mol_id & atoms(:,2) == 2, 5:7);
        
        for j = 1:size(h_atoms,1)
            plot3([oxygens(i,5) h_atoms(j,1)], ...
                  [oxygens(i,6) h_atoms(j,2)], ...
                  [oxygens(i,7) h_atoms(j,3)], 'b-', 'LineWidth', 1);
        end
    end
    
    % Configuración de la vista
    axis equal;
    xlabel('X (Å)');
    ylabel('Y (Å)');
    zlabel('Z (Å)');
    title('Estructura 3D del Sistema de Agua TIP4P');
    legend({'Caja', 'Oxígeno', 'Hidrógeno', 'Enlaces'});
    view(3);
    rotate3d on;
    hold off;
end