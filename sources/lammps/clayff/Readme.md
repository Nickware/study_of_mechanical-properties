# LAMMPS con potencial CLAYFF

------

# 📘 Simulación de C-S-H con aditivos de PE en LAMMPS

## 🧩 1. Script base en LAMMPS con potencial CLAYFF

```lammps
units real
atom_style full
boundary p p p

read_data clay_structure.data

pair_style lj/cut 10.0
pair_coeff * * 0.1 3.5

bond_style harmonic
angle_style harmonic

fix 1 all nvt temp 300.0 300.0 100.0

thermo 100
dump 1 all xyz 100 output.xyz

run 10000
```

**Descripción:**

- Define una simulación con unidades tipo `real` y estilo de átomos `full`.
- Usa el potencial Lennard-Jones (`lj/cut`) como aproximación inicial.
- Aplica condiciones periódicas en las tres direcciones.
- Ejecuta un ensamble NVT a 300 K durante 10,000 pasos.

------

## 📂 2. Archivo `clay_structure.data` (estructura genérica)

```lammps
LAMMPS data file for clay simulation

100 atoms
100 bonds
100 angles

Atoms
1 1 0.0 5.0 5.0 0.0
2 1 0.0 5.5 5.5 0.0
3 2 0.0 6.0 6.0 0.0
...
100 2 0.0 10.0 10.0 0.0

Bonds
1 1 1 2
...
100 1 99 100

Angles
1 1 1 2 3
...
100 1 98 99 100
```

**Descripción:**

- Contiene 100 átomos con tipos alternados (`1` y `2`).
- Define enlaces y ángulos básicos para simular interacciones moleculares.
- Las posiciones están distribuidas en el plano XY con Z constante.

------

## 🧪 3. Script para simular C-S-H con aditivos de PE

```lammps
units real
atom_style full
boundary p p p

read_data csh_pe_structure.data

pair_style lj/cut 10.0
pair_coeff * * 0.1 3.5

bond_style harmonic
angle_style harmonic

fix 1 all nvt temp 300.0 300.0 100.0

thermo 100
dump 1 all xyz 100 output.xyz

run 10000
```

**Descripción:**

- Similar al script base, pero orientado a simular una mezcla de C-S-H y PE.
- El archivo `csh_pe_structure.data` contiene la estructura combinada.
- Ideal para estudiar interacciones entre matriz cementicia y polímeros.

------

## 🧱 4. Archivo `csh_pe_structure.data` (modelo mixto)

```lammps
LAMMPS data file for C-S-H with PE additives

200 atoms
180 bonds
160 angles

Atoms
1 1 0.0 5.0 5.0 0.0
2 1 0.0 5.5 5.5 0.0
3 2 0.0 6.0 6.0 0.0
...
200 2 0.0 15.0 15.0 0.0

Bonds
1 1 1 2
...
180 2 179 180

Angles
1 1 1 2 3
...
160 2 159 160
```

**Descripción:**

- Contiene 200 átomos: tipo `1` para C-S-H y tipo `2` para PE.
- Define enlaces y ángulos para simular la estructura molecular.
- Las posiciones están distribuidas en el espacio tridimensional.

------

## 🧮 5. Script en Octave para generar posiciones

```octave
clc; clear;

num_atoms = 200;
x_min = 0; x_max = 20;
y_min = 0; y_max = 20;
z_min = 0; z_max = 20;
atom_types = [1, 2];

fid = fopen("generated_positions.data", "w");

fprintf(fid, "LAMMPS data file for generated atoms\n\n");
fprintf(fid, "%d atoms\n\n", num_atoms);
fprintf(fid, "Atoms\n\n");

for i = 1:num_atoms
    atom_id = i;
    atom_type = atom_types(mod(i, length(atom_types)) + 1);
    charge = 0.0;
    x = x_min + (x_max - x_min) * rand();
    y = y_min + (y_max - y_min) * rand();
    z = z_min + (z_max - z_min) * rand();
    fprintf(fid, "%d %d %.1f %.3f %.3f %.3f\n", atom_id, atom_type, charge, x, y, z);
end

fclose(fid);
disp("Archivo 'generated_positions.data' creado exitosamente!");
```

**Descripción:**

- Genera 200 posiciones aleatorias en un cubo de 20x20x20 Å.
- Alterna tipos de átomos entre C-S-H (`1`) y PE (`2`).
- Guarda el resultado en un archivo `.data` listo para usar en LAMMPS.
