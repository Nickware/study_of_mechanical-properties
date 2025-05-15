import numpy as np
from numpy.linalg import norm
import os

def load_last_frame(xyz_file):
    """Carga el último frame de un archivo XYZ de LAMMPS"""
    with open(xyz_file, 'r') as f:
        lines = f.readlines()
    
    # Encontrar el último frame
    frame_starts = [i for i, line in enumerate(lines) if line.strip().isdigit()]
    if not frame_starts:
        raise ValueError("No se encontraron frames en el archivo XYZ")
    
    last_frame_start = frame_starts[-1]
    num_atoms = int(lines[last_frame_start].strip())
    comment = lines[last_frame_start + 1].strip()
    atom_lines = lines[last_frame_start + 2: last_frame_start + 2 + num_atoms]
    
    atoms = []
    positions = []
    for line in atom_lines:
        parts = line.split()
        atoms.append(parts[0])
        positions.append([float(x) for x in parts[1:4]])
    
    return atoms, np.array(positions), comment

def calculate_angles(positions):
    """Calcula todos los ángulos H-C-H en una molécula de metano"""
    # Asumimos que el carbono es el primer átomo
    C = positions[0]
    H_positions = positions[1:]
    
    # Vectores CH
    vectors = [H - C for H in H_positions]
    
    # Calcular todos los ángulos únicos
    angles = []
    for i in range(len(vectors)):
        for j in range(i + 1, len(vectors)):
            v1 = vectors[i]
            v2 = vectors[j]
            cos_angle = np.dot(v1, v2) / (norm(v1) * norm(v2))
            angle = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
            angles.append(angle)
    
    return angles

def calculate_distances(positions):
    """Calcula distancias C-H"""
    C = positions[0]
    distances = [norm(H - C) for H in positions[1:]]
    return distances

def analyze_methane(xyz_file):
    """Función principal que analiza la geometría del metano"""
    try:
        # Cargar datos
        atoms, positions, comment = load_last_frame(xyz_file)
        
        # Verificar que sea metano (CH4)
        if len(positions) != 5 or atoms.count('C') != 1 or atoms.count('H') != 4:
            raise ValueError("El archivo no contiene una molécula de metano (CH4)")
        
        # Calcular propiedades geométricas
        angles = calculate_angles(positions)
        distances = calculate_distances(positions)
        
        # Resultados
        print("\n=== Análisis de Geometría del Metano ===")
        print(f"Archivo analizado: {os.path.basename(xyz_file)}")
        print(f"Frame analizado: {comment}")
        
        print("\nDistancias C-H (Å):")
        for i, d in enumerate(distances, 1):
            print(f"  C-H{i}: {d:.4f}")
        
        print("\nÁngulos H-C-H (grados):")
        angle_pairs = [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]
        for (i,j), angle in zip(angle_pairs, angles):
            print(f"  H{i}-C-H{j}: {angle:.2f}°")
        
        print("\nResumen estadístico:")
        print(f"  Distancia C-H promedio: {np.mean(distances):.4f} ± {np.std(distances):.4f} Å")
        print(f"  Ángulo H-C-H promedio: {np.mean(angles):.2f} ± {np.std(angles):.2f}°")
        print(f"  Desviación del tetraedro perfecto (109.47°): {np.mean(np.abs(np.array(angles) - 109.47)):.2f}°")
        
        return {
            'distances': distances,
            'angles': angles,
            'positions': positions
        }
    
    except Exception as e:
        print(f"Error: {str(e)}")
        return None

# Ejemplo de uso
if __name__ == "__main__":
    # Cambia esto por la ruta a tu archivo de trayectoria
    xyz_file = "traj_methane.xyz"
    
    if os.path.exists(xyz_file):
        results = analyze_methane(xyz_file)
        
        # Visualización adicional (requiere matplotlib)
        try:
            import matplotlib.pyplot as plt
            
            plt.figure(figsize=(10, 5))
            
            # Histograma de ángulos
            plt.subplot(1, 2, 1)
            plt.hist(results['angles'], bins=10, color='skyblue', edgecolor='black')
            plt.axvline(109.47, color='red', linestyle='--', label='Ángulo tetraédrico ideal')
            plt.xlabel('Ángulo H-C-H (grados)')
            plt.ylabel('Frecuencia')
            plt.title('Distribución de Ángulos')
            plt.legend()
            
            # Histograma de distancias
            plt.subplot(1, 2, 2)
            plt.hist(results['distances'], bins=10, color='lightgreen', edgecolor='black')
            plt.axvline(1.09, color='red', linestyle='--', label='Distancia C-H ideal')
            plt.xlabel('Distancia C-H (Å)')
            plt.ylabel('Frecuencia')
            plt.title('Distribución de Distancias')
            plt.legend()
            
            plt.tight_layout()
            plt.show()
            
        except ImportError:
            print("\nPara visualización gráfica, instala matplotlib: pip install matplotlib")
    else:
        print(f"Error: El archivo {xyz_file} no existe en el directorio actual")
