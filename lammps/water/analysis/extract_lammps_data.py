# extract_lammps_data_fixed.py
import re
import sys
import pandas as pd
from pathlib import Path

def extract_lammps_data(output_file):
    with open(output_file, 'r') as f:
        content = f.read()

    # Patrón para que coincida con los datos termodinámicos
    pattern = r'^\s*(\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
    matches = re.finditer(pattern, content, re.MULTILINE)
    
    data = []
    for match in matches:
        data.append({
            'Step': int(match.group(1)),
            'Temp': float(match.group(2)),
            'Press': float(match.group(3)),
            'TotEng': float(match.group(4)),
            'PotEng': float(match.group(5))
        })
    
    df = pd.DataFrame(data)
    
    # Crear nombre de archivo de salida garantizando extensión .csv
    input_path = Path(output_file)
    csv_file = input_path.with_suffix('.csv')
    
    df.to_csv(csv_file, index=False)
    print(f"Data extracted to {csv_file}")
    return csv_file

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python extract_lammps_data.py <lammps_output_file>")
        sys.exit(1)
    
    output_file = sys.argv[1]
    extract_lammps_data(output_file)