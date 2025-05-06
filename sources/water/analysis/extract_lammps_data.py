# save this as extract_lammps_data.py
import re
import sys
import pandas as pd

def extract_lammps_data(output_file):
    with open(output_file, 'r') as f:
        content = f.read()

    # Pattern to match thermodynamic data
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
    csv_file = output_file.replace('.out', '.csv').replace('.log', '.csv')
    df.to_csv(csv_file, index=False)
    print(f"Extrayendo datos a {csv_file}")
    return csv_file

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python extract_lammps_data.py <lammps_output_file>")
        sys.exit(1)
    
    output_file = sys.argv[1]
    extract_lammps_data(output_file)
