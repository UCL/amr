#!/usr/bin/env python3
"""
Audit resistance mechanisms - check if all drugs showing resistance have defined mechanisms
"""

import re
from collections import defaultdict

def parse_calibration_file(filepath):
    """Parse calibration summary to find drugs with resistance > 0%"""
    drugs_with_resistance = defaultdict(set)  # drug -> set of bacteria
    
    with open(filepath, 'r', encoding='utf-8') as f:
        in_resistance_table = False
        
        for line in f:
            # Look for the resistance benchmarks table
            if 'Resistance Benchmarks (percent resistant)' in line:
                in_resistance_table = True
                next(f)  # Skip header line
                continue
            
            if not in_resistance_table:
                continue
            
            # Skip lines that are too short or have '---' (empty values)
            if len(line.strip()) < 50:
                continue
                
            # Skip header-like lines
            if 'Bacteria' in line and 'Drug' in line:
                continue
            
            # Parse data lines - extract bacteria, drug, and first numeric value
            # The format is: Bacteria  Drug  DrugClass  Number1  Number2  ...
            try:
                # Split by multiple spaces to separate columns
                parts = re.split(r'\s{2,}', line.strip())
                if len(parts) < 4:
                    continue
                
                bacteria = parts[0].strip()
                drug = parts[1].strip()
                drug_class = parts[2].strip()
                
                # Skip if first field doesn't look like a bacteria name
                if not bacteria or bacteria.replace('_', '').replace(' ', '').replace('.', '').isalpha() == False:
                    continue
                
                # Find all decimal numbers in the line
                numbers = re.findall(r'\b\d+\.\d+\b', line)
                if len(numbers) >= 1:
                    # First number should be infection resistance simulation %
                    resistance_sim_str = numbers[0]
                    resistance_sim = float(resistance_sim_str)
                    
                    # Check if resistance > 0
                    if resistance_sim > 0.0:
                        drugs_with_resistance[drug].add(bacteria)
                        
            except (ValueError, IndexError, AttributeError) as e:
                continue
    
    return drugs_with_resistance


def main():
    calibration_file = r'c:\Users\w3sth\rust_amr_project\output_graphs\calibration_summary_574817.txt'
    
    print("Parsing calibration file...")
    drugs_with_resistance = parse_calibration_file(calibration_file)
    
    print(f"\n=== DRUGS WITH SIMULATED RESISTANCE > 0% ===\n")
    print(f"Total unique drugs with resistance: {len(drugs_with_resistance)}\n")
    
    for drug in sorted(drugs_with_resistance.keys()):
        bacteria_list = drugs_with_resistance[drug]
        print(f"{drug:30s} ({len(bacteria_list):2d} bacteria)")
    
    # Save to file for further analysis
    with open('drugs_with_resistance.txt', 'w') as f:
        for drug in sorted(drugs_with_resistance.keys()):
            bacteria_list = sorted(drugs_with_resistance[drug])
            f.write(f"\n{drug}:\n")
            for b in bacteria_list:
                f.write(f"  - {b}\n")
    
    print(f"\n\nDetailed list saved to drugs_with_resistance.txt")
    print(f"\nDrugs to check: {', '.join(sorted(drugs_with_resistance.keys()))}")

if __name__ == '__main__':
    main()
