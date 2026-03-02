import re
import ast

def main():
    # 1. READ RUST FILE
    with open('src/simulation/population.rs', 'r', encoding='utf-8') as f:
        rs_content = f.read()
        
    rs_match = re.search(r'pub const DRUG_SHORT_NAMES:\s*&\[&str\]\s*=\s*&\[(.*?)\];', rs_content, re.DOTALL)
    if not rs_match:
        print("Could not find DRUG_SHORT_NAMES in rust file.")
        return
        
    rs_drugs = set(re.findall(r'"([^"]+)"', rs_match.group(1)))
    
    # 2. READ PYTHON FILE
    with open('amr_simulation_output_analysis/calibration_summary.py', 'r', encoding='utf-8') as f:
        py_content = f.read()
        
    # Find the CROSS_RESISTANCE_CLASS_OVERRIDES tuple definition
    py_match = re.search(r'CROSS_RESISTANCE_CLASS_OVERRIDES.*?=\s*(\(.*?\))\n\s*@dataclass', py_content, re.DOTALL)
    if not py_match:
        print("Could not find CROSS_RESISTANCE_CLASS_OVERRIDES in python file.")
        return
        
    try:
        # Evaluate the tuple structure safely
        overrides = ast.literal_eval(py_match.group(1))
        
        py_drugs = set()
        for class_label, drugs in overrides:
            for d in drugs:
                py_drugs.add(d)
                
    except Exception as e:
        print(f"Error parsing python tuple: {e}")
        return
        
    print(f"Total drugs in population.rs: {len(rs_drugs)}")
    print(f"Total drugs in calibration_summary.py: {len(py_drugs)}")
    
    # 3. COMPARE
    in_rs_not_py = rs_drugs - py_drugs
    in_py_not_rs = py_drugs - rs_drugs
    
    if not in_rs_not_py and not in_py_not_rs:
        print("\nLists are IDENTICAL.")
    else:
        print("\nIn population.rs but NOT in calibration_summary.py:")
        for d in sorted(in_rs_not_py):
            print(f"  - {d}")
            
        print("\nIn calibration_summary.py but NOT in population.rs:")
        for d in sorted(in_py_not_rs):
            print(f"  - {d}")

if __name__ == "__main__":
    main()
