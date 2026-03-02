import re

def main():
    with open('src/simulation/population.rs', 'r', encoding='utf-8') as f:
        rs_content = f.read()
        
    rs_match = re.search(r'pub const DRUG_SHORT_NAMES:\s*&\[&str\]\s*=\s*&\[(.*?)\];', rs_content, re.DOTALL)
    if not rs_match:
        print("Could not find DRUG_SHORT_NAMES in rust file.")
        return
        
    rs_drugs = set(re.findall(r'"(.*?)"', rs_match.group(1)))
    
    with open('amr_simulation_output_analysis/calibration_summary.py', 'r', encoding='utf-8') as f:
        py_content = f.read()
        
    py_match = re.search(r'CROSS_RESISTANCE_CLASS_OVERRIDES.*?=\s*\((.*?)\)\n\s*@dataclass', py_content, re.DOTALL)
    if not py_match:
        print("Could not find CROSS_RESISTANCE_CLASS_OVERRIDES in python file.")
        return
        
    # Python file has tuples in tuples: ("Class name", ("drug1", "drug2"))
    # The class names start with capital letters or have spaces, while drugs are lowercase with underscores.
    # Actually, we can just find all strings and filter out those that are capitalized or have spaces/parentheses.
    # Let's extract all strings in the matched section:
    py_strings = re.findall(r'"(.*?)"', py_match.group(1))
    
    # Drugs are all lowercase and underscores
    py_drugs = set(s for s in py_strings if re.match(r'^[a-z_]+$', s))
    
    print(f"Drugs in population.rs: {len(rs_drugs)}")
    
    with open('amr_simulation_output_analysis/calibration_summary.py', 'r', encoding='utf-8') as f:
        py_content = f.read()
        
    py_match = re.search(r'CROSS_RESISTANCE_CLASS_OVERRIDES.*?=\s*\((.*?)\)\n\s*@dataclass', py_content, re.DOTALL)
    if not py_match:
        print("Could not find CROSS_RESISTANCE_CLASS_OVERRIDES in python file.")
        return
        
    py_strings = set(re.findall(r'"([^"]+)"', py_match.group(1)))
    
    # Exclude known headers/labels (class names)
    py_drugs = set(s for s in py_strings if s.islower() and " " not in s and "(" not in s)
    
    print(f"Drugs in calibration_summary.py: {len(py_drugs)}")
    
    in_rs_not_py = rs_drugs - py_drugs
    in_py_not_rs = py_drugs - rs_drugs
    
    print("\nIn population.rs but NOT in calibration_summary.py:")
    for d in sorted(in_rs_not_py):
        print(f"  - {d}")
        
    print("\nIn calibration_summary.py but NOT in population.rs:")
    for d in sorted(in_py_not_rs):
        print(f"  - {d}")

if __name__ == "__main__":
    main()
