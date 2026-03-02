import re
import ast

def parse_rs():
    with open('src/simulation/population.rs', 'r', encoding='utf-8') as f:
        content = f.read()
    match = re.search(r'pub const DRUG_SHORT_NAMES: &\[&str\] = &\[(.*?)\];', content, re.DOTALL)
    if not match:
        print("Could not find DRUG_SHORT_NAMES in rs")
        return set()
    drugs = re.findall(r'"([^"]+)"', match.group(1))
    return set(drugs)

def parse_py():
    with open('amr_simulation_output_analysis/calibration_summary.py', 'r', encoding='utf-8') as f:
        text = f.read()
    
    start = text.find('CROSS_RESISTANCE_CLASS_OVERRIDES')
    end = text.find('@dataclass', start)
    if end == -1:
        end = text.find('class ', start)
    block = text[start:end]
    
    drugs = set()
    matches = re.findall(r'\(\s*"[^"]+"\s*,\s*\(([^)]+)\)', block)
    for m in matches:
        d_matches = re.findall(r'"([^"]+)"', m)
        drugs.update(d_matches)
    return drugs

def main():
    rs_drugs = parse_rs()
    py_drugs = parse_py()

    print("--- Drugs comparison ---")
    
    rs_only = sorted(rs_drugs - py_drugs)
    if rs_only:
        print("\nFound in population.rs but NOT in calibration_summary.py:")
        for d in rs_only:
            print(f"  - {d}")
    else:
        print("\nNo drugs are exclusive to population.rs")

    py_only = sorted(py_drugs - rs_drugs)
    if py_only:
        print("\nFound in calibration_summary.py but NOT in population.rs:")
        for d in py_only:
            print(f"  - {d}")
    else:
        print("\nNo drugs are exclusive to calibration_summary.py")

if __name__ == '__main__':
    main()
