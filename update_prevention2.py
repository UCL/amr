import re

with open('MODEL_DESCRIPTION.md', 'r', encoding='utf-8') as f:
    lines = f.readlines()

new_lines = []
skip = False
for line in lines:
    if "### 6.8 Antibiotic infection prevention" in line:
        new_lines.append("### 6.8 Antibiotic infection prevention\n\n")
        new_lines.append("When an individual is actively taking an antibiotic, the treatment acts as prophylaxis against new incoming infections. The efficacy of this prevention represents the proportional reduction in infection acquisition probability for incoming sensitive strains.\n\n")
        new_lines.append("| Parameter | Baseline Value | Description |\n")
        new_lines.append("|-----------|----------------|-------------|\n")
        new_lines.append("| ntibiotic_infection_prevention_efficacy | 0.7 | The relative reduction in new infection establishment probability when the individual is already taking an effective antibiotic. |\n")
        skip = True
    elif skip and "## 7. Resistance Dynamics" in line:
        skip = False
        new_lines.append("\n---\n\n\n")
        new_lines.append(line)
        
    elif not skip:
        new_lines.append(line)

with open('MODEL_DESCRIPTION.md', 'w', encoding='utf-8') as f:
    f.writelines(new_lines)
    
print("Updated successfully!")

