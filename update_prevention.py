import re

with open('MODEL_DESCRIPTION.md', 'r', encoding='utf-8') as f:
    content = f.read()


new_section = """### 6.8 Antibiotic infection prevention

When an individual is actively taking an antibiotic, the treatment acts as prophylaxis against new incoming infections. The efficacy of this prevention represents the proportional reduction in infection acquisition probability for incoming sensitive strains.

| Parameter | Baseline Value | Description |
|-----------|----------------|-------------|
| ntibiotic_infection_prevention_efficacy | 0.7 | The relative reduction in new infection establishment probability when the individual is already taking an effective antibiotic. |"""

pattern = r'### 6\.8 Antibiotic infection prevention\s*\| Parameter \| Baseline Value \|\s*\|[-\s]+\|[-\s]+\|\s*\| ntibiotic_infection_prevention_efficacy \| 0\.7 \|'
if re.search(pattern, content):
    content = re.sub(pattern, new_section, content)
    with open('MODEL_DESCRIPTION.md', 'w', encoding='utf-8') as f:
        f.write(content)
    print("Markdown section 6.8 updated via regex successfully!")
else:
    print("Still couldn't find the section to replace.")

