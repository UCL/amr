import re

with open('MODEL_DESCRIPTION.md', 'r', encoding='utf-8') as f:
    content = f.read()

# Replace different variations of Default in the header of the tables

# | Variable | Default | Description | -> | Parameter | Baseline Value | Description |
content = re.sub(r'\|\s*Variable\s*\|\s*Default\s*\|\s*Description\s*\|', r'| Parameter | Baseline Value | Description |', content)

# | Variable pattern | Default | Description | -> | Parameter pattern | Baseline Value | Description |
content = re.sub(r'\|\s*Variable pattern\s*\|\s*Default\s*\|', r'| Parameter pattern | Baseline Value |', content)

# | Variable | Default | -> | Parameter | Baseline Value |
content = re.sub(r'\|\s*Variable\s*\|\s*Default\s*\|', r'| Parameter | Baseline Value |', content)

with open('MODEL_DESCRIPTION.md', 'w', encoding='utf-8') as f:
    f.write(content)

print("Markdown updated!")

