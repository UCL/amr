import re
with open('src/rules/mod.rs', 'r') as f:
    c = f.read()

c = re.sub(
    r'\|\s*\"erythromycin\"\s*\|\s*\"azithromycin\"\s*\|\s*\"clarithromycin\"\s*//\s*Extrudes bulky macrolides',
    '',
    c
)

with open('src/rules/mod.rs', 'w') as f:
    f.write(c)
