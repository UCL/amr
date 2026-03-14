import re

with open('MODEL_DESCRIPTION.md', 'r', encoding='utf-8') as f:
    content = f.read()

# Pattern to find the drug introduction dates table
# And drop the second column (Time step)
# Keep only the first and third columns

def replace_table(match):
    lines = match.group(0).strip().split('\n')
    new_lines = []
    
    # Process header
    header = lines[0].split('|')
    new_header = [header[0], header[1].strip(), header[3].strip(), header[4]]
    new_lines.append('|'.join(new_header))
    
    # Process separator
    sep = lines[1].split('|')
    new_sep = [sep[0], sep[1].strip(), sep[3].strip(), sep[4]]
    new_lines.append('|'.join(new_sep))
    
    # Process rows
    for line in lines[2:]:
        if not line.strip():
            continue
        cols = line.split('|')
        new_cols = [cols[0], cols[1].strip(), cols[3].strip(), cols[4]]
        new_lines.append('|'.join(new_cols))
        
    return '\n'.join(new_lines)

new_content = re.sub(
    r'\| Drug\s*\|\s*Time step\s*\|\s*~Year\s*\|.*?(?=\n\n|\Z)',
    replace_table,
    content,
    flags=re.DOTALL
)

# Update the preceding text to match the new structure
new_content = new_content.replace(
    'Each antibiotic becomes available at a specific time step:',
    'Each antibiotic becomes available in a specific year:'
)

with open('MODEL_DESCRIPTION.md', 'w', encoding='utf-8') as f:
    f.write(new_content)

print("Markdown updated!")

