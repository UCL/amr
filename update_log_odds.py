import re

with open('MODEL_DESCRIPTION.md', 'r', encoding='utf-8') as f:
    text = f.read()

# Let's find tables that have variables containing 'log_odds'
tables = re.finditer(r'\| Parameter \|\s*Baseline Value\s*\|.*?\n(?:\|[-\s]+\|[-\s]+\|.*?\n)(?:\|.*?\|.*?\n)+', text)

new_text = text
for t in tables:
    table_text = t.group(0)
    # Check if the table contains 'log_odds'
    if 'log_odds' in table_text and ('base_log_odds' in table_text or 'log_odds_ratio' in table_text or 'log_odds' in table_text):
        new_table = table_text.replace('| Parameter | Baseline Value | Description |', '| Parameter | Baseline (Log-odds ratio) | Description |')
        new_table = new_table.replace('| Parameter | Baseline Value |', '| Parameter | Baseline (Log-odds ratio) |')
        new_text = new_text.replace(table_text, new_table)

with open('MODEL_DESCRIPTION.md', 'w', encoding='utf-8') as f:
    f.write(new_text)

print("Updated tables with log odds!")
