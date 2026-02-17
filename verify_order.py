with open('src/config.rs', 'r') as f:
    lines = f.readlines()

# Find all bacteria in order of appearance
bacteria_order = []
for i, line in enumerate(lines):
    if 'bacteria_' in line and '_mechanism_enzyme_esbl_ctx_m' in line:
        bacteria_name = line.split('bacteria_')[1].split('_mechanism_')[0]
        bacteria_order.append((i+1, bacteria_name))

print("Order of bacteria in config.rs:\n")
for i, (line_num, bacteria) in enumerate(bacteria_order, 1):
    print(f"{i:2}. Line {line_num:5}: {bacteria}")
