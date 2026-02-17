with open('src/config.rs', 'rb') as f:
    content = f.read()

marker = b'start_mechanism_emergence_multiplier_parameters'
idx = content.find(marker)
print(f'Marker found at: {idx}')
if idx > 0:
    print(f'Context: {repr(content[idx-50:idx+100])}')
