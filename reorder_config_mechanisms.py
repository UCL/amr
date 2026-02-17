"""
Reorder bacteria mechanism blocks in config.rs to match biological classification.
"""

# Biological order from BACTERIA_MECHANISM_EXCLUSIONS.md
BIOLOGICAL_ORDER = [
    # Gram-Negative - Enterobacteriaceae
    ("escherichia_coli", "E. coli - Lower emergence multipliers reflect high clinical prevalence of resistance"),
    ("klebsiella_pneumoniae", "Klebsiella pneumoniae - The 'Plasmid Sponge' with high ESBL/Carbapenemase acquisition"),
    ("citrobacter_spp.", None),
    ("enterobacter_spp.", None),
    ("enterobacter_cloacae", None),
    ("morganella_spp.", None),
    ("proteus_spp.", None),
    ("serratia_spp.", None),
    ("p_stuartii", None),
    ("salmonella_enterica_serovar_typhi", None),
    ("salmonella_enterica_serovar_paratyphi_a", None),
    ("invasive_non-typhoidal_salmonella_spp.", None),
    ("shigella_spp.", None),
    ("yersinia_enterocolitica", None),
    
    # Gram-Negative - Non-fermenting
    ("pseudomonas_aeruginosa", "Pseudomonas aeruginosa - Multi-mechanism resistance under ICU pressure"),
    ("acinetobacter_baumannii", "Acinetobacter baumannii - Pan-drug resistance potential"),
    ("stenotrophomonas_maltophilia", None),
    ("burkholderia_cepacia_complex", None),
    
    # Gram-Negative - Other
    ("vibrio_cholerae", None),
    ("campylobacter_jejuni", None),
    ("helicobacter_pylori", None),
    ("neisseria_gonorrhoeae", "Neisseria gonorrhoeae - 'Superbug' potential with rapid resistance acquisition"),
    ("neisseria_meningitidis", None),
    ("moraxella_catarrhalis", None),
    ("haemophilus_influenzae", None),
    ("legionella_pneumophila", None),
    
    # Gram-Positive - Staphylococci
    ("staphylococcus_aureus", None),
    ("staphylococcus_epidermidis", None),
    
    # Gram-Positive - Streptococci
    ("streptococcus_pneumoniae", None),
    ("streptococcus_pyogenes", None),
    ("streptococcus_agalactiae", None),
    
    # Gram-Positive - Enterococci
    ("enterococcus_faecalis", None),
    ("enterococcus_faecium", None),
    
    # Gram-Positive - Other
    ("listeria_monocytogenes", None),
    ("clostridioides_difficile", None),
    ("bacteroides_fragilis", None),
    ("bordetella_pertussis", None),
    
    # Atypical (No Cell Wall)
    ("mycoplasma_genitalium", "Mycoplasma genitalium - High mutation rates in 23S (Macrolide) and ParC (FQ)"),
    ("mycoplasma_pneumoniae", None),
    
    # Obligate Intracellular/Special Cases
    ("chlamydia_trachomatis", None),
    ("treponema_pallidum", None),
    
    # Acid-Fast
    ("mdr_mycobacterium_tuberculosis", "MDR Mycobacterium tuberculosis - Chromosomal mutations drive XDR"),
]

GROUP_HEADERS = {
    "escherichia_coli": "\n        // Gram-Negative Bacteria - Enterobacteriaceae",
    "pseudomonas_aeruginosa": "\n        // Non-fermenting Gram-Negatives",
    "vibrio_cholerae": "\n        // Other Gram-Negatives",
    "staphylococcus_aureus": "\n        // Gram-Positive Bacteria - Staphylococci",
    "streptococcus_pneumoniae": "\n        // Streptococci",
    "enterococcus_faecalis": "\n        // Enterococci",
    "listeria_monocytogenes": "\n        // Other Gram-Positives",
    "mycoplasma_genitalium": "\n        // Atypical Bacteria (No Cell Wall)",
    "chlamydia_trachomatis": "\n        // Obligate Intracellular/Special Cases",
    "mdr_mycobacterium_tuberculosis": "\n        // Acid-Fast Bacteria",
}

def read_config_file():
    """Read config.rs preserving line endings."""
    with open('src/config.rs', 'r', newline='') as f:
        return f.readlines()

def find_section_bounds(lines):
    """Find start and end of mechanism section."""
    start_idx = None
    end_idx = None
    
    # Find first bacteria mechanism line (E. coli)
    for i, line in enumerate(lines):
        if 'bacteria_escherichia_coli_mechanism_' in line and 'map.insert(' in line:
            # Look back for comment
            for j in range(i-1, max(0, i-5), -1):
                if '// E. coli' in lines[j] or lines[j].strip() == '':
                    continue
                start_idx = j + 1  # Start after the last non-comment/non-blank
                break
            if start_idx is None:
                start_idx = i
            break
    
    # Find last mechanism line (Burkholderia is currently last)
    for i in range(len(lines)-1, -1, -1):
        if '_mechanism_modification_mcr_1_emergence_multiplier' in lines[i]:
            # End is right after this line
            end_idx = i + 1
            # Skip any trailing blank lines
            while end_idx < len(lines) and lines[end_idx].strip() == '':
                end_idx += 1
            break
    
    return start_idx, end_idx

def extract_bacteria_blocks(lines, start_idx, end_idx):
    """Extract all bacteria blocks as dictionary."""
    blocks = {}
    current_bacteria = None
    current_block = []
    current_comment = None
    
    for i in range(start_idx, end_idx):
        line = lines[i]
        
        # Check for bacteria-specific comment
        if line.strip().startswith('//') and 'mechanism' not in line.lower():
            # This might be a bacteria comment
            current_comment = line
            continue
        
        # Check if this is a mechanism line
        if 'bacteria_' in line and '_mechanism_' in line and 'map.insert(' in line:
            # Extract bacteria name
            bacteria_part = line.split('bacteria_')[1].split('_mechanism_')[0]
            
            if current_bacteria != bacteria_part:
                # Save previous block
                if current_bacteria and current_block:
                    blocks[current_bacteria] = {
                        'comment': blocks.get(current_bacteria, {}).get('comment'),
                        'lines': current_block
                    }
                
                # Start new block
                current_bacteria = bacteria_part
                current_block = [line]
                if current_comment:
                    blocks[current_bacteria] = {'comment': current_comment, 'lines': []}
                    current_comment = None
            else:
                current_block.append(line)
    
    # Save last block
    if current_bacteria and current_block:
        if current_bacteria in blocks:
            blocks[current_bacteria]['lines'] = current_block
        else:
            blocks[current_bacteria] = {'comment': None, 'lines': current_block}
    
    return blocks

def build_reordered_section(blocks):
    """Build reordered section with group headers."""
    output = []
    
    for bacteria_slug, custom_comment in BIOLOGICAL_ORDER:
        if bacteria_slug not in blocks:
            print(f"WARNING: Bacteria '{bacteria_slug}' not found in config.rs!")
            continue
        
        block_data = blocks[bacteria_slug]
        
        # Add group header if this starts a new group
        if bacteria_slug in GROUP_HEADERS:
            output.append(GROUP_HEADERS[bacteria_slug] + '\n')
        
        # Add bacteria-specific comment (use custom or extracted)
        if custom_comment:
            output.append(f"        // {custom_comment}\n")
        elif block_data['comment']:
            output.append(block_data['comment'])
        
        # Add all mechanism lines
        output.extend(block_data['lines'])
        
        # Add blank line after each bacteria
        output.append('\n')
    
    return output

def write_config_file(lines, start_idx, end_idx, new_section):
    """Write updated config.rs."""
    new_content = lines[:start_idx] + new_section + lines[end_idx:]
    
    with open('src/config.rs', 'w', newline='') as f:
        f.writelines(new_content)

def main():
    print("Step 1: Reading config.rs...")
    lines = read_config_file()
    print(f"  Total lines: {len(lines)}")
    
    print("Step 2: Finding section bounds...")
    start_idx, end_idx = find_section_bounds(lines)
    if start_idx is None or end_idx is None:
        print(f"ERROR: Could not find section bounds (start={start_idx}, end={end_idx})")
        return
    print(f"  Section: lines {start_idx+1} to {end_idx} ({end_idx-start_idx} lines)")
    
    print("Step 3: Extracting bacteria blocks...")
    blocks = extract_bacteria_blocks(lines, start_idx, end_idx)
    print(f"  Found {len(blocks)} bacteria")
    
    # Validate
    for bacteria_slug, _ in BIOLOGICAL_ORDER:
        if bacteria_slug in blocks:
            num_lines = len(blocks[bacteria_slug]['lines'])
            if num_lines != 25:
                print(f"  WARNING: {bacteria_slug} has {num_lines} lines (expected 25)")
    
    print("Step 4: Building reordered section...")
    new_section = build_reordered_section(blocks)
    print(f"  Generated {len(new_section)} lines")
    
    print("Step 5: Writing updated config.rs...")
    write_config_file(lines, start_idx, end_idx, new_section)
    
    print("\n✓ Successfully reordered mechanism blocks!")
    print("  Next: Run 'cargo check' to verify compilation")

if __name__ == '__main__':
    main()
