"""
Reorder mechanism emergence multiplier blocks to match BACTERIA_MECHANISM_EXCLUSIONS.md grouping.
"""

# Define the desired bacteria order based on BACTERIA_MECHANISM_EXCLUSIONS.md
BACTERIA_ORDER = [
    # Gram-Negative Bacteria
    # Enterobacteriaceae
    "escherichia_coli",
    "klebsiella_pneumoniae",
    "citrobacter_spp.",
    "enterobacter_spp.",
    "enterobacter_cloacae",
    "morganella_spp.",
    "proteus_spp.",
    "serratia_spp.",
    "p_stuartii",  # Providencia stuartii
    "salmonella_enterica_serovar_typhi",
    "salmonella_enterica_serovar_paratyphi_a",
    "invasive_non-typhoidal_salmonella_spp.",
    "shigella_spp.",
    "yersinia_enterocolitica",
    
    # Non-fermenting Gram-Negatives
    "pseudomonas_aeruginosa",
    "acinetobacter_baumannii",
    "stenotrophomonas_maltophilia",
    "burkholderia_cepacia_complex",
    
    # Other Gram-Negatives
    "vibrio_cholerae",
    "campylobacter_jejuni",
    "helicobacter_pylori",
    "neisseria_gonorrhoeae",
    "neisseria_meningitidis",
    "moraxella_catarrhalis",
    "haemophilus_influenzae",
    "legionella_pneumophila",
    
    # Gram-Positive Bacteria
    # Staphylococci
    "staphylococcus_aureus",
    "staphylococcus_epidermidis",
    
    # Streptococci
    "streptococcus_pneumoniae",
    "streptococcus_pyogenes",
    "streptococcus_agalactiae",
    
    # Enterococci
    "enterococcus_faecalis",
    "enterococcus_faecium",
    
    # Other Gram-Positives
    "listeria_monocytogenes",
    "clostridioides_difficile",
    "bacteroides_fragilis",
    "bordetella_pertussis",
    
    # Atypical Bacteria (No Cell Wall)
    "mycoplasma_genitalium",
    "mycoplasma_pneumoniae",
    
    # Obligate Intracellular/Special Cases
    "chlamydia_trachomatis",
    "treponema_pallidum",
    
    # Acid-Fast Bacteria
    "mdr_mycobacterium_tuberculosis",
]

# Group headers with their bacteria
GROUP_HEADERS = {
    "escherichia_coli": "// Gram-Negative Bacteria - Enterobacteriaceae",
    "pseudomonas_aeruginosa": "// Non-fermenting Gram-Negatives",
    "vibrio_cholerae": "// Other Gram-Negatives",
    "staphylococcus_aureus": "// Gram-Positive Bacteria - Staphylococci",
    "streptococcus_pneumoniae": "// Streptococci",
    "enterococcus_faecalis": "// Enterococci",
    "listeria_monocytogenes": "// Other Gram-Positives",
    "mycoplasma_genitalium": "// Atypical Bacteria (No Cell Wall)",
    "chlamydia_trachomatis": "// Obligate Intracellular/Special Cases",
    "mdr_mycobacterium_tuberculosis": "// Acid-Fast Bacteria",
}

# Bacteria-specific comments
BACTERIA_COMMENTS = {
    "escherichia_coli": "// E. coli - Lower emergence multipliers reflect high clinical prevalence of resistance",
    "klebsiella_pneumoniae": "// Klebsiella pneumoniae - The 'Plasmid Sponge' with high ESBL/Carbapenemase acquisition",
    "acinetobacter_baumannii": "// Acinetobacter baumannii - Pan-drug resistance potential",
    "pseudomonas_aeruginosa": "// Pseudomonas aeruginosa - Multi-mechanism resistance under ICU pressure",
    "neisseria_gonorrhoeae": "// Neisseria gonorrhoeae - 'Superbug' potential with rapid resistance acquisition",
    "mycoplasma_genitalium": "// Mycoplasma genitalium - High mutation rates in 23S (Macrolide) and ParC (FQ)",
    "mdr_mycobacterium_tuberculosis": "// MDR Mycobacterium tuberculosis - Chromosomal mutations drive XDR",
}

def read_config_section():
    """Read the mechanism emergence multiplier section from config.rs"""
    with open('src/config.rs', 'r') as f:
        lines = f.readlines()
    
    # Find start: look for first bacteria mechanism line (E. coli)
    start_idx = None
    for i, line in enumerate(lines):
        if 'bacteria_escherichia_coli_mechanism_' in line and 'map.insert(' in line:
            # Go back to find the comment line
            for j in range(i-1, max(0, i-10), -1):
                if '// E. coli' in lines[j]:
                    start_idx = j
                    break
            break
    
    # Find end: look for end marker
    end_idx = None
    for i, line in enumerate(lines):
        if 'end_mechanism_emergence_multiplier_parameters' in line:
            end_idx = i
            break
    
    if start_idx is None or end_idx is None:
        raise ValueError(f"Could not find boundaries: start={start_idx}, end={end_idx}")
    
    # Get content before, section, and after
    before = ''.join(lines[:start_idx])
    section = ''.join(lines[start_idx:end_idx])
    after = ''.join(lines[end_idx:])
    
    return before, section, after

def parse_blocks(section):
    """Parse all bacteria blocks from the section"""
    blocks = {}
    lines = section.split('\n')
    
    current_bacteria = None
    current_block = []
    
    for line in lines:
        # Check if this is a mechanism line
        if 'bacteria_' in line and '_mechanism_' in line and '.insert(' in line:
            # Extract bacteria name
            bacteria_part = line.split('bacteria_')[1].split('_mechanism_')[0]
            
            if current_bacteria != bacteria_part:
                # Save previous block
                if current_bacteria and current_block:
                    blocks[current_bacteria] = current_block
                
                # Start new block
                current_bacteria = bacteria_part
                current_block = [line]
            else:
                current_block.append(line)
        elif current_block and (line.strip().startswith('//') or line.strip() == ''):
            # Skip comments and blank lines between bacteria (don't include in blocks)
            continue
    
    # Save last block
    if current_bacteria and current_block:
        blocks[current_bacteria] = current_block
    
    return blocks

def reorder_blocks(blocks):
    """Reorder blocks according to BACTERIA_ORDER"""
    output_lines = []
    
    for bacteria in BACTERIA_ORDER:
        if bacteria not in blocks:
            print(f"WARNING: Bacteria '{bacteria}' not found in config.rs")
            continue
        
        # Add group header if this bacteria starts a new group
        if bacteria in GROUP_HEADERS:
            output_lines.append(f"        {GROUP_HEADERS[bacteria]}\n")
        
        # Add bacteria-specific comment if exists
        if bacteria in BACTERIA_COMMENTS:
            output_lines.append(f"        {BACTERIA_COMMENTS[bacteria]}\n")
        
        # Add all lines for this bacteria
        for line in blocks[bacteria]:
            output_lines.append(f"{line}\n")
        
        # Add blank line after each bacteria block (except last)
        if bacteria != BACTERIA_ORDER[-1]:
            output_lines.append("\n")
    
    return ''.join(output_lines)

def main():
    print("Reading config.rs...")
    before, section, after = read_config_section()
    
    print("Parsing bacteria blocks...")
    blocks = parse_blocks(section)
    print(f"Found {len(blocks)} bacteria blocks")
    
    print("Reordering blocks...")
    reordered = reorder_blocks(blocks)
    
    print("Writing reordered config.rs...")
    new_content = before + reordered + after
    
    with open('src/config.rs', 'w') as f:
        f.write(new_content)
    
    print("✓ Successfully reordered mechanism blocks")
    print(f"  Order now matches BACTERIA_MECHANISM_EXCLUSIONS.md grouping")

if __name__ == '__main__':
    main()
