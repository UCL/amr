"""
Comment out biologically inapplicable mechanism entries in config.rs
based on BACTERIA_MECHANISM_EXCLUSIONS.md
"""

import re

# Mapping from BACTERIA_MECHANISM_EXCLUSIONS.md names to config.rs slugs
MECHANISM_NAME_MAP = {
    'EnzymeEsblCtxM': 'enzyme_esbl_ctx_m',
    'EnzymeEsblTem': 'enzyme_esbl_tem',
    'EnzymeEsblShv': 'enzyme_esbl_shv',
    'EnzymeAmpcCmy': 'enzyme_ampc_cmy',
    'EnzymeAmpcDha': 'enzyme_ampc_dha',
    'EnzymeKpc': 'enzyme_kpc',
    'EnzymeNdmVim': 'enzyme_ndm_vim',
    'EnzymeOxa48': 'enzyme_oxa_48',
    'EnzymeCat': 'enzyme_cat',
    'Enzyme16sRrmt': 'enzyme_16s_rrmt',
    'TargetSitePbp2aMecA': 'target_site_pbp2a_meca',
    'TargetSiteVanA': 'target_site_van_a',
    'TargetSiteVanB': 'target_site_van_b',
    'TargetSiteErmB': 'target_site_erm_b',
    'TargetSiteCfr': 'target_site_cfr',
    'MutationGyrAPrimary': 'mutation_gyra_primary',
    'MutationGyrAParCSecondary': 'mutation_gyra_parc_secondary',
    'ProtectionQnr': 'protection_qnr',
    'EffluxAcrabTolc': 'efflux_acrab_tolc',
    'EffluxMexxyOprm': 'efflux_mexxy_oprm',
    'GlobalEffluxPump': 'global_efflux_pump',
    'PorinLossOmpk35_36': 'porin_loss_ompk35_36',
    'PorinLossOprd': 'porin_loss_oprd',
    'GlobalPorinLoss': 'global_porin_loss',
    'ModificationMcr1': 'modification_mcr_1',
}

# Bacteria name mapping from markdown to config.rs slugs
BACTERIA_NAME_MAP = {
    'Escherichia coli': 'escherichia_coli',
    'Klebsiella pneumoniae': 'klebsiella_pneumoniae',
    'Citrobacter spp.': 'citrobacter_spp.',
    'Enterobacter spp.': 'enterobacter_spp.',
    'Enterobacter cloacae': 'enterobacter_cloacae',
    'Morganella spp.': 'morganella_spp.',
    'Proteus spp.': 'proteus_spp.',
    'Serratia spp.': 'serratia_spp.',
    'Providencia stuartii': 'p_stuartii',
    'Salmonella enterica serovar Typhi': 'salmonella_enterica_serovar_typhi',
    'Salmonella enterica serovar Paratyphi A': 'salmonella_enterica_serovar_paratyphi_a',
    'Invasive non-typhoidal Salmonella spp.': 'invasive_non-typhoidal_salmonella_spp.',
    'Shigella spp.': 'shigella_spp.',
    'Yersinia enterocolitica': 'yersinia_enterocolitica',
    'Pseudomonas aeruginosa': 'pseudomonas_aeruginosa',
    'Acinetobacter baumannii': 'acinetobacter_baumannii',
    'Stenotrophomonas maltophilia': 'stenotrophomonas_maltophilia',
    'Burkholderia cepacia complex': 'burkholderia_cepacia_complex',
    'Vibrio cholerae': 'vibrio_cholerae',
    'Campylobacter jejuni': 'campylobacter_jejuni',
    'Helicobacter pylori': 'helicobacter_pylori',
    'Neisseria gonorrhoeae': 'neisseria_gonorrhoeae',
    'Neisseria meningitidis': 'neisseria_meningitidis',
    'Moraxella catarrhalis': 'moraxella_catarrhalis',
    'Haemophilus influenzae': 'haemophilus_influenzae',
    'Legionella pneumophila': 'legionella_pneumophila',
    'Staphylococcus aureus': 'staphylococcus_aureus',
    'Staphylococcus epidermidis': 'staphylococcus_epidermidis',
    'Streptococcus pneumoniae': 'streptococcus_pneumoniae',
    'Streptococcus pyogenes': 'streptococcus_pyogenes',
    'Streptococcus agalactiae': 'streptococcus_agalactiae',
    'Enterococcus faecalis': 'enterococcus_faecalis',
    'Enterococcus faecium': 'enterococcus_faecium',
    'Listeria monocytogenes': 'listeria_monocytogenes',
    'Clostridioides difficile': 'clostridioides_difficile',
    'Bacteroides fragilis': 'bacteroides_fragilis',
    'Bordetella pertussis': 'bordetella_pertussis',
    'Mycoplasma genitalium': 'mycoplasma_genitalium',
    'Mycoplasma pneumoniae': 'mycoplasma_pneumoniae',
    'Chlamydia trachomatis': 'chlamydia_trachomatis',
    'Treponema pallidum': 'treponema_pallidum',
    'MDR Mycobacterium tuberculosis': 'mdr_mycobacterium_tuberculosis',
}

def parse_exclusions():
    """Parse BACTERIA_MECHANISM_EXCLUSIONS.md to build exclusion dictionary."""
    exclusions = {}
    
    with open('BACTERIA_MECHANISM_EXCLUSIONS.md', 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Find bacteria sections with **Name**
    bacteria_pattern = r'\*\*([^*]+)\*\*\n((?:- `[^`]+`[^\n]*\n)*)'
    
    for match in re.finditer(bacteria_pattern, content):
        bacteria_name = match.group(1).strip()
        exclusion_block = match.group(2)
        
        if bacteria_name not in BACTERIA_NAME_MAP:
            continue
            
        bacteria_slug = BACTERIA_NAME_MAP[bacteria_name]
        
        # Check for "No exclusions"
        if 'No exclusions' in exclusion_block or not exclusion_block.strip():
            exclusions[bacteria_slug] = []
            continue
        
        # Extract mechanism names from - `MechanismName` lines
        mechanism_pattern = r'- `([^`]+)`'
        excluded_mechanisms = []
        
        for mech_match in re.finditer(mechanism_pattern, exclusion_block):
            mech_name = mech_match.group(1)
            if mech_name in MECHANISM_NAME_MAP:
                excluded_mechanisms.append(MECHANISM_NAME_MAP[mech_name])
        
        exclusions[bacteria_slug] = excluded_mechanisms
    
    return exclusions

def process_config_file(exclusions):
    """Process config.rs and comment out excluded mechanisms."""
    with open('src/config.rs', 'r', newline='') as f:
        lines = f.readlines()
    
    modified_lines = []
    commented_count = 0
    
    for line in lines:
        # Check if this is a mechanism insert line
        if 'bacteria_' in line and '_mechanism_' in line and 'map.insert(' in line:
            # Extract bacteria and mechanism
            match = re.search(r'bacteria_([^_]+(?:_[^_]+)*)_mechanism_([^_]+(?:_[^_]+)*)_emergence', line)
            
            if match:
                bacteria_slug = match.group(1)
                mechanism_slug = match.group(2)
                
                # Check if this mechanism should be excluded for this bacteria
                if bacteria_slug in exclusions and mechanism_slug in exclusions[bacteria_slug]:
                    # Comment out the line (preserve indentation)
                    indent = len(line) - len(line.lstrip())
                    commented_line = ' ' * indent + '// ' + line.lstrip()
                    modified_lines.append(commented_line)
                    commented_count += 1
                else:
                    modified_lines.append(line)
            else:
                modified_lines.append(line)
        else:
            modified_lines.append(line)
    
    return modified_lines, commented_count

def write_config_file(lines):
    """Write modified config.rs."""
    with open('src/config.rs', 'w', newline='') as f:
        f.writelines(lines)

def main():
    print("Step 1: Parsing BACTERIA_MECHANISM_EXCLUSIONS.md...")
    exclusions = parse_exclusions()
    
    total_exclusions = sum(len(mechs) for mechs in exclusions.values())
    print(f"  Found {len(exclusions)} bacteria with exclusions")
    print(f"  Total exclusions to apply: {total_exclusions}")
    
    # Show some examples
    print("\n  Examples:")
    for bacteria, mechs in list(exclusions.items())[:3]:
        if mechs:
            print(f"    {bacteria}: {len(mechs)} exclusions")
    
    print("\nStep 2: Processing config.rs...")
    modified_lines, commented_count = process_config_file(exclusions)
    print(f"  Commented out {commented_count} mechanism entries")
    
    print("\nStep 3: Writing updated config.rs...")
    write_config_file(modified_lines)
    
    print("\n✓ Successfully commented out biologically inapplicable mechanisms!")
    print(f"  {commented_count} lines now commented")
    print(f"  {1050 - commented_count} lines remain active")
    print("\n  Next: Run 'cargo check' to verify compilation")

if __name__ == '__main__':
    main()
