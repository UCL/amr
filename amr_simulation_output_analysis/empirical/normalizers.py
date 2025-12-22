#!/usr/bin/env python3
"""
Name normalization functionality for empirical data matching.

This module contains functions for normalizing bacteria and drug names
for matching between simulation and empirical data sources.
"""

def normalize_name_for_empirical_matching(name, entity_type='bacteria', data_source=None):
    """
    Normalize bacteria/drug names for matching between simulation and empirical data.
    
    Args:
        name: The name to normalize
        entity_type: 'bacteria' or 'drug' - determines normalization strategy
        data_source: 'drug_failure', 'mic_values', 'resistance', etc. - specific empirical data type
    
    For bacteria: Handle mixed underscore/space usage in empirical data
    For drugs: Both simulation and empirical use underscores for combination drugs.
    """
    if name is None:
        return None
    
    if entity_type == 'drug':
        # For drugs, handle region-prefixed names (e.g., 'europe_penicilling' -> 'penicilling')
        # and keep underscores as empirical data uses them for combinations
        
        # Extract base drug name from region-prefixed names
        regions = ['north_america', 'south_america', 'europe', 'asia', 'africa', 'oceania']
        normalized = name
        
        for region in regions:
            if name.startswith(f"{region}_"):
                normalized = name[len(f"{region}_"):]
                break
        
        # Handle specific drug name variations if needed
        drug_mappings = {
            # Add any specific drug name mappings here if needed
        }
        
        if normalized in drug_mappings:
            return drug_mappings[normalized]
        
        return normalized
    
    else:  # bacteria
        # CONTEXT-AWARE BACTERIA NORMALIZATION
        # Different empirical datasets use different naming conventions
        
        if data_source == 'drug_failure':
            # Drug failure data uses abbreviated names (s_aureus, e_coli, etc.)
            simulation_to_drug_failure_bacteria = {
                'staphylococcus_aureus': 's_aureus',
                'escherichia_coli': 'e_coli',
                'klebsiella_pneumoniae': 'k_pneumoniae', 
                'pseudomonas_aeruginosa': 'p_aeruginosa',
                'streptococcus_pneumoniae': 's_pneumoniae',
                'acinetobacter_baumannii': 'a_baumannii',
                'enterococcus_faecalis': 'enterococcus_faecalis',  # These stay full
                'enterococcus_faecium': 'enterococcus_faecium'
            }
            
            if name in simulation_to_drug_failure_bacteria:
                return simulation_to_drug_failure_bacteria[name]
            else:
                return name  # Return as-is if no mapping found
                
        else:
            # For other data sources (MIC, resistance, etc.) - use space/underscore mappings
            simulation_to_empirical_bacteria = {
                # Space/underscore mappings for resistance empirical data
                'treponema_pallidum': 'treponema pallidum',
                'acinetobacter_baumannii': 'acinetobacter baumannii',
                'haemophilus_influenzae': 'haemophilus influenzae',
                'chlamydia_trachomatis': 'chlamydia trachomatis',
                'enterococcus_faecalis': 'enterococcus faecalis',
                'enterococcus_faecium': 'enterococcus faecium',
                'escherichia_coli': 'escherichia coli',
                'klebsiella_pneumoniae': 'klebsiella pneumoniae',
                'mdr_mycobacterium_tuberculosis': 'mdr mycobacterium tuberculosis',
                'neisseria_gonorrhoeae': 'neisseria gonorrhoeae',
                'neisseria_meningitidis': 'neisseria meningitidis',
                'pseudomonas_aeruginosa': 'pseudomonas aeruginosa',
                'salmonella_enterica_serovar_paratyphi_a': 'salmonella enterica serovar paratyphi a',
                'salmonella_enterica_serovar_typhi': 'salmonella enterica serovar typhi',
                'staphylococcus_aureus': 'staphylococcus aureus',
                'streptococcus_agalactiae': 'streptococcus agalactiae',
                'streptococcus_pneumoniae': 'streptococcus pneumoniae',
                'streptococcus_pyogenes': 'streptococcus pyogenes',
                'vibrio_cholerae': 'vibrio cholerae',
                'yersinia_enterocolitica': 'yersinia enterocolitica',
                # Add missing bacteria mappings
                'enterobacter_spp.': 'enterobacter spp.',
                'invasive_non_typhoidal_salmonella_spp': 'invasive non-typhoidal salmonella spp.',
                'citrobacter_spp.': 'citrobacter spp.',
                'morganella_spp.': 'morganella spp.',
                'proteus_spp.': 'proteus spp.',
                'serratia_spp.': 'serratia spp.',
                'shigella_spp.': 'shigella spp.'
            }
            
            # Check if we have a direct mapping
            if name in simulation_to_empirical_bacteria:
                return simulation_to_empirical_bacteria[name]
            
            # For most cases, empirical data uses same format as simulation (underscores)
            return name