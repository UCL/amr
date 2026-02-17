#!/usr/bin/env python3
"""
Generate standardized resistance mechanism blocks for all 42 bacteria in config.rs

This script generates all 25 resistance mechanisms in standardized biological order
for all 42 bacteria, preserving custom multiplier values where they differ from 1000.0.
"""

import re

# Standard biological order for all 25 mechanisms
STANDARD_MECHANISM_ORDER = [
    # ESBL Enzymes (3)
    "enzyme_esbl_ctx_m",
    "enzyme_esbl_tem",
    "enzyme_esbl_shv",
    # AmpC Enzymes (2)
    "enzyme_ampc_cmy",
    "enzyme_ampc_dha",
    # Carbapenemases (3)
    "enzyme_kpc",
    "enzyme_ndm_vim",
    "enzyme_oxa_48",
    # Other Enzymatic (2)
    "enzyme_cat",
    "enzyme_16s_rrmt",
    # Target Sites (5)
    "target_site_pbp2a_meca",
    "target_site_van_a",
    "target_site_van_b",
    "target_site_erm_b",
    "target_site_cfr",
    # FQ Resistance (3)
    "mutation_gyra_primary",
    "mutation_gyra_parc_secondary",
    "protection_qnr",
    # Efflux (3)
    "efflux_acrab_tolc",
    "efflux_mexxy_oprm",
    "global_efflux_pump",
    # Porin (3)
    "porin_loss_ompk35_36",
    "porin_loss_oprd",
    "global_porin_loss",
    # Surface (1)
    "modification_mcr_1",
]

# List of all 42 bacteria (from BACTERIA_LIST in population.rs)
BACTERIA_LIST = [
    "acinetobacter_baumannii",
    "citrobacter_spp.",
    "enterobacter_spp.",
    "enterococcus_faecalis",
    "enterococcus_faecium",
    "escherichia_coli",
    "klebsiella_pneumoniae",
    "morganella_spp.",
    "proteus_spp.",
    "serratia_spp.",
    "p_stuartii",
    "pseudomonas_aeruginosa",
    "stenotrophomonas_maltophilia",
    "staphylococcus_aureus",
    "staphylococcus_epidermidis",
    "streptococcus_pneumoniae",
    "salmonella_enterica_serovar_typhi",
    "salmonella_enterica_serovar_paratyphi_a",
    "invasive_non-typhoidal_salmonella_spp.",
    "shigella_spp.",
    "neisseria_gonorrhoeae",
    "streptococcus_pyogenes",
    "streptococcus_agalactiae",
    "haemophilus_influenzae",
    "chlamydia_trachomatis",
    "mycoplasma_genitalium",
    "vibrio_cholerae",
    "neisseria_meningitidis",
    "listeria_monocytogenes",
    "clostridioides_difficile",
    "bacteroides_fragilis",
    "campylobacter_jejuni",
    "enterobacter_cloacae",
    "yersinia_enterocolitica",
    "moraxella_catarrhalis",
    "treponema_pallidum",
    "bordetella_pertussis",
    "helicobacter_pylori",
    "mdr_mycobacterium_tuberculosis",
    "mycoplasma_pneumoniae",
    "legionella_pneumophila",
    "burkholderia_cepacia_complex",
]

# Custom multipliers for specific bacteria (extracted from existing config.rs)
# Format: {bacteria_name: {mechanism_name: value}}
CUSTOM_MULTIPLIERS = {
    "escherichia_coli": {
        # E. coli has special values of 1.0 for many mechanisms
        "mutation_gyra_primary": 1.0,
        "efflux_acrab_tolc": 1.0,
        "porin_loss_ompk35_36": 1.0,
        "protection_qnr": 1.0,
        "target_site_erm_b": 1.0,
        "enzyme_esbl_ctx_m": 1.0,
        "enzyme_ampc_cmy": 1.0,
        "target_site_pbp2a_meca": 1.0,
        "enzyme_kpc": 1.0,
        "target_site_van_a": 1.0,
        "enzyme_16s_rrmt": 1.0,
        "modification_mcr_1": 1.0,
        "enzyme_cat": 1.0,
        "enzyme_oxa_48": 1.0,
        "enzyme_ndm_vim": 1.0,
    },
    "klebsiella_pneumoniae": {
        "mutation_gyra_primary": 10.0,
        "protection_qnr": 10.0,
        "enzyme_esbl_ctx_m": 10.0,
        "enzyme_kpc": 10.0,
    },
    "neisseria_gonorrhoeae": {
        # All mechanisms 5000.0
        **{mech: 5000.0 for mech in STANDARD_MECHANISM_ORDER}
    },
    "mycoplasma_genitalium": {
        **{mech: 10000.0 for mech in STANDARD_MECHANISM_ORDER}
    },
    "mdr_mycobacterium_tuberculosis": {
        "mutation_gyra_primary": 300000.0,
        "protection_qnr": 5000.0,
    },
    "moraxella_catarrhalis": {
        "mutation_gyra_primary": 0.3,
        "efflux_acrab_tolc": 0.3,
        "porin_loss_ompk35_36": 0.3,
        "enzyme_esbl_ctx_m": 100.0,
    },
    "mycoplasma_pneumoniae": {
        "mutation_gyra_primary": 3.0,
        "target_site_erm_b": 3.0,
        "efflux_acrab_tolc": 0.5,
        "porin_loss_ompk35_36": 0.5,
    },
    "legionella_pneumophila": {
        "efflux_acrab_tolc": 1.0,
        "porin_loss_ompk35_36": 1.0,
    },
    "listeria_monocytogenes": {
        "efflux_acrab_tolc": 0.2,
    },
    "burkholderia_cepacia_complex": {
        "mutation_gyra_primary": 50.0,
        "efflux_acrab_tolc": 100.0,
        "porin_loss_ompk35_36": 50.0,
        "enzyme_esbl_ctx_m": 20.0,
        "enzyme_kpc": 10.0,
    },
    "morganella_spp.": {
        "mutation_gyra_primary": 100.0,
        "efflux_acrab_tolc": 100.0,
        "porin_loss_ompk35_36": 100.0,
        "enzyme_16s_rrmt": 100.0,
        "modification_mcr_1": 100.0,
    },
    "serratia_spp.": {
        **{mech: 50.0 for mech in ["mutation_gyra_primary", "efflux_acrab_tolc", 
                                     "porin_loss_ompk35_36", "protection_qnr", 
                                     "enzyme_esbl_ctx_m", "enzyme_ampc_cmy", 
                                     "enzyme_kpc", "enzyme_16s_rrmt"]}
    },
    "pseudomonas_aeruginosa": {
        "mutation_gyra_primary": 50.0,
        # Will add more as we read existing config
    },
}

def get_multiplier(bacteria: str, mechanism: str) -> float:
    """Get the emergence multiplier for a bacteria-mechanism pair."""
    if bacteria in CUSTOM_MULTIPLIERS and mechanism in CUSTOM_MULTIPLIERS[bacteria]:
        return CUSTOM_MULTIPLIERS[bacteria][mechanism]
    return 1000.0  # Default

def generate_mechanism_block(bacteria: str, comment: str = "") -> str:
    """Generate a complete mechanism block for one bacteria."""
    lines = []
    
    if comment:
        lines.append(f"\n        // {comment}")
    
    for mechanism in STANDARD_MECHANISM_ORDER:
        multiplier = get_multiplier(bacteria, mechanism)
        key = f"bacteria_{bacteria}_mechanism_{mechanism}_emergence_multiplier"
        lines.append(f'        map.insert("{key}".to_string(), {multiplier});')
    
    return "\n".join(lines)

def generate_all_blocks() -> str:
    """Generate mechanism blocks for all 42 bacteria."""
    output = []
    
    # Add header
    output.append("\n        // =====================================================================================")
    output.append("        // RESISTANCE MECHANISM EMERGENCE MULTIPLIERS")
    output.append("        // =====================================================================================")
    output.append("        // All 25 resistance mechanisms in standardized biological order for all 42 bacteria.")
    output.append("        // Default multiplier: 1000.0 (rare emergence)")
    output.append("        // Lower values (1.0-100.0) indicate higher propensity for that mechanism")
    output.append("        // =====================================================================================\n")
    
    # E. coli first (special case with 1.0 values)
    output.append(generate_mechanism_block(
        "escherichia_coli",
        "E. coli - Lower emergence multipliers reflect high clinical prevalence of resistance"
    ))
    
    # All other bacteria
    for bacteria in BACTERIA_LIST:
        if bacteria == "escherichia_coli":
            continue
        
        # Add descriptive comments for notable bacteria
        comments = {
            "klebsiella_pneumoniae": "Klebsiella pneumoniae - The 'Plasmid Sponge' with high ESBL/Carbapenemase acquisition",
            "neisseria_gonorrhoeae": "Neisseria gonorrhoeae - 'Superbug' potential with rapid resistance acquisition",
            "mycoplasma_genitalium": "Mycoplasma genitalium - High mutation rates in 23S (Macrolide) and ParC (FQ)",
            "mdr_mycobacterium_tuberculosis": "MDR Mycobacterium tuberculosis - Chromosomal mutations drive XDR",
            "pseudomonas_aeruginosa": "Pseudomonas aeruginosa - Multi-mechanism resistance under ICU pressure",
            "acinetobacter_baumannii": "Acinetobacter baumannii - Pan-drug resistance potential",
        }
        
        comment = comments.get(bacteria, "")
        output.append(generate_mechanism_block(bacteria, comment))
    
    return "\n".join(output)

if __name__ == "__main__":
    result = generate_all_blocks()
    print(result)
    print(f"\n\n// Generated {len(BACTERIA_LIST)} bacteria × {len(STANDARD_MECHANISM_ORDER)} mechanisms = {len(BACTERIA_LIST) * len(STANDARD_MECHANISM_ORDER)} lines")
