#!/usr/bin/env python3
"""
Comprehensive resistance mechanism audit
Maps drugs with observed resistance to their resistance mechanisms
"""

# All resistance mechanisms from src/simulation/population.rs
MECHANISMS = {
    "enzyme_esbl_ctx_m": {
        "name": "ESBL CTX-M",
        "type": "Beta-lactamase",
        "affects": ["penicillins", "cephalosporins (gen 1-3)"],
        "bacteria": ["E. coli", "Klebsiella", "Enterobacter", "Citrobacter", "Morganella", "Proteus", "Serratia"]
    },
    "enzyme_esbl_tem": {
        "name": "ESBL TEM",
        "type": "Beta-lactamase",
        "affects": ["penicillins", "cephalosporins (gen 1-3)"],
        "bacteria": ["E. coli", "Klebsiella", "Enterobacter", "H. influenzae", "N. gonorrhoeae"]
    },
    "enzyme_esbl_shv": {
        "name": "ESBL SHV",
        "type": "Beta-lactamase",
        "affects": ["penicillins", "cephalosporins (gen 1-3)"],
        "bacteria": ["Klebsiella", "E. coli", "Enterobacter"]
    },
    "enzyme_kpc": {
        "name": "KPC Carbapenemase",
        "type": "Carbapenemase",
        "affects": ["all beta-lactams including carbapenems"],
        "bacteria": ["Klebsiella", "E. coli", "Enterobacter", "Citrobacter", "Pseudomonas"]
    },
    "enzyme_ndm_vim": {
        "name": "NDM/VIM Carbapenemase",
        "type": "Metallo-beta-lactamase",
        "affects": ["all beta-lactams including carbapenems"],
        "bacteria": ["Gram-negatives incl. Klebsiella, E. coli, Acinetobacter, Pseudomonas"]
    },
    "enzyme_oxa_48": {
        "name": "OXA-48 Carbapenemase",
        "type": "Carbapenemase",
        "affects": ["carbapenems"],
        "bacteria": ["Enterobacterales, Acinetobacter"]
    },
    "enzyme_ampc_cmy": {
        "name": "AmpC CMY",
        "type": "AmpC Beta-lactamase",
        "affects": ["penicillins", "cephalosporins (gen 1-3)", "beta-lactam/inhibitor combos"],
        "bacteria": ["Enterobacter", "Citrobacter", "Serratia", "Morganella", "E. coli (plasmid)"]
    },
    "enzyme_ampc_dha": {
        "name": "AmpC DHA",
        "type": "AmpC Beta-lactamase",
        "affects": ["penicillins", "cephalosporins (gen 1-3)", "beta-lactam/inhibitor combos"],
        "bacteria": ["Enterobacterales"]
    },
    "target_site_pbp2a_meca": {
        "name": "PBP2a (mecA)",
        "type": "Target modification",
        "affects": ["all beta-lactams"],
        "bacteria": ["S. aureus (MRSA)", "S. epidermidis"]
    },
    "target_site_van_a": {
        "name": "VanA",
        "type": "Target modification",
        "affects": ["vancomycin", "teicoplanin", "dalbavancin"],
        "bacteria": ["Enterococcus"]
    },
    "target_site_van_b": {
        "name": "VanB",
        "type": "Target modification",
        "affects": ["vancomycin"],
        "bacteria": ["Enterococcus"]
    },
    "mutation_gyra_primary": {
        "name": "gyrA Primary Mutation",
        "type": "Chromosomal mutation",
        "affects": ["fluoroquinolones (all)"],
        "bacteria": ["Most bacteria"]
    },
    "mutation_gyra_parc_secondary": {
        "name": "gyrA/parC Secondary Mutation",
        "type": "Chromosomal mutation",
        "affects": ["fluoroquinolones (high-level)"],
        "bacteria": ["Most bacteria"]
    },
    "protection_qnr": {
        "name": "Qnr Protection",
        "type": "Plasmid-mediated quinolone resistance",
        "affects": ["fluoroquinolones (low-level)"],
        "bacteria": ["Enterobacterales"]
    },
    "enzyme_16s_rrmt": {
        "name": "16S rRNA Methyltransferase",
        "type": "Ribosomal modification",
        "affects": ["aminoglycosides (all)"],
        "bacteria": ["Gram-negatives"]
    },
    "target_site_erm_b": {
        "name": "ErmB",
        "type": "Ribosomal methylation",
        "affects": ["macrolides", "lincosamides (MLSB phenotype)"],
        "bacteria": ["Gram-positives (Strep, Staph, Enterococcus)"]
    },
    "target_site_cfr": {
        "name": "Cfr",
        "type": "Ribosomal methylation",
        "affects": ["linezolid", "chloramphenicol"],
        "bacteria": ["Staphylococcus", "Enterococcus"]
    },
    "enzyme_cat": {
        "name": "Chloramphenicol Acetyltransferase",
        "type": "Enzymatic inactivation",
        "affects": ["chloramphenicol"],
        "bacteria": ["Various"]
    },
    "efflux_acrab_tolc": {
        "name": "AcrAB-TolC Efflux",
        "type": "Efflux pump",
        "affects": ["multiple drug classes"],
        "bacteria": ["E. coli", "Enterobacterales"]
    },
    "efflux_mexxy_oprm": {
        "name": "MexXY-OprM Efflux",
        "type": "Efflux pump",
        "affects": ["aminoglycosides", "fluoroquinolones"],
        "bacteria": ["Pseudomonas"]
    },
    "porin_loss_ompk35_36": {
        "name": "OmpK35/36 Porin Loss",
        "type": "Reduced permeability",
        "affects": ["carbapenems", "beta-lactams"],
        "bacteria": ["Klebsiella"]
    },
    "porin_loss_oprd": {
        "name": "OprD Porin Loss",
        "type": "Reduced permeability",
        "affects": ["carbapenems (esp. imipenem)"],
        "bacteria": ["Pseudomonas"]
    },
    "modification_mcr_1": {
        "name": "MCR-1",
        "type": "Lipid A modification",
        "affects": ["colistin", "polymyxins"],
        "bacteria": ["Enterobacterales"]
    },
    "global_efflux_pump": {
        "name": "Global Efflux Pump",
        "type": "Multi-drug efflux",
        "affects": ["multiple classes"],
        "bacteria": ["Various"]
    },
    "global_porin_loss": {
        "name": "Global Porin Loss",
        "type": "Reduced permeability",
        "affects": ["multiple beta-lactams"],
        "bacteria": ["Gram-negatives"]
    }
}

# Drug class to mechanism mapping
DRUG_MECHANISM_MAP = {
    # Beta-lactams - Penicillins
    "penicillin_g": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_ampc_cmy", "enzyme_ampc_dha", "target_site_pbp2a_meca"],
    "ampicillin": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_ampc_cmy", "enzyme_ampc_dha", "target_site_pbp2a_meca"],
    "amoxicillin": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_ampc_cmy", "enzyme_ampc_dha", "target_site_pbp2a_meca"],
    "piperacillin": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_ampc_cmy", "enzyme_ampc_dha", "enzyme_kpc", "enzyme_ndm_vim"],
    "ticarcillin": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_ampc_cmy", "enzyme_ampc_dha"],
    
    # Beta-lactam/inhibitor combinations
    "amoxicillin_clavulanate": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_ampc_cmy", "enzyme_ampc_dha", "target_site_pbp2a_meca"],
    "ampicillin_sulbactam": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_ampc_cmy", "enzyme_ampc_dha", "target_site_pbp2a_meca"],
    "piperacillin_tazobactam": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_kpc", "enzyme_ndm_vim"],
    "ticarcillin_clavulanate": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_ampc_cmy", "enzyme_ampc_dha"],
    
    # Cephalosporins
    "cefazolin": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_ampc_cmy", "enzyme_ampc_dha", "target_site_pbp2a_meca"],
    "cephalexin": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_ampc_cmy", "enzyme_ampc_dha", "target_site_pbp2a_meca"],
    "cefuroxime": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_ampc_cmy", "enzyme_ampc_dha", "target_site_pbp2a_meca"],
    "ceftriaxone": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_kpc", "enzyme_ndm_vim"],
    "ceftazidime": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_kpc", "enzyme_ndm_vim", "efflux_mexxy_oprm"],
    "cefepime": ["enzyme_kpc", "enzyme_ndm_vim", "enzyme_oxa_48", "efflux_mexxy_oprm"],
    "ceftaroline": ["enzyme_kpc", "enzyme_ndm_vim"],
    "ceftazidime_avibactam": ["enzyme_ndm_vim", "efflux_mexxy_oprm"],
    
    # Monobactams
    "aztreonam": ["enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv", "enzyme_kpc", "enzyme_oxa_48"],
    
    # Carbapenems
    "imipenem_c": ["enzyme_kpc", "enzyme_ndm_vim", "enzyme_oxa_48", "porin_loss_ompk35_36", "porin_loss_oprd", "efflux_mexxy_oprm"],
    "meropenem": ["enzyme_kpc", "enzyme_ndm_vim", "enzyme_oxa_48", "porin_loss_ompk35_36", "porin_loss_oprd", "efflux_mexxy_oprm"],
    "ertapenem": ["enzyme_kpc", "enzyme_ndm_vim", "enzyme_oxa_48", "porin_loss_ompk35_36"],
    "meropenem_vaborbactam": ["enzyme_ndm_vim", "enzyme_oxa_48", "porin_loss_oprd"],
    
    # Fluoroquinolones
    "ciprofloxacin": ["mutation_gyra_primary", "mutation_gyra_parc_secondary", "protection_qnr", "efflux_acrab_tolc", "efflux_mexxy_oprm"],
    "levofloxacin": ["mutation_gyra_primary", "mutation_gyra_parc_secondary", "protection_qnr", "efflux_acrab_tolc", "efflux_mexxy_oprm"],
    "moxifloxacin": ["mutation_gyra_primary", "mutation_gyra_parc_secondary", "protection_qnr", "efflux_acrab_tolc"],
    "ofloxacin": ["mutation_gyra_primary", "mutation_gyra_parc_secondary", "protection_qnr", "efflux_acrab_tolc"],
    
    # Aminoglycosides
    "gentamicin": ["enzyme_16s_rrmt", "efflux_mexxy_oprm"],
    "tobramycin": ["enzyme_16s_rrmt", "efflux_mexxy_oprm"],
    "amikacin": ["enzyme_16s_rrmt", "efflux_mexxy_oprm"],
    
    # Macrolides
    "erythromycin": ["target_site_erm_b", "efflux_acrab_tolc"],
    "azithromycin": ["target_site_erm_b", "efflux_acrab_tolc"],
    "clarithromycin": ["target_site_erm_b", "efflux_acrab_tolc"],
    
    # Lincosamides
    "clindamycin": ["target_site_erm_b"],
    
    # Tetracyclines
    "tetracycline": ["efflux_acrab_tolc", "global_efflux_pump"],
    "doxycycline": ["efflux_acrab_tolc", "global_efflux_pump"],
    "minocycline": ["efflux_acrab_tolc", "global_efflux_pump"],
    
    # Glycopeptides
    "vancomycin": ["target_site_van_a", "target_site_van_b"],
    "teicoplanin": ["target_site_van_a"],
    "dalbavancin": ["target_site_van_a"],
    
    # Oxazolidinones
    "linezolid": ["target_site_cfr"],
    "tedizolid": ["target_site_cfr"],
    
    # Polymyxins
    "colistin": ["modification_mcr_1"],
    
    # Others
    "chloramphenicol": ["enzyme_cat", "target_site_cfr"],
    "rifampicin": [],  # TB-specific, rpoB mutations not in current mechanism list
    "trim_sulf": [],  # Sulfonamide/DHFR mutations not in current mechanism list
    "sulfanilamide": [],  # Sulfonamide mutations not in current mechanism list
    "nitrofurantoin": [],  # Nitrofuran resistance mechanisms not in current mechanism list
}

def print_audit_report():
    """Print comprehensive audit report"""
    
    print("="*100)
    print("RESISTANCE MECHANISM AUDIT REPORT")
    print("="*100)
    print()
    
    # Drugs with resistance from calibration file (from previous parse)
    drugs_with_resistance = {
        "amikacin": 6, "amoxicillin": 14, "amoxicillin_clavulanate": 16, "ampicillin": 14,
        "ampicillin_sulbactam": 15, "azithromycin": 9, "aztreonam": 12, "cefazolin": 13,
        "cefepime": 12, "ceftaroline": 10, "ceftazidime": 15, "ceftazidime_avibactam": 13,
        "ceftriaxone": 15, "cefuroxime": 13, "cephalexin": 12, "chloramphenicol": 17,
        "ciprofloxacin": 20, "clarithromycin": 9, "clindamycin": 7, "colistin": 2,
        "dalbavancin": 1, "doxycycline": 9, "ertapenem": 12, "erythromycin": 9,
        "gentamicin": 7, "imipenem_c": 11, "levofloxacin": 19, "linezolid": 2,
        "meropenem": 15, "meropenem_vaborbactam": 11, "minocycline": 3, "moxifloxacin": 14,
        "nitrofurantoin": 1, "ofloxacin": 18, "penicillin_g": 6, "piperacillin": 14,
        "piperacillin_tazobactam": 14, "rifampicin": 2, "sulfanilamide": 1, "tedizolid": 1,
        "teicoplanin": 2, "tetracycline": 14, "ticarcillin": 14, "ticarcillin_clavulanate": 14,
        "tobramycin": 7, "trim_sulf": 3, "vancomycin": 2
    }
    
    # Categorize drugs
    drugs_with_mechanisms = []
    drugs_without_mechanisms = []
    
    for drug, bacteria_count in sorted(drugs_with_resistance.items()):
        mechanisms = DRUG_MECHANISM_MAP.get(drug, [])
        if mechanisms:
            drugs_with_mechanisms.append((drug, bacteria_count, mechanisms))
        else:
            drugs_without_mechanisms.append((drug, bacteria_count))
    
    # Print summary
    print("EXECUTIVE SUMMARY")
    print("-"*100)
    print(f"Total drugs with observed resistance: {len(drugs_with_resistance)}")
    print(f"Drugs with defined mechanisms: {len(drugs_with_mechanisms)} ({len(drugs_with_mechanisms)/len(drugs_with_resistance)*100:.1f}%)")
    print(f"Drugs WITHOUT mechanisms (GAPS): {len(drugs_without_mechanisms)} ({len(drugs_without_mechanisms)/len(drugs_with_resistance)*100:.1f}%)")
    print()
    
    # Print drugs WITH mechanisms
    print("DRUGS WITH RESISTANCE MECHANISMS DEFINED")
    print("-"*100)
    for drug, bacteria_count, mechanisms in drugs_with_mechanisms:
        print(f"\n✓ {drug.upper()}")
        print(f"  Bacteria with resistance: {bacteria_count}")
        print(f"  Mechanisms ({len(mechanisms)}):")
        for mech in mechanisms:
            info = MECHANISMS.get(mech, {})
            print(f"    - {info.get('name', mech)} ({info.get('type', 'unknown')})")
    
    # Print drugs WITHOUT mechanisms (CRITICAL GAPS)
    print("\n\n")
    print("="*100)
    print("CRITICAL GAPS - DRUGS WITH RESISTANCE BUT NO MECHANISMS")
    print("="*100)
    
    for drug, bacteria_count in drugs_without_mechanisms:
        print(f"\n✗ {drug.upper()}")
        print(f"  Bacteria showing resistance: {bacteria_count}")
        print(f"  STATUS: NO MECHANISM DEFINED")
        
        # Provide recommendations
        if drug in ["rifampicin"]:
            print("  RECOMMENDATION: Add rpoB mutation mechanism for TB resistance")
        elif drug in ["trim_sulf", "sulfanilamide"]:
            print("  RECOMMENDATION: Add folP/DHFR mutation mechanisms for sulfonamide/trimethoprim resistance")
        elif drug in ["nitrofurantoin"]:
            print("  RECOMMENDATION: Add nfsA/nfsB mutations or efflux mechanisms for nitrofuran resistance")
    
    print("\n" + "="*100)
    print("END OF AUDIT REPORT")
    print("="*100)

if __name__ == "__main__":
    print_audit_report()
