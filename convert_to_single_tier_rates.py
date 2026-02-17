#!/usr/bin/env python3
"""
Convert two-tier emergence rates (mechanism base rate × bacteria-mechanism multiplier)
to single-tier bacteria-mechanism emergence rates.
"""

import re

# Mechanism base rates from config.rs lines 7375-7405
MECHANISM_BASE_RATES = {
    "enzyme_esbl_ctx_m": 0.00001,
    "enzyme_esbl_tem": 0.00001,
    "enzyme_esbl_shv": 0.00001,
    "enzyme_ampc_cmy": 0.00005,
    "enzyme_ampc_dha": 0.000005,
    "enzyme_kpc": 0.000001,
    "enzyme_ndm_vim": 0.0000005,
    "enzyme_oxa_48": 0.000001,
    "target_site_pbp2a_meca": 0.000005,
    "target_site_van_a": 0.000001,
    "target_site_van_b": 0.000001,
    "mutation_gyra_primary": 0.00005,
    "mutation_gyra_parc_secondary": 0.000005,
    "target_site_erm_b": 0.00002,
    "target_site_cfr": 0.000001,
    "protection_qnr": 0.00005,
    "enzyme_16s_rrmt": 0.00002,
    "enzyme_cat": 0.00001,
    "modification_mcr_1": 0.000001,
    "efflux_acrab_tolc": 0.00005,
    "efflux_mexxy_oprm": 0.00005,
    "porin_loss_ompk35_36": 0.00002,
    "porin_loss_oprd": 0.00005,
    "global_efflux_pump": 0.00001,
    "global_porin_loss": 0.00001,
}

def convert_config_file(input_file, output_file):
    """Read config.rs, convert multipliers to rates, write to output."""
    
    with open(input_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Pattern to match bacteria-mechanism multiplier lines
    # Example: map.insert("bacteria_e_coli_mechanism_enzyme_esbl_tem_emergence_multiplier".to_string(), 5.0);
    pattern = r'map\.insert\("bacteria_([^"]+)_mechanism_([^"]+)_emergence_multiplier"\.to_string\(\), ([0-9.]+)\);'
    
    def replace_multiplier(match):
        bacteria = match.group(1)
        mechanism = match.group(2)
        multiplier = float(match.group(3))
        
        # Look up base rate for this mechanism
        if mechanism in MECHANISM_BASE_RATES:
            base_rate = MECHANISM_BASE_RATES[mechanism]
            new_rate = base_rate * multiplier
            
            # Format with appropriate precision (use scientific notation if < 0.001)
            if new_rate < 0.001:
                rate_str = f"{new_rate:.6e}"
            else:
                rate_str = f"{new_rate}"
            
            # Return new format with "_rate" instead of "_multiplier"
            return f'map.insert("bacteria_{bacteria}_mechanism_{mechanism}_emergence_rate".to_string(), {rate_str});'
        else:
            # Mechanism not found, return unchanged
            print(f"WARNING: Unknown mechanism '{mechanism}' for bacteria '{bacteria}'")
            return match.group(0)
    
    # Replace all multiplier lines with rate lines
    converted_content = re.sub(pattern, replace_multiplier, content)
    
    # Count conversions
    multiplier_count = len(re.findall(pattern, content))
    print(f"Converted {multiplier_count} bacteria-mechanism multipliers to rates")
    
    # Write output
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(converted_content)
    
    print(f"Wrote converted config to: {output_file}")

if __name__ == "__main__":
    convert_config_file(
        "src/config.rs",
        "src/config.rs.converted"
    )
