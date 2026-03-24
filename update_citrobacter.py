import re

with open('src/config.rs', 'r', encoding='utf-8') as f:
    text = f.read()

# Make transformations to Citrobacter rates:
# Penicillins/Cephs down 10x
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_enzyme_esbl_ctx_m_emergence_rate"\.to_string\(\),\s*)0\.000_3', r'\g<1>0.000_03', text)
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_enzyme_esbl_tem_emergence_rate"\.to_string\(\),\s*)0\.000_3', r'\g<1>0.000_03', text)
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_enzyme_esbl_shv_emergence_rate"\.to_string\(\),\s*)0\.000_1', r'\g<1>0.000_01', text)
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_enzyme_ampc_dha_emergence_rate"\.to_string\(\),\s*)0\.000_3', r'\g<1>0.000_03', text)

# Carbapenems down 10x
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_enzyme_kpc_emergence_rate"\.to_string\(\),\s*)0\.000_015', r'\g<1>0.000_001_5', text)
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_enzyme_ndm_vim_emergence_rate"\.to_string\(\),\s*)0\.000_015', r'\g<1>0.000_001_5', text)
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_enzyme_oxa_48_emergence_rate"\.to_string\(\),\s*)0\.000_015', r'\g<1>0.000_001_5', text)

# FQ boost
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_mutation_gyra_primary_emergence_rate"\.to_string\(\),\s*)0\.000_3', r'\g<1>0.004', text)
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_mutation_gyra_parc_secondary_emergence_rate"\.to_string\(\),\s*)0\.000_3', r'\g<1>0.004', text)
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_protection_qnr_emergence_rate"\.to_string\(\),\s*)0\.000_3', r'\g<1>0.002', text)

# Aminoglycosides boost
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_enzyme_aac_aph_emergence_rate"\.to_string\(\),\s*)0\.000_000_000_05', r'\g<1>0.004', text)
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_enzyme_16s_rrmt_emergence_rate"\.to_string\(\),\s*)0\.000_005', r'\g<1>0.002_5', text)

# Tetracyclines boost
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_protection_tet_m_emergence_rate"\.to_string\(\),\s*)0\.000_3', r'\g<1>0.003', text)
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_efflux_tet_abc_emergence_rate"\.to_string\(\),\s*)0\.000_01', r'\g<1>0.002_5', text)
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_efflux_acrab_tolc_emergence_rate"\.to_string\(\),\s*)0\.000_3', r'\g<1>0.002_5', text)

# Global Porin Loss / Global Efflux 
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_global_porin_loss_emergence_rate"\.to_string\(\),\s*)0\.000_05', r'\g<1>0.004', text)
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_global_efflux_pump_emergence_rate"\.to_string\(\),\s*)0\.000_03', r'\g<1>0.003', text)

# Polymyxins
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_modification_mcr_1_emergence_rate"\.to_string\(\),\s*)0\.000_03', r'\g<1>0.003', text)

# Folate pathway
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_mutation_folate_pathway_emergence_rate"\.to_string\(\),\s*)0\.000_3', r'\g<1>0.004', text)

# Nitrofurans
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_mutation_nitroreductase_emergence_rate"\.to_string\(\),\s*)0\.000_05', r'\g<1>0.002', text)

# Rifamycins
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_mutation_rpo_b_emergence_rate"\.to_string\(\),\s*)0\.000_005', r'\g<1>0.002', text)

# Chloramphenicol
text = re.sub(r'("bacteria_citrobacter_spp\._mechanism_enzyme_cat_emergence_rate"\.to_string\(\),\s*)0\.000_3', r'\g<1>0.003', text)

with open('src/config.rs', 'w', encoding='utf-8') as f:
    f.write(text)

print("Done")
