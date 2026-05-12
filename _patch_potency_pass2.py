"""
Pass 2A: Raise potency in config.rs for drugs that ARE clinically used for
the listed organisms but were incorrectly set to the default 0.10 (below
the 0.15 selection threshold), making resistance in those pairs unscoreable.
Each old value is verified before replacement to catch stale assumptions.
"""

PATH = "src/config.rs"

# key -> (expected_old, new_value, rationale)
RAISES = {
    # TIGECYCLINE — FDA/EMA approved for complicated intra-abdominal, skin, CAP;
    # EUCAST clinical breakpoints exist for all organisms below.
    "drug_tigecycline_for_bacteria_acinetobacter_baumannii_potency_when_no_r":      (0.10, 0.55, "last-resort MDR; EUCAST S≤1 R>2 mg/L"),
    "drug_tigecycline_for_bacteria_bacteroides_fragilis_potency_when_no_r":         (0.10, 0.65, "good anaerobic spectrum; FDA approved complicated intra-abdominal"),
    "drug_tigecycline_for_bacteria_citrobacter_spp._potency_when_no_r":             (0.10, 0.55, "Enterobacterales; EUCAST breakpoints"),
    "drug_tigecycline_for_bacteria_enterobacter_cloacae_potency_when_no_r":         (0.10, 0.55, "Enterobacterales; EUCAST breakpoints"),
    "drug_tigecycline_for_bacteria_enterobacter_spp._potency_when_no_r":            (0.10, 0.55, "Enterobacterales; EUCAST breakpoints"),
    "drug_tigecycline_for_bacteria_enterococcus_faecalis_potency_when_no_r":        (0.10, 0.65, "VRE/VSE infections; EUCAST breakpoints"),
    "drug_tigecycline_for_bacteria_enterococcus_faecium_potency_when_no_r":         (0.10, 0.65, "VRE treatment option; EUCAST breakpoints"),
    "drug_tigecycline_for_bacteria_escherichia_coli_potency_when_no_r":             (0.10, 0.60, "FDA approved; primary Enterobacteriaceae indication"),
    "drug_tigecycline_for_bacteria_haemophilus_influenzae_potency_when_no_r":       (0.10, 0.50, "EUCAST breakpoints exist; moderate clinical use"),
    "drug_tigecycline_for_bacteria_klebsiella_pneumoniae_potency_when_no_r":        (0.10, 0.60, "FDA approved; major Enterobacteriaceae indication"),
    "drug_tigecycline_for_bacteria_moraxella_catarrhalis_potency_when_no_r":        (0.10, 0.50, "EUCAST breakpoints; used in mixed RTI"),
    "drug_tigecycline_for_bacteria_staphylococcus_aureus_potency_when_no_r":        (0.10, 0.65, "FDA approved SSTI; active against MRSA"),
    "drug_tigecycline_for_bacteria_streptococcus_agalactiae_potency_when_no_r":     (0.10, 0.50, "EUCAST breakpoints; used in complicated infections"),
    "drug_tigecycline_for_bacteria_streptococcus_pneumoniae_potency_when_no_r":     (0.10, 0.55, "FDA approved CAP; EUCAST breakpoints"),
    "drug_tigecycline_for_bacteria_streptococcus_pyogenes_potency_when_no_r":       (0.10, 0.50, "EUCAST breakpoints; used in complicated SSTI"),
    # CEFTAROLINE — 5th-gen cephalosporin; good activity against non-AmpC
    # Enterobacteriaceae; EUCAST clinical breakpoints for Enterobacterales.
    # AmpC-producing ESKAPE organisms (Enterobacter, Citrobacter, Morganella,
    # Serratia) have intrinsic AmpC resistance → potency stays at 0.10 → dot CSV.
    "drug_ceftaroline_for_bacteria_escherichia_coli_potency_when_no_r":             (0.10, 0.65, "5th-gen ceph; EUCAST S≤0.5 R>1 mg/L; used in complicated SSTI"),
    "drug_ceftaroline_for_bacteria_klebsiella_pneumoniae_potency_when_no_r":        (0.10, 0.55, "EUCAST breakpoints; susceptible non-ESBL strains"),
    # CEFIDEROCOL — FDA/EMA approved for Gram-negative infections including
    # Acinetobacter (FDA 2019) and Pseudomonas.
    "drug_cefiderocol_for_bacteria_acinetobacter_baumannii_potency_when_no_r":      (0.10, 0.55, "FDA approved for Acinetobacter; siderophore cephalosporin"),
    "drug_cefiderocol_for_bacteria_pseudomonas_aeruginosa_potency_when_no_r":       (0.10, 0.55, "FDA approved; active against MDR Pseudomonas"),
    # CEFTOLOZANE/TAZOBACTAM — primary indication MDR Pseudomonas aeruginosa.
    "drug_ceftolozane_tazobactam_for_bacteria_pseudomonas_aeruginosa_potency_when_no_r": (0.10, 0.65, "primary MDR Pseudomonas indication; EUCAST breakpoints"),
    # CIPROFLOXACIN — used for Mycoplasma genitalium (STI treatment despite
    # high resistance; still prescribed per guidelines).
    "drug_ciprofloxacin_for_bacteria_mycoplasma_genitalium_potency_when_no_r":      (0.10, 0.55, "used in STI treatment; high acquired resistance well-documented"),
}

with open(PATH, "r", encoding="utf-8") as f:
    content = f.read()

changed = 0
for key, (expected_old, new_val, rationale) in RAISES.items():
    old_str = f'map.insert("{key}".to_string(), {expected_old:.2f});'
    new_str = f'map.insert("{key}".to_string(), {new_val:.2f}); // raised: {rationale}'
    if old_str not in content:
        print(f"  MISS (not found or already changed): {key}")
        continue
    content = content.replace(old_str, new_str, 1)
    changed += 1
    short = key.split("drug_")[1].split("_potency")[0]
    print(f"  {short}: {expected_old:.2f} -> {new_val:.2f}")

print(f"\nTotal: {changed}/{len(RAISES)} keys changed.")
with open(PATH, "w", encoding="utf-8") as f:
    f.write(content)
print("Written to", PATH)
