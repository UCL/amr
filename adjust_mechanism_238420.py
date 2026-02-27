"""
Mechanism-aware calibration adjustments based on run 238420 results.

Strategy:
  GROUP A — Revert 5 bacteria with stripped band multipliers to backup values,
            then apply per-bacterium multipliers:
    S. aureus:         revert to backup (no additional multiplier)
    P. aeruginosa:     revert to backup (no additional multiplier)
    N. gonorrhoeae:    revert to backup, then x1.5 on non-zero mechanisms
    E. faecium:        revert to backup, then x3 on non-zero mechanisms (cap 1e-3)
    S. pneumoniae:     revert to backup, then x0.7 on non-zero mechanisms

  GROUP B — Boost 9 "all-0%" bacteria with uniform multiplier on non-zero mechanisms:
    E. cloacae:        x1.5 (was at 45% before, stochastic collapse)
    All others:        x3

  GROUP C — Mechanism-selective adjustments for bacteria with partial resistance:
    K. pneumoniae:     lower carbapenemases x0.3, boost FQ/AG/Tet/etc x2-3
    A. baumannii:      lower beta-lactam/carbapenem mechs, boost TMP-SMX/rifampicin
    B. fragilis:       lower carbapenemases x0.3, boost FQ/AG/Tet/TMP-SMX
"""

import re
import math

TIER_VALUES = [
    0.0, 1e-9, 1e-8, 5e-8, 1e-7, 5e-7,
    1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3,
]
MAX_TIER = len(TIER_VALUES) - 1  # 12
MAX_RATE = 1e-3

CONFIG = r"src\config.rs"
BACKUP = r"src\config.rs.112750"

def value_to_tier(val):
    """Find the nearest tier for a given value."""
    if val <= 0:
        return 0
    best_tier = 0
    best_dist = float('inf')
    for i, tv in enumerate(TIER_VALUES):
        if tv == 0:
            continue
        dist = abs(math.log10(val) - math.log10(tv))
        if dist < best_dist:
            best_dist = dist
            best_tier = i
    return best_tier

def format_value(val):
    """Format a rate value consistently."""
    if val == 0.0:
        return "0.0"
    elif val >= 1e-3:
        return f"{val:.1e}"
    else:
        return f"{val:.1e}"

# ---------- Load backup file for GROUP A reverts ----------
with open(BACKUP, "r") as f:
    backup_lines = f.readlines()

# Build a lookup: (bacterium, mechanism) -> (value, tier) from backup
backup_pattern = re.compile(
    r'map\.insert\("bacteria_([\w.\-]+)_mechanism_([\w.]+)_emergence_rate"'
    r'\.to_string\(\),\s*([\d.eE\-+]+)\);\s*//\s*tier\s+(\d+)'
)
backup_lookup = {}
for line in backup_lines:
    m = backup_pattern.search(line)
    if m:
        bact = m.group(1)
        mech = m.group(2)
        val = float(m.group(3))
        tier = int(m.group(4))
        backup_lookup[(bact, mech)] = (val, tier)

# ---------- Define adjustment rules ----------

# GROUP A: Revert to backup with optional multiplier
GROUP_A = {
    "staphylococcus_aureus":       1.0,   # pure revert
    "pseudomonas_aeruginosa":      1.0,   # pure revert
    "neisseria_gonorrhoeae":       1.5,   # revert then x1.5
    "enterococcus_faecium":        3.0,   # revert then x3
    "streptococcus_pneumoniae":    0.7,   # revert then x0.7
}

# GROUP B: Multiply current values (uniform)
GROUP_B = {
    "enterobacter_cloacae":        1.5,
    "shigella_spp.":               3.0,
    "campylobacter_jejuni":        3.0,
    "invasive_non-typhoidal_salmonella_spp.": 3.0,
    "salmonella_enterica_serovar_paratyphi_a": 3.0,
    "enterobacter_spp.":           3.0,
    "citrobacter_spp.":            3.0,
    "serratia_spp.":               3.0,
    "morganella_spp.":             3.0,
    "escherichia_coli":            3.0,   # included — 22.78 Δ, all 0%, high PDs
}

# GROUP C: Per-mechanism multipliers (applied to current values)
GROUP_C = {
    "klebsiella_pneumoniae": {
        # Lower carbapenemases (sim 55% vs target 5-15%)
        "enzyme_kpc":                    0.3,
        "enzyme_ndm_vim":                0.3,
        "enzyme_oxa_48":                 0.3,
        # Boost FQ (targets 35-38%)
        "mutation_gyra_primary":         2.0,
        "mutation_gyra_parc_secondary":  2.0,
        "protection_qnr":               2.0,
        # Boost AG (targets 8-22%)
        "enzyme_16s_rrmt":               3.0,
        # Boost Tet/efflux (targets 25-35%)
        "protection_tet_m":              2.0,
        "efflux_acrab_tolc":             2.0,
        "global_efflux_pump":            2.0,
        # Boost TMP/SMX (target 55%)
        "mutation_folate_pathway":       2.0,
        # Boost Chloramphenicol (target 18%)
        "enzyme_cat":                    2.0,
        # Boost Colistin (target 8%)
        "modification_mcr_1":            2.0,
    },
    "acinetobacter_baumannii": {
        # Lower beta-lactam/carbapenem (sim 96% vs targets 45-70%)
        "enzyme_oxa_48":                 0.3,
        "enzyme_ampc_cmy":               0.5,
        "enzyme_ndm_vim":                0.5,
        "enzyme_kpc":                    0.5,
        "enzyme_esbl_ctx_m":             0.5,
        "enzyme_esbl_tem":               0.5,
        "enzyme_esbl_shv":               0.5,
        "enzyme_ampc_dha":               0.5,
        # Boost trim_sulf (sim 12.6% vs target 65%)
        "mutation_folate_pathway":       3.0,
        # Boost rifampicin (sim 12.6% vs target 55%)
        "mutation_rpo_b":                5.0,
    },
    "bacteroides_fragilis": {
        # Lower carbapenemases (sim 51.72% vs carbapenem targets 5-8%)
        "enzyme_ndm_vim":                0.3,
        "enzyme_kpc":                    0.3,
        "enzyme_oxa_48":                 0.3,
        # Boost ESBL/AmpC (sim 51.72% vs amox/amp targets 85-90%)
        "enzyme_esbl_ctx_m":             2.0,
        "enzyme_esbl_tem":               2.0,
        "enzyme_ampc_cmy":               1.5,
        # Boost FQ (0% vs target 35%)
        "mutation_gyra_primary":         3.0,
        "mutation_gyra_parc_secondary":  3.0,
        # Boost AG (0% vs targets 75-85%)
        "enzyme_16s_rrmt":               5.0,
        # Boost Chloramphenicol (0% vs target 20%)
        "enzyme_cat":                    2.0,
        # Boost TMP/SMX (0% vs target 30%)
        "mutation_folate_pathway":       3.0,
        # Boost Tet/efflux (0% vs targets 35-55%)
        "efflux_acrab_tolc":             2.0,
        "global_efflux_pump":            2.0,
        "protection_tet_m":              2.0,
    },
}

# ---------- Process config.rs ----------
with open(CONFIG, "r") as f:
    lines = f.readlines()

# Pattern matches emergence rate lines, using [\w.\-] for bacteria names with dots/hyphens
line_pattern = re.compile(
    r'^(\s*map\.insert\("bacteria_)([\w.\-]+)(_mechanism_)([\w.]+)(_emergence_rate"'
    r'\.to_string\(\),\s*)'
    r'([\d.eE\-+]+)'
    r'(\);\s*//\s*tier\s+)(\d+)(.*\n)$'
)

modified = 0
reverted = 0
boosted = 0
selective = 0
capped = 0
skipped_zero = 0

for i, line in enumerate(lines):
    m = line_pattern.match(line)
    if not m:
        continue

    prefix1 = m.group(1)   # "        map.insert(\"bacteria_"
    bact = m.group(2)      # bacterium name
    mid1 = m.group(3)      # "_mechanism_"
    mech = m.group(4)      # mechanism name
    mid2 = m.group(5)      # "_emergence_rate\".to_string(), "
    old_val_str = m.group(6)
    mid3 = m.group(7)      # "); // tier "
    old_tier_str = m.group(8)
    trailing = m.group(9)  # rest of line

    old_val = float(old_val_str)

    new_val = None

    # GROUP A: Revert to backup then multiply
    if bact in GROUP_A:
        multiplier = GROUP_A[bact]
        backup_key = (bact, mech)
        if backup_key not in backup_lookup:
            continue  # shouldn't happen
        backup_val, _backup_tier = backup_lookup[backup_key]

        if backup_val == 0.0:
            new_val = 0.0  # preserve disabled mechanism
        else:
            new_val = backup_val * multiplier
        reverted += 1

    # GROUP B: Multiply current values uniformly
    elif bact in GROUP_B:
        multiplier = GROUP_B[bact]
        if old_val == 0.0:
            new_val = 0.0  # preserve disabled
            skipped_zero += 1
        else:
            new_val = old_val * multiplier
        boosted += 1

    # GROUP C: Per-mechanism selective multipliers
    elif bact in GROUP_C:
        mech_rules = GROUP_C[bact]
        if mech in mech_rules:
            multiplier = mech_rules[mech]
            if old_val == 0.0:
                new_val = 0.0
                skipped_zero += 1
            else:
                new_val = old_val * multiplier
            selective += 1
        else:
            continue  # no change for this mechanism
    else:
        continue  # bacterium not in any group

    # Cap at MAX_RATE
    if new_val > MAX_RATE:
        new_val = MAX_RATE
        capped += 1

    # Compute new tier
    new_tier = value_to_tier(new_val)

    # Format the new value
    new_val_str = format_value(new_val)

    # Check if anything actually changed
    # Compare values within floating point tolerance
    if new_val == 0.0 and old_val == 0.0:
        # Still write it to ensure format consistency
        pass

    lines[i] = f"{prefix1}{bact}{mid1}{mech}{mid2}{new_val_str}{mid3}{new_tier}{trailing}"
    modified += 1

with open(CONFIG, "w") as f:
    f.writelines(lines)

print(f"Modified {modified} emergence rate lines total.")
print(f"  - GROUP A reverts:     {reverted}")
print(f"  - GROUP B boosts:      {boosted}")
print(f"  - GROUP C selective:   {selective}")
print(f"  - Capped at 1e-3:     {capped}")
print(f"  - Preserved zeros:     {skipped_zero}")
print("Done.")
