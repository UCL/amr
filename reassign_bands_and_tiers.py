#!/usr/bin/env python3
"""
Reassign bacteria emergence-rate bands (incidence-based) and apply first-pass
tier calibration using the K. pneumoniae anchor.

This script modifies src/config.rs in place.

Approach
--------
1. Bands are set mechanically: band_mult = REFERENCE_CONSTANT / target_incidence(%)
   where REFERENCE_CONSTANT is chosen so the K. pneumoniae anchor is preserved.
   K. pneumoniae (0.15% incidence) was the only bacterium showing genuine
   drug-class differentiation at stored rate ~2e-5 with ~400 PDs.

2. Tier calibration: For each bacterium, the hallmark (maximum non-zero)
   stored rate is rescaled so that:
       target_stored = 2e-5 × (mean_target_R / 30) × (400 / PDs)
   All other mechanisms for that bacterium are scaled by the same factor,
   preserving relative ordering.

3. The tier system documentation is extended to include tiers 10-12.
"""

import re
import math
import sys
import os
from collections import OrderedDict

# ═══════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════

# Reference: K. pneumoniae (0.15% incidence) works at band ×10
# So constant = 10 × 0.15 = 1.5
# band_mult = 1.5 / target_incidence(%)
REFERENCE_CONSTANT = 1.5

# Anchor formula: target_stored = ANCHOR_RATE × (mean_R / ANCHOR_R) × (ANCHOR_PDS / PDs)
ANCHOR_RATE = 2.0e-5   # K. pneumoniae hallmark stored rate
ANCHOR_R = 30.0         # Mean target R that produced this
ANCHOR_PDS = 400.0      # PDs at that working point

# Extended tier ladder (half-log spacing)
TIER_LADDER = {
    0:  0.0,
    1:  1.0e-09,
    2:  1.0e-08,
    3:  5.0e-08,
    4:  1.0e-07,
    5:  5.0e-07,
    6:  1.0e-06,
    7:  5.0e-06,
    8:  1.0e-05,
    9:  5.0e-05,
    10: 1.0e-04,
    11: 5.0e-04,
    12: 1.0e-03,
}

# ═══════════════════════════════════════════════════════════════════════════
# CALIBRATION DATA  (from calibration_summary_411919.txt / 058298.txt)
# ═══════════════════════════════════════════════════════════════════════════

# Mapping: calibration file name → config key prefix
NAME_TO_KEY = {
    "acinetobacter baumannii":                     "acinetobacter_baumannii",
    "bacteroides fragilis":                        "bacteroides_fragilis",
    "bordetella pertussis":                        "bordetella_pertussis",
    "burkholderia cepacia complex":                "burkholderia_cepacia_complex",
    "campylobacter jejuni":                        "campylobacter_jejuni",
    "chlamydia trachomatis":                       "chlamydia_trachomatis",
    "citrobacter spp.":                            "citrobacter_spp.",
    "clostridioides difficile":                    "clostridioides_difficile",
    "enterobacter cloacae":                        "enterobacter_cloacae",
    "enterobacter spp.":                           "enterobacter_spp.",
    "enterococcus faecalis":                       "enterococcus_faecalis",
    "enterococcus faecium":                        "enterococcus_faecium",
    "escherichia coli":                            "escherichia_coli",
    "haemophilus influenzae":                       "haemophilus_influenzae",
    "helicobacter pylori":                         "helicobacter_pylori",
    "invasive non-typhoidal salmonella spp.":      "invasive_non-typhoidal_salmonella_spp.",
    "klebsiella pneumoniae":                       "klebsiella_pneumoniae",
    "legionella pneumophila":                      "legionella_pneumophila",
    "listeria monocytogenes":                      "listeria_monocytogenes",
    "mdr mycobacterium tuberculosis":              "mdr_mycobacterium_tuberculosis",
    "moraxella catarrhalis":                       "moraxella_catarrhalis",
    "morganella spp.":                             "morganella_spp.",
    "mycoplasma genitalium":                       "mycoplasma_genitalium",
    "mycoplasma pneumoniae":                       "mycoplasma_pneumoniae",
    "neisseria gonorrhoeae":                       "neisseria_gonorrhoeae",
    "neisseria meningitidis":                      "neisseria_meningitidis",
    "proteus spp.":                                "proteus_spp.",
    "providencia stuartii":                        "p_stuartii",
    "pseudomonas aeruginosa":                      "pseudomonas_aeruginosa",
    "salmonella enterica serovar paratyphi a":     "salmonella_enterica_serovar_paratyphi_a",
    "salmonella enterica serovar typhi":           "salmonella_enterica_serovar_typhi",
    "serratia spp.":                               "serratia_spp.",
    "shigella spp.":                               "shigella_spp.",
    "staphylococcus aureus":                       "staphylococcus_aureus",
    "staphylococcus epidermidis":                  "staphylococcus_epidermidis",
    "stenotrophomonas maltophilia":                "stenotrophomonas_maltophilia",
    "streptococcus agalactiae":                    "streptococcus_agalactiae",
    "streptococcus pneumoniae":                    "streptococcus_pneumoniae",
    "streptococcus pyogenes":                      "streptococcus_pyogenes",
    "treponema pallidum":                          "treponema_pallidum",
    "vibrio cholerae":                             "vibrio_cholerae",
    "yersinia enterocolitica":                     "yersinia_enterocolitica",
}

# Current band multipliers (from config.rs)
CURRENT_BANDS = {
    "escherichia_coli":                          0.1,
    "klebsiella_pneumoniae":                     10.0,
    "citrobacter_spp.":                          100.0,
    "enterobacter_spp.":                         100.0,
    "enterobacter_cloacae":                      100.0,
    "morganella_spp.":                           100.0,
    "proteus_spp.":                              10.0,
    "serratia_spp.":                             100.0,
    "p_stuartii":                                10.0,
    "salmonella_enterica_serovar_typhi":         10.0,
    "salmonella_enterica_serovar_paratyphi_a":   10.0,
    "invasive_non-typhoidal_salmonella_spp.":    10.0,
    "shigella_spp.":                             1.0,
    "yersinia_enterocolitica":                   100.0,
    "pseudomonas_aeruginosa":                    10.0,
    "acinetobacter_baumannii":                   10.0,
    "stenotrophomonas_maltophilia":              10.0,
    "burkholderia_cepacia_complex":              100.0,
    "vibrio_cholerae":                           10.0,
    "campylobacter_jejuni":                      1.0,
    "helicobacter_pylori":                       1.0,
    "neisseria_gonorrhoeae":                     1.0,
    "neisseria_meningitidis":                    10.0,
    "moraxella_catarrhalis":                     10.0,
    "haemophilus_influenzae":                     10.0,
    "legionella_pneumophila":                    10.0,
    "staphylococcus_aureus":                     10.0,
    "staphylococcus_epidermidis":                10.0,
    "streptococcus_pneumoniae":                  0.1,
    "streptococcus_pyogenes":                    10.0,
    "streptococcus_agalactiae":                  10.0,
    "enterococcus_faecalis":                     10.0,
    "enterococcus_faecium":                      100.0,
    "listeria_monocytogenes":                    100.0,
    "clostridioides_difficile":                  100.0,
    "bacteroides_fragilis":                      100.0,
    "bordetella_pertussis":                      10.0,
    "mycoplasma_genitalium":                     10.0,
    "mycoplasma_pneumoniae":                     10.0,
    "chlamydia_trachomatis":                     1.0,
    "treponema_pallidum":                        10.0,
    "mdr_mycobacterium_tuberculosis":            10.0,
}


# ═══════════════════════════════════════════════════════════════════════════
# PARSE CALIBRATION FILE
# ═══════════════════════════════════════════════════════════════════════════

def parse_calibration(cal_path):
    """
    Parse calibration_summary_NNNNNN.txt to extract:
    - target_incidence (%) for each bacterium
    - PD counts (infected person-days)
    - mean target infection resistance (%)
    """
    with open(cal_path, "r", encoding="utf-8") as f:
        text = f.read()
    lines = text.split("\n")

    # ---- Parse Bacteria Burden Benchmarks ----
    target_incidences = {}
    in_burden = False
    for line in lines:
        if "Bacteria Burden Benchmarks" in line:
            in_burden = True
            continue
        if in_burden and line.strip() == "":
            continue
        if in_burden and ("Infection target" in line or "Bacteria" in line):
            continue  # header
        if in_burden and "Note:" in line:
            in_burden = False
            continue
        if in_burden:
            # Fixed-width fields: bacterium name (40 chars), then numeric fields
            stripped = line.strip()
            if not stripped:
                in_burden = False
                continue
            # Parse by finding the first number after the bacterium name
            parts = re.split(r'\s{2,}', stripped)
            if len(parts) >= 3:
                bact_name = parts[0].strip().lower()
                try:
                    target_inf = float(parts[1])
                    target_incidences[bact_name] = target_inf
                except ValueError:
                    pass

    # ---- Parse Resistance Benchmarks for PDs and target R ----
    pd_counts = {}        # config_key → int PDs
    target_r_lists = {}   # config_key → list of target R values
    in_resistance = False
    for line in lines:
        if "Resistance Benchmarks (percent resistant)" in line:
            in_resistance = True
            continue
        if in_resistance and ("Bacteria" in line and "Drug" in line and "Drug class" in line):
            continue  # header
        if not in_resistance:
            continue
        stripped = line.strip()
        if not stripped:
            continue

        # Parse the fixed-width resistance benchmark lines
        # The bacterium name is in the first ~40 chars
        parts = re.split(r'\s{2,}', stripped)
        if len(parts) < 8:
            continue

        bact_name = parts[0].strip().lower()
        if bact_name not in NAME_TO_KEY:
            continue

        config_key = NAME_TO_KEY[bact_name]

        # Find "Infection resistance target (%)" - field index 4
        # And "Infected person-days" - field index 8
        # Fields: Bacteria, Drug, Drug class, InfR_sim, InfR_target, AvgR_sim, AvgR_target,
        #         Microbiome_sim, Microbiome_target, Infected_PDs, Resistant_PDs, Microbiome_carrier_days, Note
        # But the layout can vary. Let me parse more carefully.
        # The line format is fixed-width. Let's use character positions.
        pass

    # Since fixed-width parsing is fragile, let me use a simpler approach:
    # search for PD values per bacterium by finding numeric PD values in each line.

    # Reset and re-parse with regex
    pd_counts = {}
    target_r_lists = {}
    in_resistance = False

    for line in lines:
        if "Resistance Benchmarks (percent resistant)" in line:
            in_resistance = True
            continue
        if not in_resistance:
            continue
        stripped = line.strip()
        if not stripped or stripped.startswith("Bacteria"):
            continue

        # Try to extract bacteria name (first field, left-aligned, up to 40 chars)
        # Then drug, drug class, and numeric fields
        # Use the original line (not stripped) to work with positions
        if len(line) < 100:
            continue

        # Extract bacterium name from the first ~41 chars (Drug column starts at col 41)
        bact_raw = line[:41].strip().lower()
        if not bact_raw or bact_raw.startswith("bacteria"):
            continue
        if bact_raw not in NAME_TO_KEY:
            # Try partial match
            found = False
            for cal_name in NAME_TO_KEY:
                if bact_raw.startswith(cal_name[:20]):
                    bact_raw = cal_name
                    found = True
                    break
            if not found:
                continue

        config_key = NAME_TO_KEY[bact_raw]

        # Extract Infection resistance target (%) - it's after Infection resistance simulation
        # Find all numbers in the line
        # Pattern: look for the target resistance and PDs
        # The columns are at roughly fixed positions.
        # Infection resistance simulation: chars ~100-140
        # Infection resistance target: chars ~140-175
        # Infected person-days: chars ~280-305

        # Let me extract all potential numbers from the line
        # Parse by splitting on 2+ spaces
        rest = line[42:]
        fields = re.split(r'\s{2,}', rest.strip())

        # We need to find:
        # - Infection resistance target (%)  → value or "---"
        # - Infected person-days             → integer or "---"

        # Fields after bacteria name and drug/class:
        # drug_name, drug_class, inf_r_sim, inf_r_target, avg_r_sim, avg_r_target,
        # micro_sim, micro_target, inf_pds, res_pds, micro_carrier, note

        if len(fields) < 6:
            continue

        # Find infection resistance target - it could be at various positions
        # Let me look for numeric values matching the known patterns
        inf_r_target = None
        pds = None

        # Parse the full line trying to identify field positions
        # The fixed-width format has these approximate positions:
        # Col 0-41: Bacteria name
        # Col 42-66: Drug name
        # Col 67-103: Drug class
        # Col 104-143: Infection resistance simulation (%)
        # Col 144-175: Infection resistance target (%)
        # Col 176-200: Average resistant simulation
        # Col 201-225: Average resistant target
        # Col 226-248: Microbiome simulation
        # Col 249-267: Microbiome target
        # Col 268-290: Infected person-days
        # Col 291-312: Resistant person-days
        # Col 313-340: Microbiome carrier-days
        # Col 341+: Note

        try:
            # infection resistance target (col ~140-176)
            target_str = line[140:176].strip() if len(line) > 176 else ""
            if target_str and target_str != "---":
                inf_r_target = float(target_str)

            # infected person-days (col ~271-293)
            pds_str = line[271:293].strip() if len(line) > 293 else ""
            if pds_str:
                pds_str = pds_str.replace(",", "")
                if pds_str != "---" and pds_str:
                    pds_val = int(pds_str)
                    if config_key not in pd_counts or pds_val > 0:
                        pd_counts[config_key] = pds_val
        except (ValueError, IndexError):
            pass

        if inf_r_target is not None:
            if config_key not in target_r_lists:
                target_r_lists[config_key] = []
            target_r_lists[config_key].append(inf_r_target)

    # Compute mean target R
    mean_target_r = {}
    for key, r_list in target_r_lists.items():
        if r_list:
            mean_target_r[key] = sum(r_list) / len(r_list)

    return target_incidences, pd_counts, mean_target_r


# ═══════════════════════════════════════════════════════════════════════════
# PARSE CONFIG.RS
# ═══════════════════════════════════════════════════════════════════════════

def find_nearest_tier(value, band_mult, tier_ladder=TIER_LADDER):
    """Given a stored rate and band multiplier, find the nearest tier."""
    if value == 0.0:
        return 0
    tier_base = value / band_mult
    if tier_base <= 0:
        return 0

    best_tier = 0
    best_ratio = float('inf')
    for t, base in tier_ladder.items():
        if base == 0.0:
            continue
        ratio = abs(math.log10(tier_base) - math.log10(base))
        if ratio < best_ratio:
            best_ratio = ratio
            best_tier = t
    return best_tier


def format_rate(value):
    """Format emergence rate in scientific notation matching Rust config style."""
    if value == 0.0:
        return "0.0"
    # Use e notation: X.Ye-NN
    exp = math.floor(math.log10(abs(value)))
    mantissa = value / (10 ** exp)

    # Snap to clean values: round mantissa to 1 decimal
    mantissa = round(mantissa, 1)
    if mantissa >= 10.0:
        mantissa /= 10.0
        exp += 1
    if mantissa < 1.0:
        mantissa *= 10.0
        exp -= 1

    if mantissa == int(mantissa):
        mant_str = f"{int(mantissa)}.0"
    else:
        mant_str = f"{mantissa}"

    if exp >= 0:
        return f"{mant_str}e+{exp:02d}"
    else:
        return f"{mant_str}e-{abs(exp):02d}"


def snap_to_tier_value(value, band_mult):
    """
    Snap a stored rate to the nearest tier×band value on the tier ladder.
    Returns (snapped_stored, tier_number).
    """
    if value == 0.0:
        return 0.0, 0

    tier_base = value / band_mult
    if tier_base <= 0:
        return 0.0, 0

    best_tier = 0
    best_base = 0.0
    best_dist = float('inf')
    for t, base in TIER_LADDER.items():
        if base == 0.0:
            continue
        dist = abs(math.log10(tier_base) - math.log10(base))
        if dist < best_dist:
            best_dist = dist
            best_tier = t
            best_base = base

    return best_base * band_mult, best_tier


# ═══════════════════════════════════════════════════════════════════════════
# MAIN TRANSFORMATION
# ═══════════════════════════════════════════════════════════════════════════

def main():
    config_path = os.path.join(os.path.dirname(__file__), "src", "config.rs")
    cal_path = os.path.join(os.path.dirname(__file__), "output_graphs",
                            "calibration_summary_411919.txt")

    if not os.path.exists(cal_path):
        # Fallback: try 058298
        cal_path = os.path.join(os.path.dirname(__file__), "output_graphs",
                                "calibration_summary_058298.txt")

    print(f"Parsing calibration: {cal_path}")
    target_incidences, pd_counts, mean_target_r = parse_calibration(cal_path)

    print(f"\nParsed {len(target_incidences)} target incidences")
    print(f"Parsed {len(pd_counts)} PD counts")
    print(f"Parsed {len(mean_target_r)} mean target R values")

    # Fill in any missing data with reasonable defaults
    for cal_name, config_key in NAME_TO_KEY.items():
        if config_key not in pd_counts:
            # Estimate PDs from target incidence: PDs ≈ pop × incidence × 7 days
            if cal_name in target_incidences:
                est_pds = int(63875 * target_incidences[cal_name] / 100.0 * 7)
                pd_counts[config_key] = max(est_pds, 50)
                print(f"  Estimated PDs for {config_key}: {pd_counts[config_key]}")
            else:
                pd_counts[config_key] = 200  # fallback
                print(f"  Fallback PDs for {config_key}: 200")
        if config_key not in mean_target_r:
            mean_target_r[config_key] = 25.0  # conservative default
            print(f"  Default mean_target_R for {config_key}: 25.0%")

    # ---- Compute new band multipliers ----
    new_bands = {}
    for cal_name, config_key in NAME_TO_KEY.items():
        if cal_name in target_incidences:
            incidence = target_incidences[cal_name]
            new_bands[config_key] = REFERENCE_CONSTANT / incidence
        else:
            # Keep current band
            new_bands[config_key] = CURRENT_BANDS.get(config_key, 10.0)
            print(f"  Keeping current band for {config_key}: {new_bands[config_key]}")

    # ---- Read config.rs ----
    with open(config_path, "r", encoding="utf-8") as f:
        content = f.read()
    lines = content.split("\n")

    # ---- Process each bacterium block ----
    # Pattern to match emergence rate lines
    rate_pattern = re.compile(
        r'(\s*map\.insert\("bacteria_)([a-zA-Z0-9_\-.]+)'
        r'(_mechanism_)([a-zA-Z0-9_]+)(_emergence_rate"\.to_string\(\),\s*)'
        r'([0-9eE.+\-]+)'
        r'(\);\s*//\s*tier\s*)(\d+)(.*)'
    )
    band_pattern = re.compile(r'^(\s*// Band )(\d+)( \(x)([0-9.eE+\-]+)(\)\s*)$')

    # Collect per-bacterium data for tier calibration
    bacterium_lines = {}  # config_key → [(line_idx, mechanism, stored_rate)]

    for i, line in enumerate(lines):
        m = rate_pattern.match(line)
        if m:
            config_key = m.group(2)
            mechanism = m.group(4)
            stored_rate = float(m.group(6))
            if config_key not in bacterium_lines:
                bacterium_lines[config_key] = []
            bacterium_lines[config_key].append((i, mechanism, stored_rate))

    print(f"\nFound {len(bacterium_lines)} bacteria blocks in config.rs")

    # ---- Compute scale factors and apply ----
    results_table = []

    for config_key, mech_list in bacterium_lines.items():
        old_band = CURRENT_BANDS.get(config_key, 10.0)
        new_band = new_bands.get(config_key, old_band)
        pds = pd_counts.get(config_key, 400)
        mean_r = mean_target_r.get(config_key, 25.0)

        # Find hallmark (max non-zero) stored rate
        non_zero_rates = [r for _, _, r in mech_list if r > 0]
        if not non_zero_rates:
            # All zero — skip this bacterium
            results_table.append({
                'key': config_key, 'old_band': old_band, 'new_band': new_band,
                'pds': pds, 'mean_r': mean_r, 'scale': 1.0,
                'hallmark_old': 0.0, 'hallmark_new': 0.0,
                'note': 'all-zero (no change)'
            })
            continue

        hallmark_stored = max(non_zero_rates)

        # Compute target hallmark stored rate
        if mean_r > 0 and pds > 0:
            target_stored = ANCHOR_RATE * (mean_r / ANCHOR_R) * (ANCHOR_PDS / pds)
        else:
            target_stored = hallmark_stored  # no change

        # Scale factor: applied to ALL mechanisms
        scale_factor = target_stored / hallmark_stored

        # Clamp scale factor to avoid extreme changes
        # (limit to 100× up or 0.001× down for safety)
        scale_factor = max(1e-3, min(100.0, scale_factor))

        results_table.append({
            'key': config_key, 'old_band': old_band, 'new_band': new_band,
            'pds': pds, 'mean_r': mean_r, 'scale': scale_factor,
            'hallmark_old': hallmark_stored, 'hallmark_new': hallmark_stored * scale_factor,
            'note': ''
        })

        # Apply to each mechanism line
        for line_idx, mechanism, old_rate in mech_list:
            if old_rate == 0.0:
                continue  # Preserve tier 0

            new_rate = old_rate * scale_factor

            # Snap to nearest tier×band value
            snapped_rate, new_tier = snap_to_tier_value(new_rate, new_band)

            # Format the new rate
            rate_str = format_rate(snapped_rate)

            # Reconstruct the line
            m = rate_pattern.match(lines[line_idx])
            if m:
                lines[line_idx] = (
                    f"{m.group(1)}{m.group(2)}{m.group(3)}{m.group(4)}"
                    f"{m.group(5)}{rate_str}); // tier {new_tier}"
                )

    # ---- Update band comments ----
    for i, line in enumerate(lines):
        bm = band_pattern.match(line)
        if bm:
            # Find which bacterium this band comment belongs to
            # Look at the next few lines to find a map.insert with the config key
            for j in range(i + 1, min(i + 3, len(lines))):
                rm = rate_pattern.match(lines[j])
                if rm:
                    config_key = rm.group(2)
                    new_band = new_bands.get(config_key)
                    if new_band is not None:
                        # Format band multiplier nicely
                        if new_band >= 1.0:
                            band_str = f"x{new_band:.1f}" if new_band != int(new_band) else f"x{int(new_band)}"
                        else:
                            band_str = f"x{new_band:.4f}".rstrip('0').rstrip('.')

                        # Find the nearest "band number" for documentation
                        band_num = find_band_number(new_band)
                        lines[i] = f"{bm.group(1)}{band_num} ({band_str})"
                    break

    # ---- Update tier system documentation ----
    # Find the tier documentation block and extend it
    for i, line in enumerate(lines):
        if "tier 9 = 5.0e-05" in line and "endemic" in line:
            # Add tiers 10-12 after tier 9
            indent = "        //   "
            lines.insert(i + 1, f"{indent}tier 10 = 1.0e-04   — extreme / very high frequency")
            lines.insert(i + 2, f"{indent}tier 11 = 5.0e-04   — extraordinary / maximum plausible")
            lines.insert(i + 3, f"{indent}tier 12 = 1.0e-03   — theoretical maximum")
            break

    # ---- Update band documentation ----
    for i, line in enumerate(lines):
        if "INFECTION INCIDENCE BANDS" in line:
            # Find the band documentation block and replace it
            j = i + 1
            new_band_docs = [
                "        //   Band multiplier = 1.5 / target_infection_incidence(%)",
                "        //   Reference: K. pneumoniae (0.15%) → ×10",
                "        //   Higher incidence → lower multiplier (normalizes for PD exposure)",
                "        //   Example multipliers:",
                "        //     E. coli (2.5%)        → ×0.6",
                "        //     S. pneumoniae (2.0%)   → ×0.75",
                "        //     N. gonorrhoeae (1.0%)   → ×1.5",
                "        //     K. pneumoniae (0.15%)  → ×10",
                "        //     A. baumannii (0.015%)  → ×100",
            ]
            # Find end of old band documentation (look for "Final rate" line)
            while j < len(lines) and "Final rate" not in lines[j]:
                j += 1
            # Replace old band docs (from i+1 to j-1) with new docs
            lines[i + 1:j] = new_band_docs
            break

    # ---- Write back ----
    new_content = "\n".join(lines)
    with open(config_path, "w", encoding="utf-8") as f:
        f.write(new_content)

    # ---- Print verification table ----
    print("\n" + "=" * 120)
    print("VERIFICATION TABLE")
    print("=" * 120)
    print(f"{'Bacterium':<45} {'OldBand':>8} {'NewBand':>10} {'PDs':>6} "
          f"{'MeanR%':>7} {'Scale':>8} {'HallOld':>10} {'HallNew':>10} {'Note'}")
    print("-" * 120)

    results_table.sort(key=lambda x: x['key'])
    for r in results_table:
        print(f"{r['key']:<45} {r['old_band']:>8.2f} {r['new_band']:>10.4f} "
              f"{r['pds']:>6} {r['mean_r']:>7.1f} {r['scale']:>8.4f} "
              f"{r['hallmark_old']:>10.2e} {r['hallmark_new']:>10.2e} {r['note']}")

    print(f"\nDone! Modified {config_path}")
    print(f"Total bacteria processed: {len(results_table)}")


def find_band_number(mult):
    """Find the nearest 'band number' for documentation purposes."""
    if mult <= 0:
        return 0
    # Bands are roughly log-spaced. Use log10(mult) + 6 as band number
    # But since bands are now incidence-specific, just show the multiplier
    # Map to rough band numbers for consistency:
    log_val = math.log10(mult)
    # Band 6 = ×1 (log=0), each band is ~1 decade
    band = int(round(log_val + 6))
    band = max(0, min(10, band))
    return band


if __name__ == "__main__":
    main()
