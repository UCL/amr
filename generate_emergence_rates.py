#!/usr/bin/env python3
"""Generate biologically-tiered emergence rates for AMR simulation config.rs

Tier system (half-log spacing):
  0 = 0.0       (does not emerge)
  1 = 1e-9      (near-impossible)
  2 = 1e-8      (exceedingly rare)
  3 = 5e-8      (very rare)
  4 = 1e-7      (rare)
  5 = 5e-7      (uncommon but documented)
  6 = 1e-6      (moderate / clinically significant)
  7 = 5e-6      (common / well-established)
  8 = 1e-5      (very common / hallmark)
  9 = 5e-5      (endemic / highest-frequency)

Incidence bands (multiplied into final rate):
  High     (>0.005/yr)   x0.1
  Moderate (0.001-0.005) x1.0
  Low      (0.0002-0.001) x3.0
  VeryLow  (<0.0002)     x10.0
"""

import os
import re

TIER_VALUES = {
    0: 0.0, 1: 1e-9, 2: 1e-8, 3: 5e-8, 4: 1e-7,
    5: 5e-7, 6: 1e-6, 7: 5e-6, 8: 1e-5, 9: 5e-5,
}

BAND_MULT = {"High": 0.1, "Moderate": 1.0, "Low": 3.0, "VeryLow": 10.0}
BAND_DISPLAY = {"High": "High", "Moderate": "Moderate", "Low": "Low", "VeryLow": "Very Low"}
BAND_FACTOR = {"High": "0.1", "Moderate": "1.0", "Low": "3.0", "VeryLow": "10.0"}

MECHANISMS = [
    "enzyme_esbl_ctx_m", "enzyme_esbl_tem", "enzyme_esbl_shv",
    "enzyme_ampc_cmy", "enzyme_ampc_dha",
    "enzyme_kpc", "enzyme_ndm_vim", "enzyme_oxa_48",
    "enzyme_cat", "enzyme_16s_rrmt",
    "target_site_pbp2a_meca", "target_site_van_a", "target_site_van_b",
    "target_site_erm_b", "target_site_cfr",
    "mutation_gyra_primary", "mutation_gyra_parc_secondary",
    "protection_qnr",
    "efflux_acrab_tolc", "efflux_mexxy_oprm", "global_efflux_pump",
    "porin_loss_ompk35_36", "porin_loss_oprd", "global_porin_loss",
    "modification_mcr_1",
    "mutation_folate_pathway", "mutation_nitroreductase", "enzyme_fos_a",
    "mutation_mpr_f", "mutation_rpo_b", "protection_fus_b",
    "protection_tet_m",
    "as_yet_unknown_1", "as_yet_unknown_2", "as_yet_unknown_3",
]

assert len(MECHANISMS) == 35

# =====================================================================
# ORGANISM DATA
# Each entry: (config_key, display_name, description, band, section_header_or_None, [35 tiers])
# Section headers are printed once when present, before the organism block.
# =====================================================================

ORGS = [
    # ---- Gram-Negative Enterobacterales ----
    ("escherichia_coli", "E. coli", "Gram-negative, Enterobacterales", "High",
     "Gram-Negative Bacteria — Enterobacterales",
     #ctx tem shv cmy dha kpc ndm oxa cat 16s  p2a vA  vB  erm cfr  gyA gyP qnr acr mxy efx  omk opd gpn mcr fol nit fos  mpr rpo fus tet  u1  u2  u3
     [8,  7,  5,  6,  4,  3,  4,  3,  5,  2,   0,  0,  0,  0,  0,   7,  6,  5,  7,  1,  5,   5,  2,  4,  4,  7,  4,  4,   0,  2,  0,  7,   3,  3,  3]),

    ("klebsiella_pneumoniae", "K. pneumoniae", "Gram-negative, Enterobacterales", "Moderate",
     None,
     [7,  5,  7,  4,  4,  6,  5,  5,  4,  3,   0,  0,  0,  0,  0,   6,  5,  5,  5,  1,  4,   6,  2,  5,  4,  6,  3,  5,   0,  2,  0,  5,   3,  3,  3]),

    ("citrobacter_spp.", "Citrobacter spp.", "Gram-negative, Enterobacterales", "VeryLow",
     None,
     [5,  5,  4,  6,  5,  3,  3,  3,  5,  2,   0,  0,  0,  0,  0,   6,  5,  5,  6,  1,  4,   4,  2,  4,  4,  6,  4,  4,   0,  2,  0,  6,   3,  3,  3]),

    ("enterobacter_spp.", "Enterobacter spp.", "Gram-negative, Enterobacterales", "VeryLow",
     None,
     [6,  5,  5,  7,  6,  5,  4,  3,  5,  2,   0,  0,  0,  0,  0,   6,  5,  5,  6,  1,  4,   5,  2,  4,  4,  6,  4,  4,   0,  2,  0,  6,   3,  3,  3]),

    ("enterobacter_cloacae", "E. cloacae", "Gram-negative, Enterobacterales", "VeryLow",
     None,
     [6,  5,  5,  7,  6,  5,  4,  3,  5,  2,   0,  0,  0,  0,  0,   6,  5,  5,  6,  1,  4,   5,  2,  4,  4,  6,  4,  4,   0,  2,  0,  6,   3,  3,  3]),

    ("morganella_spp.", "Morganella spp.", "Gram-negative, Enterobacterales", "VeryLow",
     None,
     [5,  5,  4,  6,  4,  3,  3,  3,  5,  2,   0,  0,  0,  0,  0,   6,  5,  4,  5,  1,  4,   4,  2,  4,  3,  6,  4,  4,   0,  2,  0,  6,   3,  3,  3]),

    ("proteus_spp.", "Proteus spp.", "Gram-negative, Enterobacterales", "Low",
     None,
     [5,  6,  4,  5,  4,  3,  3,  3,  5,  2,   0,  0,  0,  0,  0,   6,  5,  5,  6,  1,  4,   4,  2,  4,  4,  6,  1,  4,   0,  2,  0,  5,   3,  3,  3]),

    ("serratia_spp.", "Serratia spp.", "Gram-negative, Enterobacterales", "VeryLow",
     None,
     [5,  5,  4,  6,  5,  4,  3,  3,  5,  2,   0,  0,  0,  0,  0,   6,  5,  5,  6,  1,  4,   4,  2,  4,  3,  6,  4,  4,   0,  2,  0,  6,   3,  3,  3]),

    ("p_stuartii", "P. stuartii", "Gram-negative, Enterobacterales", "VeryLow",
     None,
     [5,  5,  4,  6,  5,  3,  3,  3,  5,  2,   0,  0,  0,  0,  0,   6,  5,  5,  5,  1,  4,   4,  2,  4,  3,  6,  4,  4,   0,  2,  0,  6,   3,  3,  3]),

    ("salmonella_enterica_serovar_typhi", "S. Typhi", "Gram-negative, Enterobacterales", "Moderate",
     None,
     [5,  5,  3,  2,  2,  1,  1,  1,  6,  1,   0,  0,  0,  0,  0,   7,  6,  4,  5,  1,  4,   3,  1,  3,  3,  7,  3,  3,   0,  2,  0,  6,   3,  3,  3]),

    ("salmonella_enterica_serovar_paratyphi_a", "S. Paratyphi A", "Gram-negative, Enterobacterales", "Low",
     None,
     [4,  4,  3,  2,  2,  1,  1,  1,  5,  1,   0,  0,  0,  0,  0,   6,  5,  4,  4,  1,  3,   3,  1,  3,  3,  6,  3,  3,   0,  2,  0,  5,   3,  3,  3]),

    ("invasive_non-typhoidal_salmonella_spp.", "iNTS", "Gram-negative, Enterobacterales", "Low",
     None,
     [6,  5,  4,  3,  3,  2,  2,  2,  6,  2,   0,  0,  0,  0,  0,   5,  4,  4,  5,  1,  4,   3,  1,  3,  3,  6,  3,  3,   0,  2,  0,  6,   3,  3,  3]),

    ("shigella_spp.", "Shigella spp.", "Gram-negative, Enterobacterales", "High",
     None,
     [6,  6,  4,  4,  3,  2,  2,  2,  5,  2,   0,  0,  0,  0,  0,   7,  6,  5,  6,  1,  4,   4,  1,  3,  3,  7,  3,  3,   0,  2,  0,  7,   3,  3,  3]),

    ("yersinia_enterocolitica", "Y. enterocolitica", "Gram-negative, Enterobacterales", "VeryLow",
     None,
     [4,  4,  3,  6,  5,  2,  2,  2,  5,  2,   0,  0,  0,  0,  0,   5,  4,  4,  5,  1,  4,   4,  2,  4,  3,  5,  4,  3,   0,  2,  0,  5,   3,  3,  3]),

    # ---- Non-fermenting Gram-Negatives ----
    ("pseudomonas_aeruginosa", "P. aeruginosa", "Gram-negative, NonFermenter", "Low",
     "Non-fermenting Gram-Negatives",
     [2,  2,  2,  7,  5,  4,  5,  3,  2,  3,   0,  0,  0,  0,  0,   7,  6,  3,  4,  7,  7,   2,  7,  6,  3,  3,  0,  4,   0,  2,  0,  3,   3,  3,  3]),

    ("acinetobacter_baumannii", "A. baumannii", "Gram-negative, NonFermenter", "VeryLow",
     None,
     [3,  3,  2,  6,  4,  3,  5,  7,  3,  3,   0,  0,  0,  0,  0,   6,  5,  2,  5,  4,  6,   3,  4,  5,  3,  4,  0,  2,   0,  3,  0,  5,   3,  3,  3]),

    ("stenotrophomonas_maltophilia", "S. maltophilia", "Gram-negative, NonFermenter", "VeryLow",
     None,
     [1,  1,  1,  3,  2,  1,  2,  1,  3,  2,   0,  0,  0,  0,  0,   5,  4,  2,  4,  3,  6,   3,  3,  4,  1,  5,  0,  2,   0,  2,  0,  4,   3,  3,  3]),

    ("burkholderia_cepacia_complex", "B. cepacia complex", "Gram-negative, NonFermenter", "VeryLow",
     None,
     [2,  2,  2,  5,  4,  2,  3,  2,  3,  2,   0,  0,  0,  0,  0,   5,  4,  2,  4,  4,  6,   3,  4,  5,  2,  4,  0,  3,   0,  2,  0,  3,   3,  3,  3]),

    # ---- Other Gram-Negatives ----
    ("vibrio_cholerae", "V. cholerae", "Gram-negative, EntericPathogen", "Low",
     "Other Gram-Negatives",
     [4,  5,  2,  2,  2,  1,  2,  1,  6,  1,   0,  0,  0,  0,  0,   6,  5,  4,  4,  1,  4,   2,  1,  3,  2,  7,  2,  2,   0,  2,  0,  7,   3,  3,  3]),

    ("campylobacter_jejuni", "C. jejuni", "Gram-negative, Helicobacter group", "High",
     None,
     [0,  0,  0,  0,  0,  0,  0,  0,  4,  0,   0,  0,  0,  0,  0,   7,  5,  0,  0,  0,  5,   0,  0,  3,  0,  3,  0,  0,   0,  2,  0,  6,   3,  3,  3]),

    ("helicobacter_pylori", "H. pylori", "Gram-negative, Helicobacter group", "High",
     None,
     [0,  0,  0,  0,  0,  0,  0,  0,  3,  0,   0,  0,  0,  0,  0,   6,  5,  0,  0,  0,  4,   0,  0,  2,  0,  2,  0,  0,   0,  4,  0,  4,   3,  3,  3]),

    ("neisseria_gonorrhoeae", "N. gonorrhoeae", "Gram-negative, Fastidious", "High",
     None,
     [2,  6,  1,  2,  1,  0,  1,  0,  4,  1,   0,  0,  0,  5,  1,   7,  6,  2,  3,  1,  5,   2,  1,  4,  1,  4,  2,  0,   0,  2,  0,  6,   3,  3,  3]),

    ("neisseria_meningitidis", "N. meningitidis", "Gram-negative, Fastidious", "VeryLow",
     None,
     [1,  3,  1,  1,  1,  0,  0,  0,  3,  1,   0,  0,  0,  3,  1,   4,  3,  1,  2,  1,  3,   1,  1,  2,  1,  4,  2,  0,   0,  3,  0,  4,   2,  2,  2]),

    ("moraxella_catarrhalis", "M. catarrhalis", "Gram-negative, Fastidious", "Low",
     None,
     [2,  5,  1,  3,  2,  0,  0,  0,  3,  1,   0,  0,  0,  4,  1,   4,  3,  2,  3,  1,  4,   1,  1,  2,  1,  5,  2,  0,   0,  2,  0,  5,   3,  3,  3]),

    ("haemophilus_influenzae", "H. influenzae", "Gram-negative, Fastidious", "Low",
     None,
     [2,  6,  2,  3,  2,  0,  0,  0,  4,  1,   0,  0,  0,  4,  1,   4,  3,  2,  3,  1,  4,   2,  1,  3,  1,  6,  2,  0,   0,  2,  0,  5,   3,  3,  3]),

    ("legionella_pneumophila", "L. pneumophila", "Gram-negative, Fastidious", "VeryLow",
     None,
     [1,  1,  1,  1,  1,  0,  0,  0,  2,  1,   0,  0,  0,  3,  1,   4,  3,  1,  1,  1,  3,   1,  1,  2,  1,  2,  1,  0,   0,  2,  0,  3,   2,  2,  2]),

    # ---- Gram-Positive Bacteria ----
    ("staphylococcus_aureus", "S. aureus", "Gram-positive, Staphylococcus", "High",
     "Gram-Positive Bacteria — Staphylococci",
     [0,  0,  0,  0,  0,  0,  0,  0,  4,  0,   7,  2,  1,  7,  2,   5,  4,  0,  0,  0,  5,   0,  0,  0,  0,  6,  0,  0,   4,  4,  5,  6,   3,  3,  3]),

    ("staphylococcus_epidermidis", "S. epidermidis", "Gram-positive, Staphylococcus", "Low",
     None,
     [0,  0,  0,  0,  0,  0,  0,  0,  4,  0,   7,  1,  1,  6,  2,   5,  4,  0,  0,  0,  4,   0,  0,  0,  0,  5,  0,  0,   3,  4,  4,  5,   3,  3,  3]),

    ("streptococcus_pneumoniae", "S. pneumoniae", "Gram-positive, Streptococcus", "High",
     "Streptococci",
     [0,  0,  0,  0,  0,  0,  0,  0,  4,  0,   0,  1,  1,  7,  1,   5,  5,  0,  0,  0,  5,   0,  0,  0,  0,  7,  0,  0,   1,  2,  1,  6,   3,  3,  3]),

    ("streptococcus_pyogenes", "S. pyogenes", "Gram-positive, Streptococcus", "Moderate",
     None,
     [0,  0,  0,  0,  0,  0,  0,  0,  3,  0,   0,  0,  0,  6,  1,   3,  2,  0,  0,  0,  4,   0,  0,  0,  0,  5,  0,  0,   1,  1,  1,  6,   3,  3,  3]),

    ("streptococcus_agalactiae", "S. agalactiae", "Gram-positive, Streptococcus", "Low",
     None,
     [0,  0,  0,  0,  0,  0,  0,  0,  3,  0,   1,  1,  1,  6,  1,   4,  3,  0,  0,  0,  4,   0,  0,  0,  0,  5,  0,  0,   2,  2,  2,  7,   3,  3,  3]),

    # ---- Enterococci ----
    ("enterococcus_faecalis", "E. faecalis", "Gram-positive, Enterococcus", "Low",
     "Enterococci",
     [0,  0,  0,  0,  0,  0,  0,  0,  4,  0,   1,  4,  4,  6,  2,   5,  4,  0,  0,  0,  4,   0,  0,  0,  0,  4,  0,  0,   3,  3,  2,  7,   3,  3,  3]),

    ("enterococcus_faecium", "E. faecium", "Gram-positive, Enterococcus", "VeryLow",
     None,
     [0,  0,  0,  0,  0,  0,  0,  0,  4,  0,   1,  7,  6,  6,  2,   6,  5,  0,  0,  0,  4,   0,  0,  0,  0,  4,  0,  0,   4,  3,  2,  7,   3,  3,  3]),

    # ---- Other Gram-Positives / Anaerobes ----
    ("listeria_monocytogenes", "L. monocytogenes", "Gram-positive, Listeria", "VeryLow",
     "Other Gram-Positives and Anaerobes",
     [0,  0,  0,  0,  0,  0,  0,  0,  3,  0,   0,  1,  1,  4,  1,   3,  2,  0,  0,  0,  3,   0,  0,  0,  0,  4,  0,  0,   2,  2,  2,  5,   2,  2,  2]),

    ("clostridioides_difficile", "C. difficile", "Anaerobe", "Low",
     None,
     [1,  2,  1,  2,  1,  0,  1,  0,  4,  1,   0,  0,  0,  5,  2,   5,  4,  2,  3,  1,  4,   2,  1,  3,  1,  3,  6,  0,   0,  5,  0,  6,   3,  3,  3]),

    ("bacteroides_fragilis", "B. fragilis", "Anaerobe", "Low",
     None,
     [2,  3,  2,  5,  4,  2,  4,  2,  4,  2,   0,  0,  0,  5,  2,   4,  3,  2,  3,  1,  4,   3,  1,  3,  1,  3,  5,  0,   0,  3,  0,  6,   3,  3,  3]),

    # ---- Fastidious / Atypical ----
    ("bordetella_pertussis", "B. pertussis", "Gram-negative, Fastidious", "Moderate",
     "Fastidious and Atypical Bacteria",
     [1,  2,  1,  1,  1,  0,  0,  0,  3,  1,   0,  0,  0,  5,  1,   3,  2,  1,  2,  1,  3,   1,  1,  2,  0,  4,  1,  0,   0,  2,  0,  4,   2,  2,  2]),

    ("mycoplasma_genitalium", "M. genitalium", "Atypical (no cell wall), Fastidious", "Moderate",
     None,
     [0,  0,  0,  0,  0,  0,  0,  0,  1,  0,   0,  0,  0,  7,  1,   6,  6,  0,  0,  0,  3,   0,  0,  1,  0,  1,  1,  0,   0,  2,  0,  4,   3,  3,  3]),

    ("mycoplasma_pneumoniae", "M. pneumoniae", "Atypical (no cell wall), Fastidious", "High",
     None,
     [0,  0,  0,  0,  0,  0,  0,  0,  1,  0,   0,  0,  0,  7,  1,   4,  3,  0,  0,  0,  3,   0,  0,  1,  0,  1,  1,  0,   0,  2,  0,  4,   2,  2,  2]),

    # ---- Obligate Intracellular / Special ----
    ("chlamydia_trachomatis", "C. trachomatis", "Obligate intracellular, Fastidious", "High",
     "Obligate Intracellular and Special Cases",
     [0,  0,  0,  0,  0,  0,  0,  0,  1,  0,   0,  0,  0,  3,  1,   4,  3,  0,  0,  0,  2,   0,  0,  1,  0,  1,  1,  0,   0,  2,  0,  4,   2,  2,  2]),

    ("treponema_pallidum", "T. pallidum", "Spirochete", "Moderate",
     None,
     [0,  0,  0,  0,  0,  0,  0,  0,  2,  0,   0,  0,  0,  0,  0,   4,  3,  0,  0,  0,  2,   0,  0,  2,  0,  2,  0,  0,   0,  2,  0,  4,   2,  2,  2]),

    # ---- Acid-Fast ----
    ("mdr_mycobacterium_tuberculosis", "MDR M. tuberculosis", "Acid-fast, Mycobacteria", "VeryLow",
     "Acid-Fast Bacteria",
     [0,  0,  0,  0,  0,  0,  0,  0,  3,  0,   0,  0,  0,  0,  0,   6,  5,  0,  0,  0,  5,   0,  0,  3,  0,  3,  0,  0,   0,  7,  0,  2,   3,  3,  3]),
]


def fmt_rate(rate):
    """Format rate for Rust source."""
    if rate == 0.0:
        return "0.0"
    return f"{rate:.1e}"


def generate_section():
    """Generate the full Rust code section."""
    lines = []
    indent = "        "  # 8 spaces

    # Header / legend
    lines.append(f"{indent}// --- Resistance Mechanisms Parameters ---")
    lines.append(f"{indent}// Implementation of granular resistance mechanisms (35 types)")
    lines.append(f"{indent}//")
    lines.append(f"{indent}// =====================================================================================")
    lines.append(f"{indent}// BACTERIA-MECHANISM-SPECIFIC EMERGENCE RATES")
    lines.append(f"{indent}// =====================================================================================")
    lines.append(f"{indent}// Direct emergence rates (per day when drug present) for each bacteria-mechanism pair.")
    lines.append(f"{indent}// All 35 resistance mechanisms in standardized biological order for all 42 bacteria.")
    lines.append(f"{indent}//")
    lines.append(f"{indent}// TIER SYSTEM (biological plausibility tiers, half-log spacing):")
    lines.append(f"{indent}//   tier 0 = 0.0       — does not emerge in this organism")
    lines.append(f"{indent}//   tier 1 = 1.0e-09   — near-impossible, placeholder only")
    lines.append(f"{indent}//   tier 2 = 1.0e-08   — exceedingly rare")
    lines.append(f"{indent}//   tier 3 = 5.0e-08   — very rare")
    lines.append(f"{indent}//   tier 4 = 1.0e-07   — rare / low clinical relevance")
    lines.append(f"{indent}//   tier 5 = 5.0e-07   — uncommon but documented")
    lines.append(f"{indent}//   tier 6 = 1.0e-06   — moderate / clinically significant")
    lines.append(f"{indent}//   tier 7 = 5.0e-06   — common / well-established")
    lines.append(f"{indent}//   tier 8 = 1.0e-05   — very common / hallmark mechanism")
    lines.append(f"{indent}//   tier 9 = 5.0e-05   — endemic / highest-frequency resistance")
    lines.append(f"{indent}//")
    lines.append(f"{indent}// INFECTION INCIDENCE BANDS (multiplier applied to tier base rate):")
    lines.append(f"{indent}//   High      (annual incidence > 0.005)   — x0.1  (many Bernoulli trials -> lower per-trial rate)")
    lines.append(f"{indent}//   Moderate  (0.001 - 0.005)              — x1.0  (reference)")
    lines.append(f"{indent}//   Low       (0.0002 - 0.001)             — x3.0  (fewer trials -> higher per-trial rate)")
    lines.append(f"{indent}//   Very Low  (< 0.0002)                   — x10.0 (rare infections -> highest per-trial rate)")
    lines.append(f"{indent}//")
    lines.append(f"{indent}// Final rate = tier_base_rate x band_multiplier")
    lines.append(f"{indent}// =====================================================================================")
    lines.append("")

    for org in ORGS:
        key, name, desc, band, section, tiers = org
        assert len(tiers) == 35, f"{key}: expected 35 tiers, got {len(tiers)}"
        mult = BAND_MULT[band]
        band_disp = BAND_DISPLAY[band]
        band_fac = BAND_FACTOR[band]

        # Section header (group divider)
        if section is not None:
            lines.append(f"{indent}// {'=' * 70}")
            lines.append(f"{indent}// {section}")
            lines.append(f"{indent}// {'=' * 70}")

        # Organism header
        lines.append(f"{indent}// {name} — {desc}")
        lines.append(f"{indent}// {band_disp} infection incidence band (x{band_fac})")

        # Mechanism lines
        for i, mech in enumerate(MECHANISMS):
            tier = tiers[i]
            base_rate = TIER_VALUES[tier]
            final_rate = base_rate * mult
            rate_str = fmt_rate(final_rate)
            config_key = f"bacteria_{key}_mechanism_{mech}_emergence_rate"
            line = f'{indent}map.insert("{config_key}".to_string(), {rate_str}); // tier {tier}'
            lines.append(line)

        lines.append("")  # blank line between organisms

    # Footer
    lines.append(f"// Generated 42 bacteria x 35 mechanisms = 1470 emergence rate values")
    return "\n".join(lines)


def main():
    config_path = os.path.join(os.path.dirname(__file__), "src", "config.rs")
    with open(config_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Find section boundaries
    start_marker = "        // --- Resistance Mechanisms Parameters ---"
    end_marker = "// Generated 42 bacteria"

    start_idx = content.find(start_marker)
    if start_idx == -1:
        print("ERROR: Could not find start marker")
        return

    # Find the end marker line and include the entire line
    end_idx = content.find(end_marker, start_idx)
    if end_idx == -1:
        print("ERROR: Could not find end marker")
        return
    # Extend to end of that line
    end_of_line = content.find("\n", end_idx)
    if end_of_line == -1:
        end_of_line = len(content)

    # Generate new section
    new_section = generate_section()

    # Replace
    new_content = content[:start_idx] + new_section + content[end_of_line:]

    with open(config_path, "w", encoding="utf-8") as f:
        f.write(new_content)

    # Count organisms and verify
    org_count = len(ORGS)
    print(f"SUCCESS: Replaced emergence rates section with {org_count} organisms x 35 mechanisms = {org_count * 35} values")
    print(f"File written: {config_path}")


if __name__ == "__main__":
    main()
