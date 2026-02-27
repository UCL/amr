#!/usr/bin/env python3
"""
adjust_mechanism_700021.py — Ceiling-test calibration for run 700021.

Insight from run 700021: emergence is drug-exposure-gated. A mechanism can only
emerge when the patient is treated with a drug the mechanism applies to. Beta-
lactamase mechanisms cover ~12-20 drugs and emerge easily; non-beta-lactam
mechanisms (gyrA → cipro/oflox only, TetM → 3 drugs, 16S → 3 drugs, etc.)
rarely get a chance to roll.

Strategy:
  NON-BETA-LACTAM mechanisms → set to 5e-3 (ceiling test) for most drug-
  specific resistance mechanisms, and 1e-3–3e-3 for support mechanisms.
  BETA-LACTAM mechanisms → per-bacterium adjustments (revert E. coli,
  reduce K. pneumoniae, etc.)

Applies to all 22 bacteria with mean |Δ| ≥ 20 from run 700021.
MDR TB excluded (rifampicin resistance is pre-assumed).

Safety:
  - NEVER enables tier-0 (0.0) mechanisms.
  - Caps all values at 1e-2.
  - Preserves as_yet_unknown_1/2/3 as 0.0 where they are 0.0.
"""

import re
import shutil
import sys
from pathlib import Path

CONFIG = Path("src/config.rs")
BACKUP_112750 = Path("src/config.rs.112750")
SAFETY_BACKUP = Path("src/config.rs.700021")

# ──────────────────────────────────────────────────────────────────────
# Mechanism categories
# ──────────────────────────────────────────────────────────────────────

# Non-beta-lactam mechanisms: set to ceiling test values
# Key: mechanism name → ceiling rate
CEILING_RATES = {
    # FQ mechanisms (cipro/oflox only for gyrA_primary, 4 FQs for parc)
    "mutation_gyra_primary":            5e-3,
    "mutation_gyra_parc_secondary":     5e-3,
    "protection_qnr":                   5e-3,
    # Tetracycline (3 drugs)
    "protection_tet_m":                 5e-3,
    # Aminoglycoside (3 drugs)
    "enzyme_16s_rrmt":                  5e-3,
    # Macrolide/clinda/Q-D
    "target_site_erm_b":                5e-3,
    # Folate (sulfanilamide, trim_sulf)
    "mutation_folate_pathway":          5e-3,
    # Glycopeptides
    "target_site_van_a":                5e-3,
    "target_site_van_b":                5e-3,
    # Efflux: AcrAB (tet/chlor/cipro), MexXY (AG/tet/cipro), Global efflux
    "efflux_acrab_tolc":                5e-3,
    "efflux_mexxy_oprm":                3e-3,
    "global_efflux_pump":               3e-3,
    # Chloramphenicol
    "enzyme_cat":                       3e-3,
    # Oxazolidinones/phenicols/lincos (Cfr PhLOPSA)
    "target_site_cfr":                  5e-3,
    # Colistin
    "modification_mcr_1":               3e-3,
    # Daptomycin (MprF) / Fusidic acid (FusB) — Gram-positive-only in practice
    "mutation_mpr_f":                   3e-3,
    "protection_fus_b":                 3e-3,
    # Porin loss — broad
    "porin_loss_ompk35_36":             3e-3,
    "porin_loss_oprd":                  3e-3,
    "global_porin_loss":                3e-3,
    # Metronidazole/nitrofurantoin
    "mutation_nitroreductase":          3e-3,
    # Fosfomycin
    "enzyme_fos_a":                     1e-3,
    # Fidaxomicin (NOT rifampicin)
    "mutation_rpo_b":                   1e-3,
    # PBP2a — applies to beta-lactams but only for Gram-positives
    # Handled specially per-bacterium below
}

# Beta-lactam mechanisms — these are NOT set to ceiling; handled per-bacterium
BETA_LACTAM_MECHANISMS = {
    "enzyme_esbl_ctx_m",
    "enzyme_esbl_tem",
    "enzyme_esbl_shv",
    "enzyme_ampc_cmy",
    "enzyme_ampc_dha",
    "enzyme_kpc",
    "enzyme_ndm_vim",
    "enzyme_oxa_48",
    "target_site_pbp2a_meca",  # BL mechanism for Gram-positives
}

# Mechanisms that must ALWAYS stay 0.0 if currently 0.0
MUST_STAY_ZERO = {
    "as_yet_unknown_1",
    "as_yet_unknown_2",
    "as_yet_unknown_3",
}

# ──────────────────────────────────────────────────────────────────────
# Per-bacterium config: slug, and beta-lactam strategy
# ──────────────────────────────────────────────────────────────────────
# Beta-lactam strategies:
#   "ceiling"     — set all non-zero BL mechanisms to 5e-3 (for all-0% bacteria)
#   "keep"        — don't change BL mechanisms
#   "revert_x1.5" — revert to backup values × 1.5 (E. coli)
#   "divide_N"    — divide all BL mechanisms by N
#   "custom"      — per-mechanism dict of multipliers or absolute values
#
# Special non-BL overrides:
#   "non_bl_overrides" — dict of mechanism → rate (overrides ceiling)

BACTERIA_CONFIG = {
    # ── E. coli: BL overshoot (0% → 100% at ×3). Revert to backup ×1.5 ──
    "escherichia_coli": {
        "bl_strategy": "revert_x1.5",
        # non-BL: ceiling test (default)
    },

    # ── K. pneumoniae: BL+carbapenems overshooting (42-44%). Halve. ──
    "klebsiella_pneumoniae": {
        "bl_strategy": "divide_2",
        # FQ/AG/Tet at 0% → ceiling test
    },

    # ── N. gonorrhoeae: BL at ~93% (target 3-30%). Reduce ESBL_TEM ÷5.
    #    TetM at 93% (target 25-35%) → reduce to 5e-5 instead of ceiling.
    #    FQ at 27% (target 80%) → ceiling. Macrolide at 0% → ceiling. ──
    "neisseria_gonorrhoeae": {
        "bl_strategy": "custom",
        "bl_custom": {
            "enzyme_esbl_tem": ("divide", 5),      # 2.3e-05 → 4.6e-06
            "enzyme_esbl_ctx_m": ("keep", None),
            "enzyme_esbl_shv": ("keep", None),
            "enzyme_ampc_cmy": ("keep", None),
            "enzyme_ampc_dha": ("keep", None),
            "enzyme_ndm_vim": ("keep", None),
        },
        "non_bl_overrides": {
            "protection_tet_m": 5e-5,  # Currently 2.2e-04 → 93% tet, target 25-35%
        },
    },

    # ── A. baumannii: BL at 75% (close for some drugs). Reduce OXA/NDM.
    #    FQ/AG/Tet collapsed → ceiling test. ──
    "acinetobacter_baumannii": {
        "bl_strategy": "custom",
        "bl_custom": {
            "enzyme_oxa_48": ("divide", 3),     # 3.0e-04 → 1.0e-04
            "enzyme_ndm_vim": ("divide", 3),    # 5.0e-05 → 1.7e-05
            "enzyme_kpc": ("divide", 2),        # 5.0e-06 → 2.5e-06
            "enzyme_esbl_ctx_m": ("keep", None),
            "enzyme_esbl_tem": ("keep", None),
            "enzyme_esbl_shv": ("keep", None),
            "enzyme_ampc_cmy": ("keep", None),
            "enzyme_ampc_dha": ("keep", None),
        },
    },

    # ── S. aureus: BL at 3.67% (target 35-60%). Boost PBP2a ×10.
    #    FQ/mac/AG/tet all 0% → ceiling test. ──
    "staphylococcus_aureus": {
        "bl_strategy": "custom",
        "bl_custom": {
            "target_site_pbp2a_meca": ("multiply", 10),  # 7.5e-05 → 7.5e-04
        },
    },

    # ── B. fragilis: carbapenems +60pp overshoot. Reduce carbapenemases ÷5.
    #    FQ/AG at 0% → ceiling. ──
    "bacteroides_fragilis": {
        "bl_strategy": "custom",
        "bl_custom": {
            "enzyme_kpc": ("divide", 10),       # 1.1e-07 → 1.1e-08
            "enzyme_ndm_vim": ("divide", 5),    # 5.7e-06 → 1.1e-06
            "enzyme_oxa_48": ("divide", 10),    # 1.1e-07 → 1.1e-08
            "enzyme_esbl_ctx_m": ("keep", None),
            "enzyme_esbl_tem": ("keep", None),
            "enzyme_esbl_shv": ("keep", None),
            "enzyme_ampc_cmy": ("keep", None),
            "enzyme_ampc_dha": ("keep", None),
        },
    },

    # ── E. cloacae: BL at 7.46% (improving). Keep BL. ──
    "enterobacter_cloacae": {
        "bl_strategy": "keep",
    },

    # ── E. faecium: all 0%, 56 PDs. Max everything. ──
    "enterococcus_faecium": {
        "bl_strategy": "ceiling",  # PBP2a (Gram-pos BL mechanism) → ceiling
    },

    # ── Shigella spp.: all 0%, 3701 PDs. Ceiling everything. ──
    "shigella_spp.": {
        "bl_strategy": "ceiling",
    },

    # ── C. jejuni: all 0%, 3312 PDs. No BL mechanisms (all zero). ──
    "campylobacter_jejuni": {
        "bl_strategy": "keep",  # All BL already 0.0
    },

    # ── iNTS: all 0%, 445 PDs. Ceiling everything. ──
    "invasive_non-typhoidal_salmonella_spp.": {
        "bl_strategy": "ceiling",
    },

    # ── S. Paratyphi A: all 0%, 171 PDs. Ceiling everything. ──
    "salmonella_enterica_serovar_paratyphi_a": {
        "bl_strategy": "ceiling",
    },

    # ── S. Typhi: all 0%, 673 PDs. Ceiling everything. ──
    "salmonella_enterica_serovar_typhi": {
        "bl_strategy": "ceiling",
    },

    # ── Citrobacter spp.: all 0%, 150 PDs. Ceiling everything. ──
    "citrobacter_spp.": {
        "bl_strategy": "ceiling",
    },

    # ── Serratia spp.: all 0%, 57 PDs. Ceiling everything. ──
    "serratia_spp.": {
        "bl_strategy": "ceiling",
    },

    # ── Morganella spp.: all 0%, 77 PDs. Ceiling everything. ──
    "morganella_spp.": {
        "bl_strategy": "ceiling",
    },

    # ── Proteus spp.: all 0%, 400 PDs. Ceiling everything. ──
    "proteus_spp.": {
        "bl_strategy": "ceiling",
    },

    # ── E. faecalis: all 0%, 214 PDs. Ceiling everything. ──
    "enterococcus_faecalis": {
        "bl_strategy": "ceiling",
    },

    # ── M. genitalium: all 0%, ~200 PDs. No BL (all zero). ──
    "mycoplasma_genitalium": {
        "bl_strategy": "keep",  # All BL already 0.0
    },

    # ── H. pylori: all 0%, ~1000+ PDs. No BL (all zero). ──
    "helicobacter_pylori": {
        "bl_strategy": "keep",  # All BL already 0.0
    },
}

# Maximum allowed emergence rate (safety cap)
MAX_RATE = 1e-2

# ──────────────────────────────────────────────────────────────────────
# Regex pattern for emergence rate lines
# ──────────────────────────────────────────────────────────────────────
RATE_PATTERN = re.compile(
    r'(map\.insert\("bacteria_)([^"]+?)(_mechanism_)([^"]+?)(_emergence_rate"\.to_string\(\),\s*)'
    r'([0-9eE.+-]+)'
    r'(\);.*)'
)


def read_backup_rates(backup_path: Path) -> dict:
    """Read all emergence rates from backup file into a dict keyed by (bacteria, mechanism)."""
    rates = {}
    with open(backup_path, "r") as f:
        for line in f:
            m = RATE_PATTERN.search(line)
            if m:
                bacteria = m.group(2)
                mechanism = m.group(4)
                rate = float(m.group(6))
                rates[(bacteria, mechanism)] = rate
    return rates


def get_tier_comment(rate: float) -> str:
    """Return a tier comment string for the given rate."""
    if rate == 0.0:
        return "tier 0"
    tiers = [
        (1e-9, 1), (1e-8, 2), (5e-8, 3), (1e-7, 4), (5e-7, 5),
        (1e-6, 6), (5e-6, 7), (1e-5, 8), (5e-5, 9), (1e-4, 10),
        (5e-4, 11), (1e-3, 12),
    ]
    # Find nearest tier
    best_tier = 12
    best_dist = abs(rate - 1e-3)
    for threshold, tier in tiers:
        dist = abs(rate - threshold)
        if dist < best_dist:
            best_dist = dist
            best_tier = tier
    if rate > 1e-3:
        return f">T12 ({rate:.1e})"
    return f"tier {best_tier}"


def format_rate(rate: float) -> str:
    """Format a rate value in scientific notation matching config.rs style."""
    if rate == 0.0:
        return "0.0"
    # Use engineering-style: N.Ne-M
    exp = 0
    val = rate
    if val >= 1.0:
        return f"{val:.1f}"
    while val < 1.0:
        val *= 10
        exp += 1
    # val is now between 1.0 and 10.0
    return f"{val:.1f}e-{exp:02d}"


def compute_new_rate(bacteria: str, mechanism: str, current_rate: float,
                     config: dict, backup_rates: dict) -> float:
    """Compute the new rate for a given bacteria/mechanism combo."""
    # Never change zero rates
    if current_rate == 0.0:
        return 0.0

    # Never change as_yet_unknown
    if mechanism in MUST_STAY_ZERO:
        return current_rate  # Should already be 0.0, but safety

    # Check if this is a BL mechanism
    is_bl = mechanism in BETA_LACTAM_MECHANISMS

    if is_bl:
        strategy = config.get("bl_strategy", "keep")

        if strategy == "keep":
            return current_rate

        elif strategy == "ceiling":
            new_rate = 5e-3
            return min(new_rate, MAX_RATE)

        elif strategy == "revert_x1.5":
            backup_rate = backup_rates.get((bacteria, mechanism))
            if backup_rate is not None and backup_rate > 0.0:
                new_rate = backup_rate * 1.5
                return min(new_rate, MAX_RATE)
            return current_rate

        elif strategy.startswith("divide_"):
            divisor = float(strategy.split("_")[1])
            new_rate = current_rate / divisor
            return max(new_rate, 1e-12)  # Don't go to zero

        elif strategy == "custom":
            bl_custom = config.get("bl_custom", {})
            if mechanism in bl_custom:
                action, param = bl_custom[mechanism]
                if action == "keep":
                    return current_rate
                elif action == "divide":
                    return max(current_rate / param, 1e-12)
                elif action == "multiply":
                    return min(current_rate * param, MAX_RATE)
                elif action == "set":
                    return min(param, MAX_RATE)
            # BL mechanism not in custom dict → keep
            return current_rate

        return current_rate

    else:
        # Non-BL mechanism: check for per-bacterium overrides first
        overrides = config.get("non_bl_overrides", {})
        if mechanism in overrides:
            new_rate = overrides[mechanism]
            return min(new_rate, MAX_RATE) if new_rate > 0.0 else 0.0

        # Otherwise, apply ceiling rate from the category table
        if mechanism in CEILING_RATES:
            new_rate = CEILING_RATES[mechanism]
            return min(new_rate, MAX_RATE)

        # Unknown mechanism — keep current
        return current_rate


def main():
    if not CONFIG.exists():
        print(f"ERROR: {CONFIG} not found.")
        sys.exit(1)
    if not BACKUP_112750.exists():
        print(f"WARNING: Backup {BACKUP_112750} not found. Revert strategies will use current values.")
        backup_rates = {}
    else:
        backup_rates = read_backup_rates(BACKUP_112750)
        print(f"Read {len(backup_rates)} backup rates from {BACKUP_112750}")

    # Safety backup
    if not SAFETY_BACKUP.exists():
        shutil.copy2(CONFIG, SAFETY_BACKUP)
        print(f"Created safety backup: {SAFETY_BACKUP}")
    else:
        print(f"Safety backup already exists: {SAFETY_BACKUP}")

    # Read current config
    with open(CONFIG, "r") as f:
        lines = f.readlines()

    modified_count = 0
    unchanged_count = 0
    skipped_zero = 0
    bacteria_modified = set()

    new_lines = []
    for line in lines:
        m = RATE_PATTERN.search(line)
        if not m:
            new_lines.append(line)
            continue

        bacteria = m.group(2)
        mechanism = m.group(4)
        current_rate = float(m.group(6))

        # Only process bacteria in our config
        if bacteria not in BACTERIA_CONFIG:
            new_lines.append(line)
            continue

        config = BACTERIA_CONFIG[bacteria]

        # Skip zero-rate mechanisms
        if current_rate == 0.0:
            new_lines.append(line)
            skipped_zero += 1
            continue

        new_rate = compute_new_rate(bacteria, mechanism, current_rate,
                                    config, backup_rates)

        if abs(new_rate - current_rate) / max(current_rate, 1e-15) < 0.001:
            # No meaningful change
            new_lines.append(line)
            unchanged_count += 1
            continue

        # Reconstruct line with new rate and updated tier comment
        prefix = m.group(1) + m.group(2) + m.group(3) + m.group(4) + m.group(5)
        suffix_raw = m.group(7)

        # Update tier comment in suffix
        tier_str = get_tier_comment(new_rate)
        # Replace existing tier comment
        suffix_updated = re.sub(
            r'//\s*tier\s+\d+|//\s*>T\d+\s*\([^)]+\)',
            f'// {tier_str}',
            suffix_raw
        )

        new_rate_str = format_rate(new_rate)
        new_line_content = f"{prefix}{new_rate_str}{suffix_updated}"

        # Preserve original indentation
        leading_ws = line[:len(line) - len(line.lstrip())]
        new_line = leading_ws + new_line_content.lstrip() + "\n"

        new_lines.append(new_line)
        modified_count += 1
        bacteria_modified.add(bacteria)

    # Write out
    with open(CONFIG, "w") as f:
        f.writelines(new_lines)

    print(f"\n{'='*60}")
    print(f"RESULTS")
    print(f"{'='*60}")
    print(f"Lines modified:   {modified_count}")
    print(f"Lines unchanged:  {unchanged_count}")
    print(f"Zero-rate skipped: {skipped_zero}")
    print(f"Bacteria touched: {len(bacteria_modified)}")
    print(f"\nBacteria modified:")
    for b in sorted(bacteria_modified):
        print(f"  - {b}")

    # Spot-check critical values
    print(f"\n{'='*60}")
    print(f"SPOT-CHECK VERIFICATION")
    print(f"{'='*60}")
    spot_checks = [
        ("escherichia_coli", "enzyme_esbl_ctx_m", "~4.5e-08 (backup×1.5)"),
        ("escherichia_coli", "mutation_gyra_primary", "5.0e-03 (ceiling)"),
        ("escherichia_coli", "as_yet_unknown_1", "0.0 (must stay zero)"),
        ("neisseria_gonorrhoeae", "protection_tet_m", "5.0e-05 (override)"),
        ("neisseria_gonorrhoeae", "enzyme_esbl_tem", "~4.6e-06 (÷5)"),
        ("neisseria_gonorrhoeae", "mutation_gyra_primary", "5.0e-03 (ceiling)"),
        ("staphylococcus_aureus", "target_site_pbp2a_meca", "7.5e-04 (×10)"),
        ("staphylococcus_aureus", "enzyme_esbl_ctx_m", "0.0 (stays zero)"),
        ("klebsiella_pneumoniae", "enzyme_esbl_ctx_m", "~2.5e-05 (÷2)"),
        ("klebsiella_pneumoniae", "mutation_gyra_primary", "5.0e-03 (ceiling)"),
        ("shigella_spp.", "mutation_gyra_primary", "5.0e-03 (ceiling)"),
        ("shigella_spp.", "enzyme_esbl_ctx_m", "5.0e-03 (ceiling BL)"),
        ("campylobacter_jejuni", "mutation_gyra_primary", "5.0e-03 (ceiling)"),
        ("campylobacter_jejuni", "enzyme_esbl_ctx_m", "0.0 (stays zero)"),
        ("acinetobacter_baumannii", "enzyme_oxa_48", "~1.0e-04 (÷3)"),
        ("bacteroides_fragilis", "enzyme_kpc", "~1.1e-08 (÷10)"),
        ("enterococcus_faecium", "target_site_van_a", "5.0e-03 (ceiling)"),
        ("enterococcus_faecalis", "target_site_erm_b", "5.0e-03 (ceiling)"),
    ]

    with open(CONFIG, "r") as f:
        final_content = f.read()

    for bact, mech, expected in spot_checks:
        pattern = re.compile(
            rf'bacteria_{re.escape(bact)}_mechanism_{re.escape(mech)}_emergence_rate.*?,\s*([0-9eE.+-]+)\)'
        )
        match = pattern.search(final_content)
        if match:
            actual = match.group(1)
            print(f"  {bact:45s} {mech:35s} → {actual:12s}  (expected: {expected})")
        else:
            print(f"  {bact:45s} {mech:35s} → NOT FOUND")


if __name__ == "__main__":
    main()
