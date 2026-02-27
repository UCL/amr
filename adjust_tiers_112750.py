"""
Tier adjustments based on run 112750 calibration results.

Changes:
  S. pneumoniae  -1 tier  (58% flat everywhere, targets ~20%)
  N. gonorrhoeae +2 tiers (0% everywhere, selective boost failed)
  E. faecium     +2 tiers (0% everywhere, +1 wasn't enough)
  P. aeruginosa  +1 tier  (lost differentiation from 180383)
  S. aureus      +1 tier  (great shape 0.3-25%, but targets 35-45%)
"""

import re

TIER_VALUES = [
    0.0, 1e-9, 1e-8, 5e-8, 1e-7, 5e-7,
    1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3,
]
MAX_TIER = len(TIER_VALUES) - 1  # 12

# bacterium_key -> tier delta
ADJUSTMENTS = {
    "streptococcus_pneumoniae": -1,
    "neisseria_gonorrhoeae":    +2,
    "enterococcus_faecium":     +2,
    "pseudomonas_aeruginosa":   +1,
    "staphylococcus_aureus":    +1,
}

CONFIG = r"src\config.rs"

with open(CONFIG, "r") as f:
    lines = f.readlines()

# Pattern matches emergence rate lines, using [\w.] for bacteria names with dots
pattern = re.compile(
    r'^(\s*map\.insert\("bacteria_)([\w.]+)(_mechanism_[\w.]+_emergence_rate"'
    r'\.to_string\(\),\s*)'
    r'([\d.eE\-+]+)'
    r'(\);?\s*//\s*tier\s+)(\d+)(.*\n)$'
)

modified = 0
skipped_capped = 0

for i, line in enumerate(lines):
    m = pattern.match(line)
    if not m:
        continue

    bact = m.group(2)
    if bact not in ADJUSTMENTS:
        continue

    delta = ADJUSTMENTS[bact]
    old_tier = int(m.group(6))
    new_tier = old_tier + delta

    # Clamp to valid range
    if new_tier < 0:
        new_tier = 0
        skipped_capped += 1
    elif new_tier > MAX_TIER:
        new_tier = MAX_TIER
        skipped_capped += 1

    if new_tier == old_tier:
        continue

    new_val = TIER_VALUES[new_tier]
    # Format value consistently
    if new_val == 0.0:
        val_str = "0.0"
    else:
        val_str = f"{new_val:.0e}"

    lines[i] = f"{m.group(1)}{bact}{m.group(3)}{val_str}{m.group(5)}{new_tier}{m.group(7)}"
    modified += 1

with open(CONFIG, "w") as f:
    f.writelines(lines)

print(f"Modified {modified} emergence rate lines.")
if skipped_capped:
    print(f"  ({skipped_capped} lines hit tier floor/ceiling)")
print("Done.")
