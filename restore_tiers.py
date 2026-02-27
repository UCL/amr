"""
Hybrid tier restoration script.

Strategy:
  1. K. pneumoniae, E. coli  — KEEP current (run 263034) rates unchanged
  2. All other 40 bacteria   — restore OLD tier assignments from .bak,
                                apply NEW incidence-based band multipliers
  3. Targeted bumps:
       S. aureus          +2 tiers  (was 0.2-1.1% IR at old rates; new band 6.7x lower)
       S. pneumoniae      +2 tiers  (was 0% IR at old rates; needs push above threshold)
       A. baumannii       +2 tiers  (was 0% IR at old rates; very low PDs)

  new_rate = TIER_VALUES[old_tier + bump] × new_band_multiplier
"""

import re
import math
import shutil

# ── TIER LADDER (half-log spacing) ──────────────────────────────────────
TIER_VALUES = [
    0.0,      # tier 0
    1e-09,    # tier 1
    1e-08,    # tier 2
    5e-08,    # tier 3
    1e-07,    # tier 4
    5e-07,    # tier 5
    1e-06,    # tier 6
    5e-06,    # tier 7
    1e-05,    # tier 8
    5e-05,    # tier 9
    1e-04,    # tier 10
    5e-04,    # tier 11
    1e-03,    # tier 12
]

# ── WHICH BACTERIA TO LEAVE ALONE ──────────────────────────────────────
KEEP_CURRENT = {'klebsiella_pneumoniae', 'escherichia_coli'}

# ── TIER BUMPS (applied BEFORE multiplication by new band) ─────────────
BUMPS = {
    'staphylococcus_aureus':   2,
    'streptococcus_pneumoniae': 2,
    'acinetobacter_baumannii':  2,
}

CONFIG  = 'src/config.rs'
BACKUP  = 'src/config.rs.bak'
SAFECOPY = 'src/config.rs.263034'   # safety copy of run-263034 config

# ── HELPERS ─────────────────────────────────────────────────────────────
def snap_to_tier(value):
    """Return (tier_number, snapped_tier_value) for the nearest tier."""
    if value <= 0:
        return 0, 0.0
    log_val = math.log10(value)
    best_i, best_dist = 1, float('inf')
    for i in range(1, len(TIER_VALUES)):
        d = abs(log_val - math.log10(TIER_VALUES[i]))
        if d < best_dist:
            best_dist = d
            best_i = i
    return best_i, TIER_VALUES[best_i]


def format_rate(value):
    """Format a rate in Rust-compatible scientific notation (e.g. 7.5e-05)."""
    if value == 0.0:
        return "0.0"
    exp = math.floor(math.log10(abs(value)))
    man = value / (10 ** exp)
    man = round(man, 1)
    if man >= 10.0:
        man /= 10.0
        exp += 1
    if man < 1.0:
        man *= 10.0
        exp -= 1
    if exp >= 0:
        return f"{man:.1f}e+{exp:02d}"
    else:
        return f"{man:.1f}e-{abs(exp):02d}"


# ── PARSE ───────────────────────────────────────────────────────────────
band_re = re.compile(r'//\s*Band\s+\d+\s*\(x([\d.]+)\)')
rate_re = re.compile(
    r'bacteria_(\w+?)_mechanism_(\w+)_emergence_rate.*?,\s*([\d.eE+-]+)\)'
)

def parse_bands_and_rates(text):
    """Return {bname: {band: float, rates: {mech: float}}}."""
    data = {}
    cur_band = None
    for line in text.split('\n'):
        bm = band_re.search(line)
        if bm:
            cur_band = float(bm.group(1))
        rm = rate_re.search(line)
        if rm:
            bname, mech, rate = rm.group(1), rm.group(2), float(rm.group(3))
            if bname not in data:
                data[bname] = {'band': cur_band, 'rates': {}}
            data[bname]['rates'][mech] = rate
    return data


# ── MAIN ────────────────────────────────────────────────────────────────
def main():
    # 0. Safety copy
    shutil.copy2(CONFIG, SAFECOPY)
    print(f"Safety copy saved → {SAFECOPY}")

    # 1. Read files
    with open(CONFIG,  'r') as f:
        config_lines = f.readlines()
    with open(BACKUP, 'r') as f:
        backup_text = f.read()
    config_text = ''.join(config_lines)

    old_data = parse_bands_and_rates(backup_text)
    new_data = parse_bands_and_rates(config_text)

    # 2. Compute replacements
    replacements = {}   # bname -> mech -> (new_rate, tier_num)
    summary = []
    warnings = []

    for bname in sorted(old_data):
        if bname in KEEP_CURRENT:
            cur_max = max(new_data[bname]['rates'].values()) if bname in new_data else 0
            summary.append({
                'name': bname, 'action': 'KEEP', 'ob': 0, 'nb': 0,
                'bump': 0, 'old_max': 0, 'new_max': cur_max, 'ratio': 0,
            })
            continue

        ob = old_data[bname]['band']
        nb = new_data[bname]['band'] if bname in new_data else None
        bump = BUMPS.get(bname, 0)

        if ob is None or nb is None:
            warnings.append(f"  SKIP {bname}: missing band (old={ob}, new={nb})")
            continue

        replacements[bname] = {}
        old_max = max(old_data[bname]['rates'].values())
        new_max = 0.0

        for mech, old_rate in old_data[bname]['rates'].items():
            if old_rate == 0.0:
                replacements[bname][mech] = (0.0, 0)
                continue

            raw_base = old_rate / ob
            tier_num, _ = snap_to_tier(raw_base)
            tier_num = min(tier_num + bump, 12)
            new_rate = TIER_VALUES[tier_num] * nb

            if new_rate > 0.05:
                warnings.append(
                    f"  WARN {bname}.{mech}: rate={new_rate:.2e} (capped at 5e-02)"
                )
                new_rate = 0.05

            replacements[bname][mech] = (new_rate, tier_num)
            new_max = max(new_max, new_rate)

        summary.append({
            'name': bname, 'action': 'RESTORE', 'ob': ob, 'nb': nb,
            'bump': bump, 'old_max': old_max, 'new_max': new_max,
            'ratio': new_max / old_max if old_max > 0 else 0,
        })

    # 3. Apply to config.rs
    out_lines = []
    modified = 0
    for line in config_lines:
        stripped = line.rstrip('\n')
        rm = rate_re.search(stripped)
        if rm:
            bname, mech = rm.group(1), rm.group(2)
            if bname in replacements and mech in replacements[bname]:
                new_rate, tier_num = replacements[bname][mech]
                rate_str = format_rate(new_rate)
                # Replace the numeric value between ", " and ");"
                new_stripped = re.sub(
                    r'(,\s*)([\d.eE+-]+)(\);)',
                    lambda m: m.group(1) + rate_str + m.group(3),
                    stripped,
                )
                # Update tier comment
                new_stripped = re.sub(
                    r'// tier \d+', f'// tier {tier_num}', new_stripped
                )
                out_lines.append(new_stripped + '\n')
                modified += 1
                continue
        out_lines.append(line)

    with open(CONFIG, 'w') as f:
        f.writelines(out_lines)

    # 4. Summary
    print(f"\nModified {modified} emergence-rate lines in {CONFIG}")
    if warnings:
        print("\nWarnings:")
        for w in warnings:
            print(w)

    print(f"\n{'Bacterium':<45} {'Action':<8} {'OldBand':>8} {'NewBand':>8} "
          f"{'Bump':>4} {'OldHallmark':>12} {'NewHallmark':>12} {'Ratio':>8}")
    print("─" * 120)
    for s in summary:
        if s['action'] == 'KEEP':
            print(f"{s['name']:<45} {'KEEP':<8} {'—':>8} {'—':>8} "
                  f"{'—':>4} {'—':>12} {s['new_max']:>12.1e} {'—':>8}")
        else:
            print(f"{s['name']:<45} {'RESTORE':<8} {s['ob']:>8.1f} {s['nb']:>8.1f} "
                  f"{s['bump']:>+4d} {s['old_max']:>12.1e} {s['new_max']:>12.1e} "
                  f"{s['ratio']:>8.2f}")


if __name__ == '__main__':
    main()
