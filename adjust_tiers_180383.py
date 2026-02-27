"""
Targeted tier adjustments based on run 180383 calibration results.

For each rate R at tier T with band multiplier B:
  R = TIER_VALUES[T] × B
  New R = TIER_VALUES[T + delta] × B = R × TIER_VALUES[T+delta] / TIER_VALUES[T]

Tier 0 (rate=0.0) rates are never modified (biologically excluded).
"""
import re
import math
import shutil

CONFIG  = 'src/config.rs'
SAFECOPY = 'src/config.rs.180383'   # safety copy of current config

TIER_VALUES = [
    0.0,      # 0
    1e-09,    # 1
    1e-08,    # 2
    5e-08,    # 3
    1e-07,    # 4
    5e-07,    # 5
    1e-06,    # 6
    5e-06,    # 7
    1e-05,    # 8
    5e-05,    # 9
    1e-04,    # 10
    5e-04,    # 11
    1e-03,    # 12
]

# ── GLOBAL TIER ADJUSTMENTS ────────────────────────────────────────────
# Applied to ALL mechanisms for these bacteria
GLOBAL_ADJUSTMENTS = {
    'streptococcus_pneumoniae': -2,   # 98.6% → undo the +2 bump
    'shigella_spp.':            -2,   # 99%   → same saturation problem
    'escherichia_coli':         -2,   # 100%  → still saturated, cut harder
    'staphylococcus_aureus':    -1,   # 82% flat → try to find middle
    'acinetobacter_baumannii':  -1,   # 97%/46% → reduce overshoot gently
    'klebsiella_pneumoniae':    +1,   # 0%   → boost above threshold
    'enterococcus_faecium':     +1,   # 0%   → cautious boost
}

# ── N. GONORRHOEAE SELECTIVE BOOST ──────────────────────────────────────
# Beta-lactam resistance is 77% (working) — leave enzyme/porin/efflux alone.
# FQ/macrolide/tet/aminoglycoside are at 0% vs targets 35–80% — boost those.
NGON_SELECTIVE = {
    'mutation_gyra_primary':            +2,   # FQ
    'mutation_gyra_parc_secondary':     +2,   # FQ
    'protection_qnr':                   +2,   # FQ
    'target_site_erm_b':                +2,   # macrolide
    'target_site_cfr':                  +2,   # oxazolidinone
    'protection_tet_m':                 +2,   # tetracycline
    'enzyme_16s_rrmt':                  +2,   # aminoglycoside
    'enzyme_cat':                       +2,   # chloramphenicol
    'mutation_folate_pathway':          +2,   # trimethoprim/sulfonamide
    'modification_mcr_1':               +2,   # colistin
    'mutation_nitroreductase':          +2,   # nitrofurantoin
    'mutation_rpo_b':                   +2,   # rifampicin
}

# ── HELPERS ─────────────────────────────────────────────────────────────
rate_re = re.compile(
    r'(bacteria_(\w+?)_mechanism_(\w+)_emergence_rate.*?,\s*)([\d.eE+-]+)(\);.*?//\s*tier\s+)(\d+)'
)

def format_rate(value):
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
    return f"{man:.1f}e-{abs(exp):02d}" if exp < 0 else f"{man:.1f}e+{exp:02d}"


def get_delta(bname, mech):
    """Return tier delta for a given bacterium-mechanism pair."""
    if bname in GLOBAL_ADJUSTMENTS:
        return GLOBAL_ADJUSTMENTS[bname]
    if bname == 'neisseria_gonorrhoeae' and mech in NGON_SELECTIVE:
        return NGON_SELECTIVE[mech]
    return 0  # no change


# ── MAIN ────────────────────────────────────────────────────────────────
def main():
    shutil.copy2(CONFIG, SAFECOPY)
    print(f"Safety copy saved → {SAFECOPY}")

    with open(CONFIG, 'r') as f:
        lines = f.readlines()

    modified = 0
    summary = {}  # bname -> list of (mech, old_tier, new_tier, old_rate, new_rate)

    out_lines = []
    for line in lines:
        m = rate_re.search(line)
        if m:
            prefix    = m.group(1)
            bname     = m.group(2)
            mech      = m.group(3)
            old_rate  = float(m.group(4))
            mid       = m.group(5)
            old_tier  = int(m.group(6))

            delta = get_delta(bname, mech)

            if delta != 0 and old_tier > 0 and old_rate > 0:
                new_tier = max(1, min(12, old_tier + delta))
                ratio = TIER_VALUES[new_tier] / TIER_VALUES[old_tier]
                new_rate = old_rate * ratio

                new_rate_str = format_rate(new_rate)

                # Rebuild the line
                new_line = re.sub(
                    r'(,\s*)([\d.eE+-]+)(\);)',
                    lambda x: x.group(1) + new_rate_str + x.group(3),
                    line.rstrip('\n')
                )
                new_line = re.sub(
                    r'// tier \d+',
                    f'// tier {new_tier}',
                    new_line
                )
                out_lines.append(new_line + '\n')
                modified += 1

                if bname not in summary:
                    summary[bname] = []
                summary[bname].append((mech, old_tier, new_tier, old_rate, new_rate))
            else:
                out_lines.append(line)
        else:
            out_lines.append(line)

    with open(CONFIG, 'w') as f:
        f.writelines(out_lines)

    # Report
    print(f"\nModified {modified} emergence-rate lines\n")
    print(f"{'Bacterium':<45} {'Delta':>5} {'Mechs':>5} "
          f"{'Old Hallmark':>13} {'New Hallmark':>13} {'Ratio':>8}")
    print("─" * 100)
    for bname in sorted(summary):
        entries = summary[bname]
        old_max = max(e[3] for e in entries)
        new_max = max(e[4] for e in entries)
        delta = entries[0][2] - entries[0][1]  # same delta for all
        ratio = new_max / old_max if old_max > 0 else 0
        print(f"{bname:<45} {delta:>+5d} {len(entries):>5} "
              f"{old_max:>13.2e} {new_max:>13.2e} {ratio:>8.2f}")

    # Scope verification
    print("\n── Scope Verification ──")
    with open(SAFECOPY, 'r') as f:
        orig = f.readlines()
    sec_start = sec_end = None
    for i, ln in enumerate(orig):
        if 'BACTERIA-MECHANISM-SPECIFIC EMERGENCE RATES' in ln:
            sec_start = i
            break
    for i in range(len(orig)-1, -1, -1):
        if 'as_yet_unknown_3_emergence_rate' in orig[i]:
            sec_end = i
            break
    diffs_outside = 0
    for i in range(sec_start):
        if orig[i] != out_lines[i]:
            diffs_outside += 1
    offset = len(out_lines) - len(orig)
    for i in range(sec_end+1, len(orig)):
        j = i + offset
        if j < len(out_lines) and orig[i] != out_lines[j]:
            diffs_outside += 1
    print(f"Lines: {len(orig)} → {len(out_lines)}, delta={offset}")
    print(f"Diffs outside emergence section: {diffs_outside}")
    if diffs_outside == 0:
        print("CONFIRMED: Changes only in emergence rates section.")


if __name__ == '__main__':
    main()
