import re

with open(r'c:\Users\w3sth\rust_amr_project\output_graphs\calibration_summary_411919.txt', 'r') as f:
    lines = f.readlines()

# Find the Resistance Benchmarks section (not Microbiome)
start = None
for i, line in enumerate(lines):
    if 'Resistance Benchmarks (percent resistant)' in line:
        start = i + 2  # skip header line
        break

bacteria_data = {}

for i in range(start, len(lines)):
    line = lines[i].rstrip()
    if not line or line.startswith('Note:') or line.startswith('---') or line.startswith('='):
        continue
    parts = re.split(r'\s{2,}', line.strip())
    if len(parts) < 5:
        continue
    bact = parts[0].strip()
    drug = parts[1].strip()
    drug_class = parts[2].strip()
    sim_str = parts[3].strip()
    target_str = parts[4].strip()

    infected_pd = '---'
    if len(parts) > 8:
        infected_pd = parts[8].strip()

    bact_lower = bact.lower()
    if bact_lower not in bacteria_data:
        bacteria_data[bact_lower] = {'name': bact, 'drugs': [], 'infected_pd': infected_pd}

    if infected_pd != '---':
        bacteria_data[bact_lower]['infected_pd'] = infected_pd

    bacteria_data[bact_lower]['drugs'].append({
        'drug': drug, 'drug_class': drug_class,
        'sim': sim_str, 'target': target_str
    })

print(f'Total bacteria found: {len(bacteria_data)}')
print()

# Summary table
print("SUMMARY TABLE")
print("=" * 120)
fmt = "{:<45} {:>12} {:>15} {:>20} {:>10}"
print(fmt.format("Bacteria", "Infected PDs", "Pattern", "Sim Range", "vs Target"))
print("-" * 120)

key_bacteria = {
    'escherichia coli', 'streptococcus pneumoniae', 'neisseria gonorrhoeae',
    'staphylococcus aureus', 'klebsiella pneumoniae', 'acinetobacter baumannii',
    'bacteroides fragilis', 'pseudomonas aeruginosa', 'enterococcus faecium',
    'haemophilus influenzae', 'shigella spp.', 'campylobacter jejuni'
}

for bact_key in sorted(bacteria_data.keys()):
    info = bacteria_data[bact_key]
    sims = []
    targets = []
    for d in info['drugs']:
        if d['sim'] != '---':
            try:
                sims.append(float(d['sim']))
            except:
                pass
        if d['target'] != '---':
            try:
                targets.append(float(d['target']))
            except:
                pass

    all_zero = all(s == 0.0 for s in sims) if sims else True
    if not sims:
        pattern = 'ALL ---'
        sim_range = 'N/A'
    elif all_zero:
        pattern = 'All 0%'
        sim_range = '0%'
    elif len(set(sims)) == 1:
        pattern = 'Uniform ' + f'{sims[0]:.2f}%'
        sim_range = f'{sims[0]:.2f}%'
    else:
        non_zero = [s for s in sims if s > 0]
        if non_zero:
            pattern = 'Varied'
            sim_range = f'{min(sims):.1f}-{max(sims):.1f}%'
        else:
            pattern = 'All 0%'
            sim_range = '0%'

    if sims and targets:
        deltas = [s - t for s, t in zip(sims, targets)]
        avg_delta = sum(deltas) / len(deltas)
        if avg_delta < -5:
            direction = 'UNDER'
        elif avg_delta > 5:
            direction = 'OVER'
        else:
            direction = 'CLOSE'
    else:
        direction = 'N/A'

    pd_val = info['infected_pd']
    print(fmt.format(info['name'], pd_val, pattern, sim_range, direction))

# Detailed breakdown for key bacteria
print("\n\n")
print("DETAILED DRUG-BY-DRUG BREAKDOWN FOR KEY BACTERIA")
print("=" * 120)

for bact_key in sorted(bacteria_data.keys()):
    if bact_key not in key_bacteria:
        continue
    info = bacteria_data[bact_key]
    print(f"\n{'='*80}")
    print(f"  {info['name']}  (Infected person-days: {info['infected_pd']})")
    print(f"{'='*80}")
    dfmt = "  {:<28} {:<35} {:>8} {:>8} {:>8}"
    print(dfmt.format("Drug", "Drug Class", "Sim%", "Target%", "Delta"))
    print("  " + "-" * 95)
    for d in info['drugs']:
        sim_s = d['sim']
        tgt_s = d['target']
        if sim_s == '---' and tgt_s == '---':
            delta_s = '---'
        else:
            try:
                delta_s = f'{float(sim_s) - float(tgt_s):+.1f}'
            except:
                delta_s = '---'
        print(dfmt.format(d['drug'], d['drug_class'], sim_s, tgt_s, delta_s))
