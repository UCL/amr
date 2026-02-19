"""
Analyze individual drug usage within Polymyxins & Others category
Using streaming to handle large CSV files efficiently
"""
import csv

# Read the CSV file
csv_path = 'amr_simulation_output_analysis_outputs/simulation_summary_051228.csv'
print(f"Loading {csv_path}...")

# First, read just the header
with open(csv_path, 'r') as f:
    reader = csv.reader(f)
    headers = next(reader)
    print(f"Found {len(headers)} columns")

# Define the target drugs
target_drugs = [
    'colistin_currently_on_drug',
    'rifampicin_currently_on_drug',
    'chloramphenicol_currently_on_drug',
    'fusidic_a_currently_on_drug',
    'retapamulin_currently_on_drug',
    'quinu_dalfo_currently_on_drug'
]

# Map drug names to column indices
drug_indices = {}
year_idx = None
total_drug_idx = None

for i, header in enumerate(headers):
    if header in target_drugs:
        drug_indices[header] = i
    if header == 'year':
        year_idx = i
    if header == 'currently_taking_drug_count':
        total_drug_idx = i

print(f"\nFound {len(drug_indices)} drug columns:")
for drug in drug_indices:
    print(f"  - {drug}")

missing_drugs = [d for d in target_drugs if d not in drug_indices]
if missing_drugs:
    print(f"\nMissing {len(missing_drugs)} drugs:")
    for drug in missing_drugs:
        print(f"  - {drug}")

# Stream through the file and accumulate statistics
print("\nStreaming data...")
drug_sums = {drug: 0.0 for drug in drug_indices}
drug_counts = {drug: 0 for drug in drug_indices}
drug_mins = {drug: float('inf') for drug in drug_indices}
drug_maxs = {drug: float('-inf') for drug in drug_indices}

total_sum = 0.0
total_count = 0

row_count = 0
filtered_count = 0

with open(csv_path, 'r') as f:
    reader = csv.reader(f)
    next(reader)  # Skip header
    
    for row in reader:
        row_count += 1
        
        # Check year filter if available
        if year_idx is not None:
            try:
                year = float(row[year_idx])
                if year < 2023 or year > 2025:
                    continue
            except (ValueError, IndexError):
                continue
        
        filtered_count += 1
        
        # Process drug columns
        for drug, idx in drug_indices.items():
            try:
                value = float(row[idx])
                drug_sums[drug] += value
                drug_counts[drug] += 1
                drug_mins[drug] = min(drug_mins[drug], value)
                drug_maxs[drug] = max(drug_maxs[drug], value)
            except (ValueError, IndexError):
                pass
        
        # Process total drug count
        if total_drug_idx is not None:
            try:
                value = float(row[total_drug_idx])
                total_sum += value
                total_count += 1
            except (ValueError, IndexError):
                pass
        
        # Progress indicator
        if row_count % 100000 == 0:
            print(f"  Processed {row_count:,} rows...")

print(f"\nTotal rows: {row_count:,}")
print(f"Filtered rows (2023-2025): {filtered_count:,}")

# Calculate results
print("\n" + "="*70)
print("INDIVIDUAL DRUG USAGE ANALYSIS")
print("="*70)

total_drug_usage = total_sum / total_count if total_count > 0 else None
if total_drug_usage:
    print(f"\nTotal drug usage (currently_taking_drug_count): {total_drug_usage:.2f}")
else:
    print("\nWarning: 'currently_taking_drug_count' not found or no data")

# Calculate means and create results
results = []
for drug in drug_indices:
    if drug_counts[drug] > 0:
        mean_usage = drug_sums[drug] / drug_counts[drug]
        min_usage = drug_mins[drug]
        max_usage = drug_maxs[drug]
    else:
        mean_usage = min_usage = max_usage = 0
    
    drug_name = drug.replace('_currently_on_drug', '')
    
    result = {
        'drug': drug_name,
        'mean': mean_usage,
        'max': max_usage,
        'min': min_usage
    }
    
    if total_drug_usage:
        pct_of_total = (mean_usage / total_drug_usage) * 100
        result['pct_of_total'] = pct_of_total
    
    results.append(result)

# Sort by mean usage descending
results.sort(key=lambda x: x['mean'], reverse=True)

print("\n")
print(f"{'Drug':<20} {'Mean Usage':<12} {'% of Total':<12} {'Min':<10} {'Max':<10}")
print("-" * 70)

for r in results:
    drug_display = r['drug']
    mean_display = f"{r['mean']:.2f}"
    pct_display = f"{r['pct_of_total']:.4f}%" if 'pct_of_total' in r else "N/A"
    min_display = f"{r['min']:.2f}"
    max_display = f"{r['max']:.2f}"
    
    print(f"{drug_display:<20} {mean_display:<12} {pct_display:<12} {min_display:<10} {max_display:<10}")

# Calculate total for Polymyxins & Others category
if results:
    total_category = sum(r['mean'] for r in results)
    if total_drug_usage:
        category_pct = (total_category / total_drug_usage) * 100
        print("\n" + "-" * 70)
        print(f"{'TOTAL CATEGORY':<20} {total_category:.2f} {'':>5} {category_pct:.4f}%")
        print("="*70)

print("\n\nTOP CONTRIBUTORS:")
print("-" * 70)
for i, r in enumerate(results[:3], 1):
    if 'pct_of_total' in r:
        print(f"{i}. {r['drug']}: {r['mean']:.2f} patients ({r['pct_of_total']:.4f}% of total)")
    else:
        print(f"{i}. {r['drug']}: {r['mean']:.2f} patients")

print("\n\nINTERPRETATION:")
print("-" * 70)
if results and 'pct_of_total' in results[0]:
    top_drug = results[0]
    if top_drug['pct_of_total'] > 5:
        print(f"⚠️  {top_drug['drug'].upper()} is the primary driver of")
        print(f"   'Polymyxins & Others' category usage at {top_drug['pct_of_total']:.2f}%")
    print(f"\nTo address the 8.83% vs 0.4% target discrepancy:")
    print(f"Focus calibration efforts on reducing {results[0]['drug']} usage")
    if len(results) > 1:
        print(f"Secondary focus: {results[1]['drug']} at {results[1]['pct_of_total']:.4f}%")
