"""
Analyze individual drug usage within Polymyxins & Others category
"""
import csv
from collections import defaultdict

# Read the CSV file
csv_path = 'amr_simulation_output_analysis_outputs/simulation_summary_051228.csv'
print(f"Loading {csv_path}...")

# Load data using csv module
data = defaultdict(list)
with open(csv_path, 'r') as f:
    reader = csv.DictReader(f)
    headers = reader.fieldnames
    for row in reader:
        for key, value in row.items():
            try:
                data[key].append(float(value))
            except (ValueError, TypeError):
                data[key].append(value)

print(f"Loaded {len(next(iter(data.values())))} rows")
print(f"Found {len(data)} columns")

# Define the target drugs
target_drugs = [
    'colistin_currently_on_drug',
    'rifampicin_currently_on_drug',
    'chloramphenicol_currently_on_drug',
    'fusidic_a_currently_on_drug',
    'retapamulin_currently_on_drug',
    'quinu_dalfo_currently_on_drug'
]

# Check which columns exist
existing_drugs = [col for col in target_drugs if col in data]
missing_drugs = [col for col in target_drugs if col not in data]

print(f"\nFound {len(existing_drugs)} drug columns:")
for drug in existing_drugs:
    print(f"  - {drug}")

if missing_drugs:
    print(f"\nMissing {len(missing_drugs)} drug columns:")
    for drug in missing_drugs:
        print(f"  - {drug}")

# Filter to 2023-2025 time window
if 'year' in data:
    # Find indices where year is between 2023 and 2025
    year_data = data['year']
    indices = [i for i, year in enumerate(year_data) if isinstance(year, (int, float)) and 2023 <= year <= 2025]
    print(f"\nFiltered to years 2023-2025: {len(indices)} rows")
    
    # Create filtered data
    filtered_data = {key: [values[i] for i in indices] for key, values in data.items()}
else:
    filtered_data = data
    print("\nNo year column found, using all data")
    indices = list(range(len(next(iter(data.values())))))

# Calculate means for each drug
print("\n" + "="*70)
print("INDIVIDUAL DRUG USAGE ANALYSIS")
print("="*70)

# Get total drug count
if 'currently_taking_drug_count' in filtered_data:
    total_values = [v for v in filtered_data['currently_taking_drug_count'] if isinstance(v, (int, float))]
    total_drug_usage = sum(total_values) / len(total_values) if total_values else 0
    print(f"\nTotal drug usage (currently_taking_drug_count): {total_drug_usage:.2f}")
else:
    print("\nWarning: 'currently_taking_drug_count' column not found")
    total_drug_usage = None

# Analyze each drug
results = []
for drug_col in existing_drugs:
    values = [v for v in filtered_data[drug_col] if isinstance(v, (int, float))]
    
    if values:
        mean_usage = sum(values) / len(values)
        max_usage = max(values)
        min_usage = min(values)
    else:
        mean_usage = max_usage = min_usage = 0
    
    drug_name = drug_col.replace('_currently_on_drug', '')
    
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
