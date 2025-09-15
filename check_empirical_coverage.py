import pandas as pd

# Load the updated empirical resistance data
df = pd.read_csv('calibration_resistance_empirical.csv')

print("📊 EMPIRICAL RESISTANCE DATA COVERAGE")
print("=" * 50)
print(f"Total records: {len(df):,}")
print(f"Years covered: {df['year'].min()}-{df['year'].max()}")
print(f"Drugs: {sorted(df['drug'].unique())}")
print(f"Bacteria: {sorted(df['bacteria'].unique())}")
print(f"Drug-bacteria combinations: {len(df.groupby(['drug', 'bacteria']))}")

# Check source quality distribution
print(f"\nSource quality breakdown:")
quality_counts = df['source_quality'].value_counts()
for quality, count in quality_counts.items():
    print(f"  {quality}: {count:,} records")

# Check sample empirical vs synthetic data
empirical_sample = df[df['source_quality'] == 'who_glass_empirical']
print(f"\nEmpirical WHO GLASS data points: {len(empirical_sample)}")
if not empirical_sample.empty:
    print("Sample empirical records:")
    for _, row in empirical_sample.head(3).iterrows():
        print(f"  {row['year']}: {row['drug']} + {row['bacteria']} = {row['mean']:.1%}")

# Verify all combinations exist
expected_combinations = 8 * 3  # 8 drugs × 3 bacteria
years_covered = len(df['year'].unique())
expected_total = expected_combinations * years_covered

print(f"\nCoverage verification:")
print(f"  Expected combinations: {expected_combinations}")
print(f"  Years covered: {years_covered}")
print(f"  Expected total records: {expected_total:,}")
print(f"  Actual records: {len(df):,}")
print(f"  Coverage: {'✅ Complete' if len(df) >= expected_total else '❌ Incomplete'}")
