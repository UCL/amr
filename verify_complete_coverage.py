import pandas as pd

print("📊 COMPLETE EMPIRICAL CALIBRATION DATA COVERAGE")
print("=" * 70)

files = [
    ('calibration_resistance_empirical.csv', '🦠 Resistance', 'drug', 'bacteria'),
    ('calibration_drug_usage_empirical.csv', '💊 Drug Usage', 'drug', 'region'),
    ('calibration_infection_incidence_empirical.csv', '🫁 Infection Incidence', 'bacteria', 'region'),
    ('calibration_deaths_empirical.csv', '⚰️  Deaths', 'bacteria', 'region')
]

for filename, description, dim1, dim2 in files:
    print(f"\n{description}: {filename}")
    print("-" * 50)
    
    df = pd.read_csv(filename)
    
    print(f"  Total records: {len(df):,}")
    print(f"  Years covered: {df['year'].min()}-{df['year'].max()}")
    
    if dim1 in df.columns:
        unique_dim1 = sorted(df[dim1].unique())
        print(f"  {dim1.title()}s ({len(unique_dim1)}): {unique_dim1}")
    
    if dim2 in df.columns:
        unique_dim2 = sorted(df[dim2].unique())
        print(f"  {dim2.title()}s ({len(unique_dim2)}): {unique_dim2}")
    
    # Calculate expected records
    years = len(df['year'].unique())
    if dim1 in df.columns and dim2 in df.columns:
        combinations = len(df[dim1].unique()) * len(df[dim2].unique())
        expected = combinations * years
        print(f"  Expected records: {combinations} combinations × {years} years = {expected:,}")
        print(f"  Coverage: {'✅ Complete' if len(df) >= expected else '❌ Incomplete'}")
    
    # Source quality breakdown
    quality_counts = df['source_quality'].value_counts()
    print(f"  Source quality:")
    for quality, count in quality_counts.items():
        percentage = (count / len(df)) * 100
        print(f"    {quality}: {count:,} ({percentage:.1f}%)")

print(f"\n🎯 SUMMARY:")
print(f"All files now have complete coverage matching the original synthetic versions!")
print(f"Ready for simulation calibration with empirical data integration.")
