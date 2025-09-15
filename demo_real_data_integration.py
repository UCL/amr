#!/usr/bin/env python3
"""
Real Data Integration Demo
Shows how to integrate actual ECDC, WHO GLASS, IQVIA, and CDC data files
"""

import pandas as pd
import numpy as np
from empirical_data_parsers import EmpiricalDataLoader
import os

def create_sample_real_data_files():
    """
    Create sample data files in the expected formats from real sources
    This demonstrates the file structures you would download from actual sources
    """
    # Create data directory
    os.makedirs('./data', exist_ok=True)
    
    print("📁 Creating sample real data files...")
    
    # 1. Sample ECDC resistance data format
    ecdc_resistance = pd.DataFrame({
        'Country': ['Germany', 'France', 'Spain', 'Germany', 'France'],
        'Year': [2022, 2022, 2022, 2021, 2021],
        'Bacteria': ['Escherichia coli', 'Escherichia coli', 'E. coli', 'S. aureus', 'MRSA'],
        'Antibiotic': ['Ciprofloxacin', 'Ciprofloxacin', 'Tetracycline', 'Methicillin', 'Methicillin'],
        'Resistance_percentage': [25.3, 28.7, 35.2, 8.1, 12.4],
        'Number_tested': [1205, 987, 756, 2134, 1876]
    })
    ecdc_resistance.to_csv('./data/ecdc_resistance_2023.csv', index=False)
    print("   ✓ ECDC resistance data: ./data/ecdc_resistance_2023.csv")
    
    # 2. Sample ECDC consumption data format
    ecdc_consumption = pd.DataFrame({
        'Country': ['Germany', 'Germany', 'France', 'France', 'Spain'],
        'Year': [2022, 2022, 2022, 2022, 2022],
        'ATC_code': ['J01MA02', 'J01CA04', 'J01MA02', 'J01AA02', 'J01CA04'],
        'DDD_per_1000_inhabitants_per_day': [1.8, 4.2, 2.1, 1.5, 3.8]
    })
    ecdc_consumption.to_csv('./data/ecdc_consumption_2023.csv', index=False)
    print("   ✓ ECDC consumption data: ./data/ecdc_consumption_2023.csv")
    
    # 3. Sample IQVIA sales data format
    iqvia_sales = pd.DataFrame({
        'Country': ['United States', 'United States', 'Germany', 'Germany'],
        'Year': [2023, 2023, 2023, 2023],
        'Product': ['CIPRO', 'AMOXIL', 'CIPROBAY', 'AMOXI'],
        'Molecule': ['CIPROFLOXACIN', 'AMOXICILLIN', 'CIPROFLOXACIN', 'AMOXICILLIN'],
        'Units': [28700000, 45200000, 8400000, 12800000],  # Standard units
        'Value': [156000000, 89000000, 42000000, 28000000]  # USD
    })
    iqvia_sales.to_csv('./data/iqvia_sales_2023.csv', index=False)
    print("   ✓ IQVIA sales data: ./data/iqvia_sales_2023.csv")
    
    # 4. Sample CDC AR threats data format
    cdc_threats = pd.DataFrame({
        'Pathogen': ['Carbapenem-resistant Enterobacteriaceae', 'Methicillin-resistant Staphylococcus aureus'],
        'Cases_2019': [13100, 120000],
        'Deaths_2019': [1100, 9700],
        'Resistance_mechanism': ['Carbapenemase', 'mecA gene']
    })
    cdc_threats.to_csv('./data/cdc_ar_threats_2019.csv', index=False)
    print("   ✓ CDC AR threats data: ./data/cdc_ar_threats_2019.csv")
    
    # 5. Sample national mortality data (USA format)
    usa_mortality = pd.DataFrame({
        'ICD_code': ['A41.9', 'B96.2', 'A15.0', 'B95.0'],
        'Deaths': [25600, 8900, 542, 15200],
        'Population': [331000000, 331000000, 331000000, 331000000],
        'Year': [2022, 2022, 2022, 2022]
    })
    usa_mortality.to_csv('./data/usa_mortality_icd10.csv', index=False)
    print("   ✓ USA mortality data: ./data/usa_mortality_icd10.csv")

def demonstrate_real_data_integration():
    """
    Demonstrate loading and processing real data files
    """
    print("\n🔬 REAL DATA INTEGRATION DEMONSTRATION")
    print("=" * 60)
    
    # Create sample files first
    create_sample_real_data_files()
    
    # Load the data using the empirical data loader
    loader = EmpiricalDataLoader()
    
    print("\n📊 Loading empirical data sources...")
    loaded_data = loader.load_all_available_data()
    
    if loaded_data:
        print(f"\n✅ Successfully loaded {len(loaded_data)} real data sources:")
        
        for source_name, df in loaded_data.items():
            print(f"\n📈 {source_name.upper()}:")
            print(f"   • Records: {len(df)}")
            print(f"   • Columns: {list(df.columns)}")
            print(f"   • Sample data:")
            print(df.head(2).to_string(index=False))
    
    # Show how to combine with the calibration system
    print("\n🔗 INTEGRATION WITH CALIBRATION SYSTEM")
    print("=" * 60)
    
    print("To integrate this real data with your simulation:")
    print("1. Download actual data files from sources")
    print("2. Place them in ./data/ directory with correct names")
    print("3. Run: python generate_empirical_calibration.py")
    print("4. Use generated *_empirical.csv files for model calibration")
    
    return loaded_data

def show_data_source_comparison():
    """
    Compare synthetic vs empirical calibration data
    """
    print("\n📊 SYNTHETIC vs EMPIRICAL DATA COMPARISON")
    print("=" * 60)
    
    # Load both types if available
    try:
        synthetic_resistance = pd.read_csv('calibration_resistance.csv')
        empirical_resistance = pd.read_csv('calibration_resistance_empirical.csv')
        
        print("SYNTHETIC DATA:")
        print(f"  Records: {len(synthetic_resistance):,}")
        print(f"  Source quality: {synthetic_resistance['source_quality'].value_counts().to_dict()}")
        
        print("\nEMPIRICAL DATA:")
        print(f"  Records: {len(empirical_resistance):,}")
        print(f"  Source quality: {empirical_resistance['source_quality'].value_counts().to_dict()}")
        
        # Show resistance rate comparison for ciprofloxacin in E. coli
        if not empirical_resistance.empty:
            cip_ecoli = empirical_resistance[
                (empirical_resistance['drug'] == 'ciprofloxacin') & 
                (empirical_resistance['bacteria'] == 'escherichia_coli')
            ]
            
            if not cip_ecoli.empty:
                print(f"\n🦠 Ciprofloxacin resistance in E. coli (empirical):")
                print(f"  WHO GLASS USA (2019): {cip_ecoli.iloc[0]['mean']:.1%}")
                print(f"  WHO GLASS India (2019): {cip_ecoli.iloc[1]['mean']:.1%} (if available)")
                print("  ➡️  This shows real-world variation that synthetic data approximates")
        
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        print("Run calibration generators first to create comparison files")

def main():
    """Main demonstration"""
    print("🌍 EMPIRICAL DATA INTEGRATION FOR AMR CALIBRATION")
    print("=" * 70)
    print("This demo shows how to integrate real data from:")
    print("  • ECDC: European antimicrobial resistance surveillance")
    print("  • WHO GLASS: Global resistance surveillance system")
    print("  • IQVIA: Pharmaceutical sales and usage data")
    print("  • CDC: US antimicrobial resistance surveillance")
    print("  • National statistics: Country-specific mortality data")
    
    # Run the demonstration
    loaded_data = demonstrate_real_data_integration()
    
    # Show comparison with synthetic data
    show_data_source_comparison()
    
    print("\n🎯 NEXT STEPS FOR REAL IMPLEMENTATION:")
    print("=" * 60)
    print("1. 📋 Follow setup instructions in empirical_data_config.py")
    print("2. 🔑 Obtain API keys/access credentials for data sources")
    print("3. 📥 Download actual surveillance data files")
    print("4. 🔧 Configure file paths in empirical_data_config.py")
    print("5. ▶️  Run generate_empirical_calibration.py")
    print("6. ✅ Use *_empirical.csv files for model validation")

if __name__ == "__main__":
    main()
