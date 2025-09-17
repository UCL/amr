#!/usr/bin/env python3
"""
Empirical Data Parsers for AMR Calibration
Handles real data formats from ECDC, WHO GLASS, IQVIA, CDC, and national sources
"""

import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
from typing import Dict, List, Tuple, Optional, Union
import requests
import json
import os
from pathlib import Path
from enhanced_empirical_data_config import EnhancedDataSourceConfig

class ECDCParser:
    """Parser for ECDC surveillance and consumption data"""
    
    def __init__(self, config: EnhancedDataSourceConfig):
        self.config = config
        
    def parse_resistance_data(self, file_path: str) -> pd.DataFrame:
        """
        Parse ECDC resistance surveillance data
        Expected format: Country, Year, Bacteria, Antibiotic, Resistance_percentage
        """
        try:
            df = pd.read_csv(file_path)
            
            # Standardize column names
            column_mapping = {
                'Country': 'country',
                'Year': 'year', 
                'Bacteria': 'bacteria',
                'Antibiotic': 'drug',
                'Resistance_percentage': 'resistance_rate',
                'Number_tested': 'sample_size'
            }
            
            df = df.rename(columns=column_mapping)
            
            # Convert resistance percentage to proportion
            if 'resistance_rate' in df.columns:
                df['resistance_rate'] = df['resistance_rate'] / 100.0
            
            # Standardize bacteria names
            bacteria_mapping = {
                'Escherichia coli': 'escherichia_coli',
                'E. coli': 'escherichia_coli',
                'Staphylococcus aureus': 'staphylococcus_aureus',
                'S. aureus': 'staphylococcus_aureus',
                'MRSA': 'staphylococcus_aureus',  # Methicillin-resistant S. aureus
                'Mycobacterium tuberculosis': 'mycobacterium_tuberculosis'
            }
            
            df['bacteria'] = df['bacteria'].map(bacteria_mapping).fillna(df['bacteria'])
            
            # Standardize drug names
            drug_mapping = {
                'Ciprofloxacin': 'ciprofloxacin',
                'Penicillin': 'penicillin',
                'Methicillin': 'methicillin',
                'Erythromycin': 'erythromycin',
                'Tetracycline': 'tetracycline',
                'Rifampicin': 'rifampicin',
                'Streptomycin': 'streptomycin'
            }
            
            df['drug'] = df['drug'].map(drug_mapping).fillna(df['drug'])
            
            return df
            
        except Exception as e:
            print(f"Error parsing ECDC resistance data: {e}")
            return pd.DataFrame()
    
    def parse_consumption_data(self, file_path: str) -> pd.DataFrame:
        """
        Parse ECDC antibiotic consumption data  
        Expected format: Country, Year, ATC_code, DDD_per_1000_inhabitants_per_day
        """
        try:
            df = pd.read_csv(file_path)
            
            # Map ATC codes to drug names
            atc_to_drug = {
                'J01CA01': 'penicillin',      # Ampicillin
                'J01CA04': 'penicillin',      # Amoxicillin  
                'J01MA02': 'ciprofloxacin',   # Ciprofloxacin
                'J01AA02': 'tetracycline',    # Doxycycline
                'J01FA01': 'erythromycin',    # Erythromycin
                'J01XD01': 'rifampicin',      # Rifampicin
                'J01GA01': 'streptomycin'     # Streptomycin
            }
            
            df['drug'] = df['ATC_code'].map(atc_to_drug)
            df = df.dropna(subset=['drug'])  # Remove unmapped ATC codes
            
            # Convert DDD to courses per 100k (approximate)
            # Assumption: 1 DDD/1000 inhabitants/day ≈ 365 courses/1000 inhabitants/year
            df['courses_per_100k'] = df['DDD_per_1000_inhabitants_per_day'] * 36.5
            
            return df
            
        except Exception as e:
            print(f"Error parsing ECDC consumption data: {e}")
            return pd.DataFrame()

class WHOGLASSParser:
    """Parser for WHO GLASS resistance data"""
    
    def __init__(self, config: EnhancedDataSourceConfig):
        self.config = config
        
    def parse_glass_excel(self, file_path: str) -> pd.DataFrame:
        """
        Parse WHO GLASS Excel surveillance reports
        Multiple sheets with different bacteria/drug combinations
        """
        try:
            # WHO GLASS typically has multiple sheets
            excel_file = pd.ExcelFile(file_path)
            
            combined_df = pd.DataFrame()
            
            for sheet_name in excel_file.sheet_names:
                if 'resistance' in sheet_name.lower() or 'amr' in sheet_name.lower():
                    sheet_df = pd.read_excel(file_path, sheet_name=sheet_name)
                    
                    # Look for standard GLASS columns
                    if all(col in sheet_df.columns for col in ['Country', 'Year', 'Bacteria', 'Antibiotic', '%R']):
                        sheet_df = sheet_df.rename(columns={
                            'Country': 'country',
                            'Year': 'year',
                            'Bacteria': 'bacteria', 
                            'Antibiotic': 'drug',
                            '%R': 'resistance_rate',
                            'No tested': 'sample_size'
                        })
                        
                        # Convert percentage to proportion
                        sheet_df['resistance_rate'] = sheet_df['resistance_rate'] / 100.0
                        
                        combined_df = pd.concat([combined_df, sheet_df], ignore_index=True)
            
            return combined_df
            
        except Exception as e:
            print(f"Error parsing WHO GLASS data: {e}")
            return pd.DataFrame()

class IQVIAParser:
    """Parser for IQVIA pharmaceutical sales data"""
    
    def __init__(self, config: EnhancedDataSourceConfig):
        self.config = config
        
    def parse_sales_data(self, file_path: str) -> pd.DataFrame:
        """
        Parse IQVIA MIDAS sales data
        Expected format: Country, Year, Product, Molecule, Units, Value
        """
        try:
            df = pd.read_csv(file_path)
            
            # Map IQVIA molecules to standard drug names
            molecule_mapping = {
                'AMOXICILLIN': 'penicillin',
                'AMPICILLIN': 'penicillin', 
                'CIPROFLOXACIN': 'ciprofloxacin',
                'DOXYCYCLINE': 'tetracycline',
                'TETRACYCLINE': 'tetracycline',
                'ERYTHROMYCIN': 'erythromycin',
                'AZITHROMYCIN': 'erythromycin',  # Macrolide class
                'RIFAMPICIN': 'rifampicin',
                'STREPTOMYCIN': 'streptomycin',
                'METHICILLIN': 'methicillin'
            }
            
            df['drug'] = df['Molecule'].str.upper().map(molecule_mapping)
            df = df.dropna(subset=['drug'])
            
            # Convert units to standardized measure
            # IQVIA typically provides "Standard Units" which can be converted
            df['standardized_units'] = df['Units']  # Placeholder for conversion logic
            
            return df
            
        except Exception as e:
            print(f"Error parsing IQVIA sales data: {e}")
            return pd.DataFrame()

class CDCParser:
    """Parser for CDC surveillance data"""
    
    def __init__(self, config: EnhancedDataSourceConfig):
        self.config = config
        
    def parse_ar_threats_report(self, file_path: str) -> pd.DataFrame:
        """
        Parse CDC AR Threats Report data
        May require manual extraction from PDF to CSV
        """
        try:
            df = pd.read_csv(file_path)
            
            # Expected columns: Pathogen, Cases_2019, Deaths_2019, Resistance_mechanism
            pathogen_mapping = {
                'Carbapenem-resistant Enterobacteriaceae': 'escherichia_coli',
                'Methicillin-resistant Staphylococcus aureus': 'staphylococcus_aureus',
                'Drug-resistant tuberculosis': 'mycobacterium_tuberculosis'
            }
            
            df['bacteria'] = df['Pathogen'].map(pathogen_mapping)
            
            return df
            
        except Exception as e:
            print(f"Error parsing CDC AR threats data: {e}")
            return pd.DataFrame()

class NationalStatsParser:
    """Parser for national mortality statistics"""
    
    def __init__(self, config: EnhancedDataSourceConfig):
        self.config = config
        
    def parse_mortality_data(self, file_path: str, country: str) -> pd.DataFrame:
        """
        Parse national mortality statistics
        Format varies by country but typically includes ICD-10 codes
        """
        try:
            df = pd.read_csv(file_path)
            
            # Map ICD-10 codes to bacterial causes
            icd10_mapping = {
                'A41': 'sepsis_unspecified',      # Sepsis
                'A49.9': 'bacterial_infection',   # Bacterial infection, unspecified
                'B95.0': 'staphylococcus_aureus', # Streptococcus, group A
                'B95.1': 'staphylococcus_aureus', # Streptococcus, group B
                'B96.2': 'escherichia_coli',      # Escherichia coli
                'A15': 'mycobacterium_tuberculosis', # Respiratory tuberculosis
                'A16': 'mycobacterium_tuberculosis'  # Tuberculosis of other organs
            }
            
            # Process based on country-specific format
            if country.lower() == 'united_states':
                # US format typically has ICD_code, Deaths, Population
                df['bacteria'] = df['ICD_code'].map(icd10_mapping)
                df['death_rate_per_100k'] = (df['Deaths'] / df['Population']) * 100000
                
            elif country.lower() in ['united_kingdom', 'germany']:
                # European format may have different column names
                df = self._harmonize_european_format(df)
            
            df['country'] = country
            return df.dropna(subset=['bacteria'])
            
        except Exception as e:
            print(f"Error parsing {country} mortality data: {e}")
            return pd.DataFrame()
    
    def _harmonize_european_format(self, df: pd.DataFrame) -> pd.DataFrame:
        """Harmonize European mortality data formats"""
        # Implementation depends on specific national formats
        # This is a placeholder for country-specific processing
        return df

class EmpiricalDataLoader:
    """Main class to coordinate all data source parsers"""
    
    def __init__(self):
        self.config = EnhancedDataSourceConfig()
        self.ecdc_parser = ECDCParser(self.config)
        self.who_parser = WHOGLASSParser(self.config)
        self.iqvia_parser = IQVIAParser(self.config)
        self.cdc_parser = CDCParser(self.config)
        self.stats_parser = NationalStatsParser(self.config)
        
    def check_data_availability(self) -> Dict[str, bool]:
        """Check which empirical data files are available"""
        availability = {}
        
        # Check ECDC files
        ecdc_files = self.config.ECDC_CONFIG['file_paths']
        availability['ecdc_resistance'] = os.path.exists(ecdc_files['resistance_csv'])
        availability['ecdc_consumption'] = os.path.exists(ecdc_files['consumption_csv'])
        
        # Check WHO GLASS files
        glass_files = self.config.WHO_GLASS_CONFIG['file_paths']
        availability['who_glass_2022'] = os.path.exists(glass_files['glass_data_2022'])
        availability['who_glass_2021'] = os.path.exists(glass_files['glass_data_2021'])
        
        # Check IQVIA files
        iqvia_files = self.config.IQVIA_CONFIG['file_paths']
        availability['iqvia_sales'] = os.path.exists(iqvia_files['sales_data_2023'])
        
        # Check CDC files
        cdc_files = self.config.CDC_CONFIG['file_paths']
        availability['cdc_ar_threats'] = os.path.exists(cdc_files['ar_threats_2019'])
        
        return availability
    
    def load_all_available_data(self) -> Dict[str, pd.DataFrame]:
        """Load all available empirical data sources"""
        print("🔍 Checking data availability...")
        availability = self.check_data_availability()
        
        loaded_data = {}
        
        # Load ECDC data if available
        if availability['ecdc_resistance']:
            print("📊 Loading ECDC resistance data...")
            ecdc_resistance = self.ecdc_parser.parse_resistance_data(
                self.config.ECDC_CONFIG['file_paths']['resistance_csv']
            )
            loaded_data['ecdc_resistance'] = ecdc_resistance
            
        if availability['ecdc_consumption']:
            print("💊 Loading ECDC consumption data...")
            ecdc_consumption = self.ecdc_parser.parse_consumption_data(
                self.config.ECDC_CONFIG['file_paths']['consumption_csv']
            )
            loaded_data['ecdc_consumption'] = ecdc_consumption
        
        # Load WHO GLASS data if available
        for year_file, available in [('who_glass_2022', availability['who_glass_2022']), 
                                   ('who_glass_2021', availability['who_glass_2021'])]:
            if available:
                print(f"🌍 Loading WHO GLASS data: {year_file}...")
                glass_data = self.who_parser.parse_glass_excel(
                    self.config.WHO_GLASS_CONFIG['file_paths'][year_file.replace('who_glass_', 'glass_data_')]
                )
                loaded_data[year_file] = glass_data
        
        # Load IQVIA data if available
        if availability['iqvia_sales']:
            print("💰 Loading IQVIA sales data...")
            iqvia_data = self.iqvia_parser.parse_sales_data(
                self.config.IQVIA_CONFIG['file_paths']['sales_data_2023']
            )
            loaded_data['iqvia_sales'] = iqvia_data
        
        # Load CDC data if available
        if availability['cdc_ar_threats']:
            print("🇺🇸 Loading CDC AR threats data...")
            cdc_data = self.cdc_parser.parse_ar_threats_report(
                self.config.CDC_CONFIG['file_paths']['ar_threats_2019']
            )
            loaded_data['cdc_ar_threats'] = cdc_data
        
        # Report what was loaded
        print(f"\n✅ Successfully loaded {len(loaded_data)} data sources")
        for source, df in loaded_data.items():
            print(f"   • {source}: {len(df)} records")
        
        # Show setup instructions if no data was loaded
        if not loaded_data:
            print("\n⚠️  No empirical data files found.")
            print("📋 Please run the enhanced data collection strategy:")
            print("   python implementation_guide.py")
            print("   python enhanced_empirical_data_collection.py")
        
        return loaded_data

def main():
    """Test the empirical data loaders"""
    loader = EmpiricalDataLoader()
    data = loader.load_all_available_data()
    
    if data:
        print("\n📈 Data Summary:")
        for source, df in data.items():
            print(f"  {source}: {df.shape[0]} rows × {df.shape[1]} columns")
    else:
        print("\n📁 To use empirical data, please:")
        print("   1. Create ./data/ directory")
        print("   2. Download data files from sources")
        print("   3. Update file paths in empirical_data_config.py")

if __name__ == "__main__":
    main()
