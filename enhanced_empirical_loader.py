#!/usr/bin/env python3
"""
Enhanced empirical data loading with integrated surveillance sources.

This module extends the original empirical data loader to integrate:
- WHO GLASS surveillance data  
- ECDC EARS-Net surveillance data
- Australian NNDSS surveillance data
- CDDEP ResistanceMap data

Usage:
    from enhanced_empirical_loader import load_integrated_empirical_data
    data = load_integrated_empirical_data()
"""

import pandas as pd
from pathlib import Path
import logging
from typing import Dict, Optional, List
import numpy as np

logger = logging.getLogger(__name__)

class IntegratedEmpiricalLoader:
    """Enhanced loader that integrates multiple surveillance data sources."""
    
    def __init__(self, data_dir: str = "data"):
        self.data_dir = Path(data_dir)
        self.who_dir = self.data_dir / "who"
        self.ecdc_dir = self.data_dir / "ecdc"
        self.australia_dir = self.data_dir / "australia"
        self.cddep_dir = self.data_dir / "cddep"
        
    def load_integrated_empirical_data(self) -> Dict[str, pd.DataFrame]:
        """
        Load and integrate empirical data from all surveillance sources.
        
        Returns:
            Dictionary with integrated empirical data for overlay on simulation plots
        """
        logger.info("Loading integrated empirical surveillance data...")
        
        integrated_data = {
            'resistance': None,
            'drug_usage': None,
            'incidence': None,
            'deaths': None,
            'drug_failure': None,
            'mic_values': None,
            'hospital_incidence': None
        }
        
        # Load new surveillance sources
        surveillance_data = self._load_surveillance_sources()
        
        # Load existing calibration data
        calibration_data = self._load_calibration_data()
        
        # Integrate surveillance data with calibration data
        integrated_data = self._integrate_data_sources(surveillance_data, calibration_data)
        
        logger.info("Integrated empirical data loading completed")
        return integrated_data
    
    def _load_surveillance_sources(self) -> Dict[str, Optional[pd.DataFrame]]:
        """Load data from the four new surveillance sources."""
        
        surveillance_sources = {
            'who_glass': self.who_dir / "glass_amr_surveillance.csv",
            'ecdc_ears_net': self.ecdc_dir / "ears_net_surveillance.csv", 
            'australia_nndss': self.australia_dir / "nndss_surveillance.csv",
            'cddep_resistancemap': self.cddep_dir / "resistancemap_surveillance.csv"
        }
        
        surveillance_data = {}
        
        for source_name, file_path in surveillance_sources.items():
            if file_path.exists():
                try:
                    df = pd.read_csv(file_path)
                    surveillance_data[source_name] = df
                    logger.info(f"   Loaded {len(df):,} records from {source_name}")
                except Exception as e:
                    logger.warning(f"   Failed to load {source_name}: {e}")
                    surveillance_data[source_name] = None
            else:
                logger.warning(f"   Surveillance file not found: {file_path}")
                surveillance_data[source_name] = None
                
        return surveillance_data
    
    def _load_calibration_data(self) -> Dict[str, Optional[pd.DataFrame]]:
        """Load existing calibration data files."""
        
        calibration_files = {
            'resistance': 'calibration_resistance_empirical.csv',
            'drug_usage': 'calibration_drug_usage_empirical.csv',
            'incidence': 'calibration_infection_incidence_empirical.csv',
            'deaths': 'calibration_deaths_empirical.csv',
            'drug_failure': 'calibration_drug_failure_empirical.csv',
            'mic_values': 'calibration_mic_empirical.csv',
            'hospital_incidence': 'calibration_hospital_incidence_empirical.csv'
        }
        
        calibration_data = {}
        
        for data_type, filename in calibration_files.items():
            file_path = Path(filename)
            if file_path.exists():
                try:
                    df = pd.read_csv(file_path)
                    calibration_data[data_type] = df
                    logger.info(f"   Loaded existing calibration data: {data_type} ({len(df):,} records)")
                except Exception as e:
                    logger.warning(f"   Failed to load calibration {data_type}: {e}")
                    calibration_data[data_type] = None
            else:
                logger.info(f"   Calibration file not found: {filename}")
                calibration_data[data_type] = None
                
        return calibration_data
    
    def _integrate_data_sources(self, surveillance_data: Dict, calibration_data: Dict) -> Dict[str, pd.DataFrame]:
        """Integrate surveillance and calibration data sources."""
        
        integrated = {}
        
        # Integrate resistance data
        integrated['resistance'] = self._integrate_resistance_data(surveillance_data, calibration_data)
        
        # Integrate drug usage data 
        integrated['drug_usage'] = self._integrate_drug_usage_data(surveillance_data, calibration_data)
        
        # Pass through other calibration data (enhanced later)
        for data_type in ['incidence', 'deaths', 'drug_failure', 'mic_values', 'hospital_incidence']:
            integrated[data_type] = calibration_data.get(data_type)
        
        return integrated
    
    def _integrate_resistance_data(self, surveillance_data: Dict, calibration_data: Dict) -> pd.DataFrame:
        """Integrate resistance surveillance data with existing calibration data."""
        
        # Start with existing calibration data
        base_resistance_data = calibration_data.get('resistance')
        
        if base_resistance_data is not None:
            integrated_resistance = base_resistance_data.copy()
        else:
            # Create empty DataFrame with proper structure
            integrated_resistance = pd.DataFrame(columns=[
                'year', 'drug', 'bacteria', 'mean', 'std', 'p5', 'p25', 'p50', 'p75', 'p95', 
                'units', 'source_quality', 'notes'
            ])
        
        # Process each surveillance source
        surveillance_sources = ['who_glass', 'ecdc_ears_net', 'australia_nndss', 'cddep_resistancemap']
        
        for source in surveillance_sources:
            source_data = surveillance_data.get(source)
            if source_data is not None:
                processed_data = self._process_resistance_surveillance_data(source_data, source)
                if processed_data is not None and len(processed_data) > 0:
                    integrated_resistance = pd.concat([integrated_resistance, processed_data], ignore_index=True)
        
        logger.info(f"   Integrated resistance data: {len(integrated_resistance):,} total records")
        return integrated_resistance
    
    def _process_resistance_surveillance_data(self, source_data: pd.DataFrame, source_name: str) -> Optional[pd.DataFrame]:
        """Process surveillance data into calibration format."""
        
        try:
            processed_records = []
            
            # Group by year, pathogen, antibiotic
            if source_name == 'who_glass':
                group_cols = ['year', 'pathogen', 'antibiotic']
            elif source_name == 'ecdc_ears_net':
                group_cols = ['year', 'pathogen', 'antibiotic']
            elif source_name == 'australia_nndss':
                group_cols = ['year', 'pathogen', 'antibiotic'] 
            elif source_name == 'cddep_resistancemap':
                group_cols = ['year', 'pathogen', 'antibiotic']
            else:
                logger.warning(f"Unknown surveillance source: {source_name}")
                return None
            
            # Process each group
            for group_key, group_data in source_data.groupby(group_cols):
                year, pathogen, antibiotic = group_key
                
                # Normalize names for matching
                normalized_bacteria = self._normalize_bacteria_name(pathogen)
                normalized_drug = self._normalize_drug_name(antibiotic)
                
                # Calculate statistics from resistance percentages
                resistance_values = group_data['resistance_percentage'].values / 100.0  # Convert to proportion
                
                if len(resistance_values) > 0:
                    mean_resistance = np.mean(resistance_values)
                    std_resistance = np.std(resistance_values) if len(resistance_values) > 1 else 0.0
                    
                    # Calculate percentiles
                    p5 = np.percentile(resistance_values, 5) if len(resistance_values) >= 2 else mean_resistance * 0.5
                    p25 = np.percentile(resistance_values, 25) if len(resistance_values) >= 2 else mean_resistance * 0.75
                    p50 = np.percentile(resistance_values, 50) if len(resistance_values) >= 2 else mean_resistance
                    p75 = np.percentile(resistance_values, 75) if len(resistance_values) >= 2 else mean_resistance * 1.25
                    p95 = np.percentile(resistance_values, 95) if len(resistance_values) >= 2 else mean_resistance * 1.5
                    
                    processed_records.append({
                        'year': year,
                        'drug': normalized_drug,
                        'bacteria': normalized_bacteria,
                        'mean': mean_resistance,
                        'std': std_resistance,
                        'p5': p5,
                        'p25': p25,
                        'p50': p50,
                        'p75': p75,
                        'p95': p95,
                        'units': 'proportion',
                        'source_quality': self._get_source_quality(source_name),
                        'notes': f'{source_name}_surveillance'
                    })
            
            if processed_records:
                return pd.DataFrame(processed_records)
            
        except Exception as e:
            logger.error(f"Error processing {source_name} resistance data: {e}")
            
        return None
    
    def _normalize_bacteria_name(self, bacteria_name: str) -> str:
        """Normalize bacteria names for consistent matching."""
        
        name = bacteria_name.lower().strip()
        
        # Common normalizations
        normalizations = {
            'escherichia coli': 'escherichia coli',
            'e. coli': 'escherichia coli',
            'e.coli': 'escherichia coli',
            'staphylococcus aureus': 'staphylococcus aureus', 
            's. aureus': 'staphylococcus aureus',
            'klebsiella pneumoniae': 'klebsiella pneumoniae',
            'k. pneumoniae': 'klebsiella pneumoniae',
            'pseudomonas aeruginosa': 'pseudomonas aeruginosa',
            'p. aeruginosa': 'pseudomonas aeruginosa',
            'enterococcus faecium': 'enterococcus faecium',
            'enterococcus faecalis': 'enterococcus faecalis',
            'enterococcus species': 'enterococcus faecium',  # Default to faecium
            'acinetobacter species': 'acinetobacter baumannii',
            'acinetobacter spp.': 'acinetobacter baumannii',
            'streptococcus pneumoniae': 'streptococcus pneumoniae',
            's. pneumoniae': 'streptococcus pneumoniae'
        }
        
        return normalizations.get(name, name)
    
    def _normalize_drug_name(self, drug_name: str) -> str:
        """Normalize drug names for consistent matching."""
        
        name = drug_name.lower().strip()
        
        # Common normalizations
        normalizations = {
            'ciprofloxacin': 'ciprofloxacin',
            'fluoroquinolones': 'ciprofloxacin',  # Use ciprofloxacin as representative
            'methicillin': 'methicillin',
            'ceftriaxone': 'ceftriaxone', 
            'ceftazidime': 'ceftazidime',
            '3rd-gen cephalosporins': 'ceftriaxone',  # Use ceftriaxone as representative
            'carbapenem': 'meropenem',  # Use meropenem as representative carbapenem
            'vancomycin': 'vancomycin',
            'ampicillin': 'ampicillin',
            'gentamicin': 'gentamicin',
            'colistin': 'colistin',
            'penicillin': 'penicillin',
            'erythromycin': 'erythromycin',
            'levofloxacin': 'levofloxacin',
            'tetracycline': 'tetracycline',
            'various': 'mixed'  # For surveillance data with multiple antibiotics
        }
        
        return normalizations.get(name, name)
    
    def _get_source_quality(self, source_name: str) -> str:
        """Get source quality rating."""
        
        quality_ratings = {
            'who_glass': 'high_quality_surveillance',
            'ecdc_ears_net': 'high_quality_surveillance', 
            'australia_nndss': 'medium_quality_surveillance',
            'cddep_resistancemap': 'aggregated_surveillance'
        }
        
        return quality_ratings.get(source_name, 'unknown_quality')
    
    def _integrate_drug_usage_data(self, surveillance_data: Dict, calibration_data: Dict) -> pd.DataFrame:
        """Integrate drug usage data from surveillance sources."""
        
        # Start with existing calibration data
        base_usage_data = calibration_data.get('drug_usage')
        
        if base_usage_data is not None:
            integrated_usage = base_usage_data.copy()
        else:
            integrated_usage = pd.DataFrame(columns=[
                'year', 'drug', 'bacteria', 'mean', 'std', 'p5', 'p25', 'p50', 'p75', 'p95',
                'units', 'source_quality', 'notes'
            ])
        
        # For now, pass through existing calibration data
        # Future enhancement: integrate AMU data from WHO GLASS and ECDC
        
        logger.info(f"   Drug usage data: {len(integrated_usage):,} records")
        return integrated_usage

def load_integrated_empirical_data() -> Dict[str, pd.DataFrame]:
    """
    Main function to load integrated empirical data.
    
    This function can be used as a drop-in replacement for the original
    load_empirical_calibration_data() function.
    """
    loader = IntegratedEmpiricalLoader()
    return loader.load_integrated_empirical_data()

# Compatibility function for existing code
def load_empirical_calibration_data():
    """
    Enhanced version of original function with integrated surveillance data.
    """
    return load_integrated_empirical_data()

if __name__ == "__main__":
    # Test the integrated loader
    data = load_integrated_empirical_data()
    
    print("Integrated Empirical Data Summary:")
    for data_type, df in data.items():
        if df is not None:
            print(f"  {data_type}: {len(df):,} records")
        else:
            print(f"  {data_type}: No data")