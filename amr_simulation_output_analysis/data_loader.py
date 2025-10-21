#!/usr/bin/env python3
"""
Data Loading and Caching for AMR Simulation Output Analysis

This module handles loading and caching of simulation data to eliminate 
repeated CSV reads that were causing performance issues in the original
analyze_simulation.py script.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any
import logging
from .config import DataConfig

logger = logging.getLogger(__name__)

class DataCache:
    """
    Singleton cache for simulation and empirical data.
    
    Eliminates repeated CSV reads by caching loaded data in memory.
    Provides methods to reload data when needed.
    """
    
    _instance: Optional['DataCache'] = None
    
    def __new__(cls) -> 'DataCache':
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance
    
    def __init__(self):
        if self._initialized:
            return
            
        self._simulation_data: Optional[pd.DataFrame] = None
        self._preprocessed_data: Optional[pd.DataFrame] = None
        self._empirical_data: Dict[str, Optional[pd.DataFrame]] = {}
        self._bacteria_list: Optional[list] = None
        self._drug_list: Optional[list] = None
        self._resistance_mechanisms: Optional[list] = None
        self._initialized = True
        
        logger.info("DataCache initialized")
    
    def get_simulation_data(self, csv_file: str = None, force_reload: bool = False) -> Optional[pd.DataFrame]:
        """
        Get cached simulation data, loading if necessary.
        
        Args:
            csv_file: Path to CSV file (uses default if None)
            force_reload: Force reload even if cached
            
        Returns:
            DataFrame with simulation data or None if loading failed
        """
        if self._simulation_data is None or force_reload:
            if csv_file is None:
                csv_file = str(DataConfig().simulation_file)
                
            self._simulation_data = load_simulation_data(csv_file)
            
            # Clear dependent cached data when simulation data reloads
            if self._simulation_data is not None:
                self._preprocessed_data = None
                self._bacteria_list = None
                self._drug_list = None 
                self._resistance_mechanisms = None
                logger.info(f"Simulation data loaded and cached: {len(self._simulation_data)} rows")
        
        return self._simulation_data
    
    def get_preprocessed_data(self, force_reload: bool = False) -> Optional[pd.DataFrame]:
        """
        Get cached preprocessed data, processing if necessary.
        
        Args:
            force_reload: Force reprocessing even if cached
            
        Returns:
            DataFrame with preprocessed simulation data or None if failed
        """
        if self._preprocessed_data is None or force_reload:
            sim_data = self.get_simulation_data()
            if sim_data is not None:
                self._preprocessed_data = preprocess_data(sim_data.copy())
                logger.info("Data preprocessing completed and cached")
        
        return self._preprocessed_data
    
    def get_bacteria_list(self, force_reload: bool = False) -> list:
        """Get cached bacteria list extracted from CSV headers."""
        if self._bacteria_list is None or force_reload:
            sim_data = self.get_simulation_data()
            if sim_data is not None:
                self._bacteria_list = extract_bacteria_list_from_csv(sim_data)
                logger.info(f"Extracted {len(self._bacteria_list)} bacteria from CSV headers")
        
        return self._bacteria_list or []
    
    def get_drug_list(self, force_reload: bool = False) -> list:
        """Get cached drug list extracted from CSV headers."""
        if self._drug_list is None or force_reload:
            sim_data = self.get_simulation_data()
            if sim_data is not None:
                self._drug_list = extract_drug_list_from_csv(sim_data)
                logger.info(f"Extracted {len(self._drug_list)} drugs from CSV headers")
        
        return self._drug_list or []
    
    def get_resistance_mechanisms(self, force_reload: bool = False) -> list:
        """Get cached resistance mechanisms extracted from CSV headers.""" 
        if self._resistance_mechanisms is None or force_reload:
            sim_data = self.get_simulation_data()
            if sim_data is not None:
                self._resistance_mechanisms = extract_resistance_mechanisms_from_csv(sim_data)
                logger.info(f"Extracted {len(self._resistance_mechanisms)} resistance mechanisms")
        
        return self._resistance_mechanisms or []
    
    def clear_cache(self):
        """Clear all cached data to free memory."""
        self._simulation_data = None
        self._preprocessed_data = None
        self._empirical_data.clear()
        self._bacteria_list = None
        self._drug_list = None
        self._resistance_mechanisms = None
        logger.info("DataCache cleared")

# Global cache instance
_cache = DataCache()

def get_cache() -> DataCache:
    """Get the global data cache instance."""
    return _cache

def load_simulation_data(csv_file: str) -> Optional[pd.DataFrame]:
    """
    Load simulation data from CSV file.
    
    Args:
        csv_file: Path to the simulation summary CSV file
        
    Returns:
        DataFrame with simulation data or None if loading failed
    """
    csv_path = Path(csv_file)
    
    if not csv_path.exists():
        logger.error(f"CSV file not found: {csv_file}")
        print(f"Error: {csv_file} not found. Run the Rust simulation first.")
        return None
    
    try:
        df = pd.read_csv(csv_file)
        logger.info(f"Loaded {len(df)} time steps from {csv_file}")
        print(f"Loaded {len(df)} time steps of simulation data")
        return df
        
    except Exception as e:
        logger.error(f"Error loading {csv_file}: {e}")
        print(f"Error loading {csv_file}: {e}")
        return None

def safe_divide(numerator, denominator, default=0):
    """Safe division avoiding division by zero."""
    return np.where(denominator > 0, numerator / denominator, default)

def preprocess_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add calculated columns and prepare data for analysis.
    
    Args:
        df: Raw simulation data DataFrame
        
    Returns:
        DataFrame with additional calculated columns
    """
    logger.info("Starting data preprocessing")
    
    # Age group proportions
    if 'num_age_0_5' in df.columns and 'total_population' in df.columns:
        df['prop_age_0_5'] = safe_divide(df['num_age_0_5'], df['total_population'])
        df['prop_age_6_14'] = safe_divide(df['num_age_6_14'], df['total_population'])
        df['prop_age_15_49'] = safe_divide(df['num_age_15_49'], df['total_population'])
        df['prop_age_50_79'] = safe_divide(df['num_age_50_79'], df['total_population'])
        df['prop_age_80plus'] = safe_divide(df['num_age_80plus'], df['total_population'])
        
    # Proportion of currently infected who are on drug
    if 'currently_infected_and_on_drug_count' in df.columns and 'total_currently_infected' in df.columns:
        df['infected_and_on_drug_proportion'] = safe_divide(
            df['currently_infected_and_on_drug_count'], 
            df['total_currently_infected']
        )
        
    # Calculate rolling past-year newly infected proportion
    if 'newly_infected_past_year' in df.columns and 'total_population' in df.columns:
        df['newly_infected_past_year_proportion'] = safe_divide(
            df['newly_infected_past_year'], 
            df['total_population']
        )
        
    # Calculate rolling past-year death proportions
    death_year_cols = [
        ('deaths_past_year', 'deaths_past_year_proportion'),
        ('deaths_background_past_year', 'deaths_background_past_year_proportion'),
        ('deaths_sepsis_past_year', 'deaths_sepsis_past_year_proportion'),
        ('deaths_drug_toxicity_past_year', 'deaths_drug_toxicity_past_year_proportion')
    ]
    
    for death_col, prop_col in death_year_cols:
        if death_col in df.columns and 'total_population' in df.columns:
            df[prop_col] = safe_divide(df[death_col], df['total_population'])

    # Convert time step to years
    df['time_in_years'] = df['time_step'] / 365
    
    # Calculate basic proportions
    df['infection_proportion'] = safe_divide(df['total_currently_infected'], df['total_population'])
    df['death_proportion'] = safe_divide(df['total_deaths'], df['total_population'])
    
    # Calculate resistance proportion among infected
    df['resistance_among_infected'] = safe_divide(df['total_with_resistance'], df['total_currently_infected'])
    
    # Calculate infection duration proportions
    df['infected_10_days_proportion'] = safe_divide(df['infected_10_days_count'], df['total_currently_infected'])
    df['infected_30_days_proportion'] = safe_divide(df['infected_30_days_count'], df['total_currently_infected'])
    
    # Calculate sepsis proportion among infected
    if 'number_with_sepsis' in df.columns:
        df['sepsis_among_infected_proportion'] = safe_divide(
            df['number_with_sepsis'], 
            df['total_currently_infected']
        )

    # Derive carrier vs non-carrier infection metrics for each bacteria
    carrier_suffix = '_infected_carrier_count'
    for carrier_col in [col for col in df.columns if col.endswith(carrier_suffix)]:
        slug = carrier_col[:-len(carrier_suffix)]
        non_carrier_col = f"{slug}_infected_non_carrier_count"
        res_carrier_col = f"{slug}_resistant_infected_carrier_count"
        res_non_carrier_col = f"{slug}_resistant_infected_non_carrier_count"

        if not all(col in df.columns for col in [non_carrier_col, res_carrier_col, res_non_carrier_col]):
            logger.debug(
                "Skipping derived carrier metrics for %s due to missing columns", slug
            )
            continue

        carrier_total = df[carrier_col] + df[non_carrier_col]
        df[f"{slug}_carrier_share"] = safe_divide(
            df[carrier_col], carrier_total, default=np.nan
        )
        df[f"{slug}_carrier_resistance_rate"] = safe_divide(
            df[res_carrier_col], df[carrier_col], default=np.nan
        )
        df[f"{slug}_non_carrier_resistance_rate"] = safe_divide(
            df[res_non_carrier_col], df[non_carrier_col], default=np.nan
        )
    
    # Calculate death cause proportions (if available)
    death_cause_cols = ['deaths_background', 'deaths_sepsis', 'deaths_drug_toxicity']
    if all(col in df.columns for col in death_cause_cols):
        df['prop_deaths_background'] = safe_divide(df['deaths_background'], df['total_deaths'])
        df['prop_deaths_sepsis'] = safe_divide(df['deaths_sepsis'], df['total_deaths'])
        df['prop_deaths_drug_toxicity'] = safe_divide(df['deaths_drug_toxicity'], df['total_deaths'])
    
    logger.info(f"Data preprocessing complete. Shape: {df.shape}")
    return df

def extract_bacteria_list_from_csv(df: pd.DataFrame) -> list:
    """
    Dynamically extract the list of bacteria from CSV column headers.
    This replaces hardcoded bacteria lists and automatically adapts to any BACTERIA_LIST configuration.
    """
    bacteria_list = []
    for col in df.columns:
        if col.endswith('_currently_infected'):
            bacteria_name = col.replace('_currently_infected', '')
            bacteria_list.append(bacteria_name)
    
    bacteria_list.sort()  # For consistent ordering
    return bacteria_list

def extract_drug_list_from_csv(df: pd.DataFrame) -> list:
    """
    Dynamically extract the list of drugs from CSV column headers.
    """
    drugs = []
    for col in df.columns:
        if col.endswith('_currently_on_drug'):
            drug_name = col.replace('_currently_on_drug', '')
            drugs.append(drug_name)
    
    drugs.sort()  # For consistent ordering
    return drugs

def extract_resistance_mechanisms_from_csv(df: pd.DataFrame) -> list:
    """
    Dynamically extract resistance mechanisms from CSV column headers.
    """
    mechanisms = set()
    for col in df.columns:
        if '_infected_with_' in col:
            # Extract mechanism name from column like "escherichia_coli_infected_with_esbl"
            parts = col.split('_infected_with_')
            if len(parts) == 2:
                mechanism = parts[1]
                mechanisms.add(mechanism)

    mechanism_list = sorted(list(mechanisms))
    return mechanism_list