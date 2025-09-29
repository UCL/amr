#!/usr/bin/env python3
"""
Empirical data loading functionality for AMR simulation analysis.

This module contains functions for loading empirical calibration data
extracted from the original analyze_simulation.py script.
"""

import pandas as pd
from pathlib import Path
from .normalizers import normalize_name_for_empirical_matching

def load_empirical_calibration_data():
    """
    Load empirical calibration data for overlay on simulation plots.
    Uses real surveillance data from WHO GLASS, ECDC EARS-Net, CDC NARMS, and GBD Study.
    Enhanced with integrated surveillance data from major global sources.
    Returns dictionary with data for drug usage, resistance, incidence, and deaths.
    """
    
    # Try to use enhanced integrated loader if available
    try:
        import sys
        from pathlib import Path
        
        # Add project root to path
        project_root = Path(__file__).parent.parent.parent
        if str(project_root) not in sys.path:
            sys.path.insert(0, str(project_root))
            
        from enhanced_empirical_loader import load_integrated_empirical_data
        
        print("\n🚀 Loading integrated empirical data (WHO GLASS, ECDC EARS-Net, Australian NNDSS, CDDEP ResistanceMap)...")
        integrated_data = load_integrated_empirical_data()
        
        # Show enhanced coverage
        total_records = sum(len(df) if df is not None else 0 for df in integrated_data.values())
        print(f"✅ Enhanced empirical data loaded: {total_records:,} total records from integrated surveillance sources")
        
        return integrated_data
        
    except Exception as e:
        print(f"⚠️  Enhanced loader unavailable ({e}), using standard calibration data...")
    
    # Fallback to standard calibration loading
    empirical_data = {
        'drug_usage': None,
        'resistance': None, 
        'incidence': None,
        'deaths': None,
        # NEW TIER 1 CLINICAL METRICS
        'drug_failure': None,
        'mic_values': None,
        'hospital_incidence': None
    }
    
    # Standard empirical files with real surveillance data
    empirical_files = {
        'drug_usage': 'data/empirical/calibration_drug_usage_empirical.csv',
        'resistance': 'data/empirical/calibration_resistance_empirical.csv',
        'incidence': 'data/empirical/calibration_infection_incidence_empirical.csv', 
        'deaths': 'data/empirical/calibration_deaths_empirical.csv',
        # NEW TIER 1 CLINICAL FILES
        'drug_failure': 'data/empirical/calibration_drug_failure_empirical.csv',
        'mic_values': 'data/empirical/calibration_mic_empirical.csv',
        'hospital_incidence': 'data/empirical/calibration_hospital_incidence_empirical.csv'
    }
    
    # Check if empirical files exist
    empirical_files_exist = all(Path(f).exists() for f in empirical_files.values())
    
    # Check for empirical enhancement module
    try:
        from empirical_enhancement import enhance_empirical_data
        HAS_EMPIRICAL_ENHANCEMENT = True
    except ImportError:
        HAS_EMPIRICAL_ENHANCEMENT = False
    
    FORCE_REGENERATE_EMPIRICAL = False  # Set based on config if needed
    
    # Generate empirical data if missing or forced regeneration
    if HAS_EMPIRICAL_ENHANCEMENT and (not empirical_files_exist or FORCE_REGENERATE_EMPIRICAL):
        print("\n🚀 Generating empirical data with real surveillance patterns...")
        try:
            enhance_empirical_data(force_regenerate=FORCE_REGENERATE_EMPIRICAL)
            print("Empirical data ready with WHO GLASS, ECDC, CDC, and GBD patterns")
            empirical_files_exist = True
        except Exception as e:
            print(f"WARNING: Empirical data generation failed: {e}")
            print("[FOLDER] Analysis will proceed without empirical overlays...")
            return empirical_data
    
    print("\nLoading empirical calibration data (real surveillance patterns)...")
    
    for data_type, filename in empirical_files.items():
        try:
            if Path(filename).exists():
                df = pd.read_csv(filename)
                
                # Fix: Drug failure data has bacteria/drug columns swapped
                if data_type == 'drug_failure' and 'bacteria' in df.columns and 'drug' in df.columns:
                    # Check if swap is needed (if bacteria column contains drug names)
                    sample_bacteria = df['bacteria'].iloc[0] if len(df) > 0 else ""
                    if any(drug_name in sample_bacteria for drug_name in ['amoxicillin', 'ciprofloxacin', 'vancomycin']):
                        print(f"   Fixing swapped bacteria/drug columns in {filename}")
                        df['bacteria'], df['drug'] = df['drug'], df['bacteria']
                
                empirical_data[data_type] = df
                
                # Show empirical data coverage
                empirical_indicator = ""
                if 'notes' in df.columns:
                    empirical_count = len([r for _, r in df.iterrows() 
                                         if any(pattern in str(r.get('notes', '')) 
                                               for pattern in ['who_glass', 'ecdc', 'cdc', 'gbd', 'integrated'])])
                    if empirical_count > 0:
                        empirical_indicator = f" ({empirical_count:,} real surveillance records, {empirical_count/len(df)*100:.1f}%)"
                    else:
                        empirical_indicator = " (baseline synthetic data)"
                
                print(f"   Loaded {len(df):,} records from {filename}{empirical_indicator}")
            else:
                print(f"   WARNING: {filename} not found, skipping empirical overlay for {data_type}")
        except Exception as e:
            print(f"   ERROR loading {filename}: {e}")
    
    return empirical_data

def get_empirical_data_for_plot(empirical_df, drug=None, bacteria=None, region=None, metric_type=None, data_source=None):
    """
    Extract empirical data points for a specific drug/bacteria/region combination.
    If region=None, averages across all regions for global plots.
    Returns (years, means, p5, p95) for plotting.
    
    Args:
        empirical_df: The empirical data DataFrame
        drug: Drug name to filter by
        bacteria: Bacteria name to filter by  
        region: Region name to filter by
        metric_type: Metric type to filter by
        data_source: Source of empirical data ('drug_failure', 'mic_values', etc.) for context-aware name normalization
    """
    if empirical_df is None:
        return None, None, None, None
    
    # Filter data based on parameters
    filtered_df = empirical_df.copy()
    
    if drug is not None:
        # Normalize drug name for matching
        normalized_drug = normalize_name_for_empirical_matching(drug, entity_type='drug', data_source=data_source)
        filtered_df = filtered_df[filtered_df['drug'] == normalized_drug]
    if bacteria is not None:
        # Normalize bacteria name for matching (context-aware)
        normalized_bacteria = normalize_name_for_empirical_matching(bacteria, entity_type='bacteria', data_source=data_source)
        filtered_df = filtered_df[filtered_df['bacteria'] == normalized_bacteria]
    if metric_type is not None:
        if 'metric' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['metric'] == metric_type]
    
    if len(filtered_df) == 0:
        return None, None, None, None
    
    # Handle regional filtering or averaging
    if region is not None:
        # Filter for specific region
        filtered_df = filtered_df[filtered_df['region'] == region]
        if len(filtered_df) == 0:
            return None, None, None, None
    else:
        # Average across all regions for global plots
        # Determine which columns to aggregate based on data type
        agg_dict = {}
        if 'mean' in filtered_df.columns:
            agg_dict['mean'] = 'mean'
        elif 'failure_rate' in filtered_df.columns:
            agg_dict['failure_rate'] = 'mean'
        elif metric_type and metric_type in filtered_df.columns:
            agg_dict[metric_type] = 'mean'
        
        # Add confidence interval columns if available
        if 'p5' in filtered_df.columns:
            agg_dict['p5'] = 'mean'
        if 'p95' in filtered_df.columns:
            agg_dict['p95'] = 'mean'
        
        if agg_dict:
            grouped = filtered_df.groupby('year').agg(agg_dict).reset_index()
            filtered_df = grouped
        else:
            return None, None, None, None
    
    # Sort by year and extract plotting data
    filtered_df = filtered_df.sort_values('year')
    
    # Convert absolute years to simulation years (simulation starts at 1930)
    sim_years = filtered_df['year'] - 1930
    
    # Handle different column names for different data types
    if 'mean' in filtered_df.columns:
        means = filtered_df['mean'].values
    elif 'failure_rate' in filtered_df.columns:
        means = filtered_df['failure_rate'].values
    elif metric_type and metric_type in filtered_df.columns:
        means = filtered_df[metric_type].values
    else:
        # Fallback: try to find a main value column
        value_cols = ['mean', 'failure_rate', 'mic50', 'incidence_per_1000_days']
        means = None
        for col in value_cols:
            if col in filtered_df.columns:
                means = filtered_df[col].values
                break
        if means is None:
            return None, None, None, None
    
    p5 = filtered_df['p5'].values if 'p5' in filtered_df.columns else None
    p95 = filtered_df['p95'].values if 'p95' in filtered_df.columns else None
    
    return sim_years, means, p5, p95