#!/usr/bin/env python3
"""
Polars-based Data Loading for AMR Simulation Output Analysis

This module provides high-performance CSV loading and preprocessing using Polars.
It delivers 2-5x speedup over pandas for large simulation files, while maintaining
backward compatibility by converting to pandas DataFrames for downstream plotting.

Key optimizations:
- Lazy evaluation for efficient query planning
- Parallel CSV parsing
- Columnar processing with SIMD
- Memory-efficient data types
"""

import logging
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np

logger = logging.getLogger(__name__)

# Try to import polars, fall back gracefully if not available
try:
    import polars as pl
    POLARS_AVAILABLE = True
    logger.info("Polars %s loaded successfully", pl.__version__)
except ImportError:
    POLARS_AVAILABLE = False
    pl = None
    logger.warning("Polars not available; falling back to pandas loader")


def load_csv_with_polars(csv_path: Path) -> Optional["pl.DataFrame"]:
    """
    Load a CSV file using Polars for maximum performance.
    
    Args:
        csv_path: Path to the CSV file
        
    Returns:
        Polars DataFrame or None if loading failed
    """
    if not POLARS_AVAILABLE:
        return None
        
    if not csv_path.exists():
        logger.error(f"CSV file not found: {csv_path}")
        return None
    
    try:
        # Use lazy loading with optimized settings
        # Note: Polars automatically uses all available threads
        df = pl.scan_csv(
            csv_path,
            infer_schema_length=50000,  # Sample more rows for better type inference
            low_memory=True,  # Prefer memory efficiency over speed
            rechunk=False,  # Skip rechunking to reduce memory usage
            ignore_errors=True,  # Be lenient with parsing errors
        ).collect(streaming=True)  # Use streaming mode for large files
        
        logger.info(f"Loaded {len(df)} rows from {csv_path} using Polars")
        print(f"Loaded {len(df)} time steps of simulation data (Polars)")
        return df
        
    except Exception as e:
        logger.error(f"Polars CSV load failed: {e}")
        return None


def polars_safe_divide(
    df: "pl.DataFrame",
    numerator_col: str,
    denominator_col: str,
    result_col: str,
    default: float = float('nan'),
) -> "pl.DataFrame":
    """
    Safe division in Polars, returning default value when denominator is zero.
    """
    return df.with_columns(
        pl.when(pl.col(denominator_col) != 0)
        .then(pl.col(numerator_col) / pl.col(denominator_col))
        .otherwise(pl.lit(default))
        .alias(result_col)
    )


def preprocess_with_polars(df: "pl.DataFrame", enable_microbiome_aggregates: bool = True) -> "pl.DataFrame":
    """
    Preprocess simulation data using Polars for maximum performance.
    
    This mirrors the preprocessing logic in data_loader.py but uses Polars
    for 2-5x speedup on large datasets.
    
    Args:
        df: Polars DataFrame with raw simulation data
        enable_microbiome_aggregates: Whether to compute microbiome acquisition/clearance totals
        
    Returns:
        Polars DataFrame with additional calculated columns
    """
    if not POLARS_AVAILABLE or df is None:
        return df
    
    logger.info("Starting Polars preprocessing")
    
    # Convert to lazy frame for query optimization
    lf = df.lazy()
    
    # Age group proportions
    age_cols = ['num_age_0_5', 'num_age_6_14', 'num_age_15_49', 'num_age_50_79', 'num_age_80plus']
    if all(col in df.columns for col in age_cols) and 'total_population' in df.columns:
        lf = lf.with_columns([
            (pl.col('num_age_0_5') / pl.col('total_population')).alias('prop_age_0_5'),
            (pl.col('num_age_6_14') / pl.col('total_population')).alias('prop_age_6_14'),
            (pl.col('num_age_15_49') / pl.col('total_population')).alias('prop_age_15_49'),
            (pl.col('num_age_50_79') / pl.col('total_population')).alias('prop_age_50_79'),
            (pl.col('num_age_80plus') / pl.col('total_population')).alias('prop_age_80plus'),
        ])
    
    # Infected and on drug proportion
    if 'currently_infected_and_on_drug_count' in df.columns and 'total_currently_infected' in df.columns:
        lf = lf.with_columns(
            pl.when(pl.col('total_currently_infected') != 0)
            .then(pl.col('currently_infected_and_on_drug_count') / pl.col('total_currently_infected'))
            .otherwise(pl.lit(float('nan')))
            .alias('infected_and_on_drug_proportion')
        )
    
    # Acquisition-person count over the rolling 365-day window.
    if 'infection_acquisition_people_past_year' in df.columns and 'total_population' in df.columns:
        lf = lf.with_columns(
            (pl.col('infection_acquisition_people_past_year') / pl.col('total_population'))
            .alias('infection_acquisition_people_past_year_proportion')
        )
    
    # Death proportions (past year)
    death_year_cols = [
        ('deaths_past_year', 'deaths_past_year_proportion'),
        ('deaths_background_past_year', 'deaths_background_past_year_proportion'),
        ('deaths_sepsis_past_year', 'deaths_sepsis_past_year_proportion'),
        ('deaths_infection_non_sepsis_past_year', 'deaths_infection_non_sepsis_past_year_proportion'),
        ('deaths_drug_toxicity_past_year', 'deaths_drug_toxicity_past_year_proportion'),
    ]
    
    for death_col, prop_col in death_year_cols:
        if death_col in df.columns and 'total_population' in df.columns:
            lf = lf.with_columns(
                (pl.col(death_col) / pl.col('total_population')).alias(prop_col)
            )
    
    # Time in years
    lf = lf.with_columns(
        (pl.col('time_step') / 365.0).alias('time_in_years')
    )
    
    # Basic proportions
    lf = lf.with_columns([
        pl.when(pl.col('total_population') != 0)
        .then(pl.col('total_currently_infected') / pl.col('total_population'))
        .otherwise(pl.lit(float('nan')))
        .alias('infection_proportion'),
        
        pl.when(pl.col('total_population') != 0)
        .then(pl.col('total_deaths') / pl.col('total_population'))
        .otherwise(pl.lit(float('nan')))
        .alias('death_proportion'),
    ])
    
    # Resistance among infected (excluding MDR-TB)
    tb_slug = "mdr_mycobacterium_tuberculosis"
    tb_infected_col = f"{tb_slug}_currently_infected"
    tb_res_carrier_col = f"{tb_slug}_resistant_infected_carrier_count"
    tb_res_non_carrier_col = f"{tb_slug}_resistant_infected_non_carrier_count"
    
    if all(col in df.columns for col in [tb_infected_col, tb_res_carrier_col, tb_res_non_carrier_col]):
        lf = lf.with_columns([
            (pl.col(tb_infected_col).fill_null(0)).alias('_tb_infected'),
            (pl.col(tb_res_carrier_col).fill_null(0) + pl.col(tb_res_non_carrier_col).fill_null(0)).alias('_tb_resistant'),
        ]).with_columns([
            (pl.col('total_currently_infected') - pl.col('_tb_infected')).alias('_infected_excl_tb'),
            (pl.col('total_with_resistance') - pl.col('_tb_resistant')).alias('_resistance_excl_tb'),
        ]).with_columns(
            pl.when(pl.col('_infected_excl_tb') != 0)
            .then(pl.col('_resistance_excl_tb') / pl.col('_infected_excl_tb'))
            .otherwise(pl.lit(float('nan')))
            .alias('resistance_among_infected')
        ).drop(['_tb_infected', '_tb_resistant', '_infected_excl_tb', '_resistance_excl_tb'])
    else:
        lf = lf.with_columns(
            pl.when(pl.col('total_currently_infected') != 0)
            .then(pl.col('total_with_resistance') / pl.col('total_currently_infected'))
            .otherwise(pl.lit(float('nan')))
            .alias('resistance_among_infected')
        )
    
    # Infection duration proportions
    if 'infected_10_days_count' in df.columns and 'infected_21_days_count' in df.columns and 'total_currently_infected' in df.columns:
        lf = lf.with_columns([
            pl.when(pl.col('total_currently_infected') != 0)
            .then(pl.col('infected_10_days_count') / pl.col('total_currently_infected'))
            .otherwise(pl.lit(float('nan')))
            .alias('infected_10_days_proportion'),
            
            pl.when(pl.col('total_currently_infected') != 0)
            .then(pl.col('infected_21_days_count') / pl.col('total_currently_infected'))
            .otherwise(pl.lit(float('nan')))
            .alias('infected_21_days_proportion'),
        ])
    
    # Sepsis proportion
    if 'number_with_sepsis' in df.columns and 'total_currently_infected' in df.columns:
        lf = lf.with_columns(
            pl.when(pl.col('total_currently_infected') != 0)
            .then(pl.col('number_with_sepsis') / pl.col('total_currently_infected'))
            .otherwise(pl.lit(float('nan')))
            .alias('sepsis_among_infected_proportion')
        )
    
    # Collect intermediate result for column-based operations
    df = lf.collect()
    
    # Process carrier/non-carrier metrics
    carrier_suffix = '_infected_carrier_count'
    carrier_cols = [col for col in df.columns if col.endswith(carrier_suffix)]
    
    new_columns = []
    for carrier_col in carrier_cols:
        slug = carrier_col[:-len(carrier_suffix)]
        non_carrier_col = f"{slug}_infected_non_carrier_count"
        res_carrier_col = f"{slug}_resistant_infected_carrier_count"
        res_non_carrier_col = f"{slug}_resistant_infected_non_carrier_count"
        
        if not all(col in df.columns for col in [non_carrier_col, res_carrier_col, res_non_carrier_col]):
            continue
        
        carrier_total = df[carrier_col] + df[non_carrier_col]
        
        new_columns.extend([
            (pl.when(carrier_total != 0)
             .then(pl.col(carrier_col) / carrier_total)
             .otherwise(pl.lit(float('nan')))
             .alias(f"{slug}_carrier_share")),
            
            (pl.when(pl.col(carrier_col) != 0)
             .then(pl.col(res_carrier_col) / pl.col(carrier_col))
             .otherwise(pl.lit(float('nan')))
             .alias(f"{slug}_carrier_resistance_rate")),
            
            (pl.when(pl.col(non_carrier_col) != 0)
             .then(pl.col(res_non_carrier_col) / pl.col(non_carrier_col))
             .otherwise(pl.lit(float('nan')))
             .alias(f"{slug}_non_carrier_resistance_rate")),
        ])
    
    if new_columns:
        df = df.with_columns(new_columns)
    
    # Resistant microbiome shares
    resistant_microbiome_suffix = '_presence_microbiome_resistant'
    base_microbiome_suffix = '_presence_microbiome'
    resistant_cols = [col for col in df.columns if col.endswith(resistant_microbiome_suffix)]
    
    new_columns = []
    for resistant_col in resistant_cols:
        slug = resistant_col[:-len(resistant_microbiome_suffix)]
        base_col = f"{slug}{base_microbiome_suffix}"
        
        if base_col not in df.columns:
            continue
        
        new_columns.append(
            pl.when(pl.col(base_col) != 0)
            .then(pl.col(resistant_col) / pl.col(base_col))
            .otherwise(pl.lit(float('nan')))
            .alias(f"{slug}_resistant_microbiome_share")
        )
        
        # Resistant infection share
        infected_carrier_col = f"{slug}_infected_carrier_count"
        infected_non_carrier_col = f"{slug}_infected_non_carrier_count"
        res_infected_carrier_col = f"{slug}_resistant_infected_carrier_count"
        res_infected_non_carrier_col = f"{slug}_resistant_infected_non_carrier_count"
        
        if all(col in df.columns for col in [infected_carrier_col, infected_non_carrier_col,
                                              res_infected_carrier_col, res_infected_non_carrier_col]):
            infected_total_expr = pl.col(infected_carrier_col) + pl.col(infected_non_carrier_col)
            res_infected_total_expr = pl.col(res_infected_carrier_col) + pl.col(res_infected_non_carrier_col)
            
            new_columns.append(
                pl.when(infected_total_expr != 0)
                .then(res_infected_total_expr / infected_total_expr)
                .otherwise(pl.lit(float('nan')))
                .alias(f"{slug}_resistant_infection_share")
            )
    
    if new_columns:
        df = df.with_columns(new_columns)
    
    # Carriage duration distributions
    duration_labels = ["0_29", "30_89", "90_179", "180_359", "360_plus"]
    duration_prefix = "_carriage_duration_days_"
    base_suffix = f"{duration_prefix}{duration_labels[0]}"
    duration_base_cols = [col for col in df.columns if col.endswith(base_suffix)]
    
    for base_col in duration_base_cols:
        slug = base_col[:-len(base_suffix)]
        duration_cols = [f"{slug}{duration_prefix}{label}" for label in duration_labels]
        
        if not all(col in df.columns for col in duration_cols):
            continue
        
        total_col = f"{slug}_carriage_duration_total"
        total_expr = sum(pl.col(col) for col in duration_cols)
        
        new_columns = [total_expr.alias(total_col)]
        for label, col_name in zip(duration_labels, duration_cols):
            share_col = f"{slug}_carriage_duration_share_{label}"
            new_columns.append(
                pl.when(total_expr != 0)
                .then(pl.col(col_name) / total_expr)
                .otherwise(pl.lit(float('nan')))
                .alias(share_col)
            )
        
        df = df.with_columns(new_columns)
    
    # Death cause proportions
    death_cause_cols = [
        ('deaths_background', 'prop_deaths_background'),
        ('deaths_sepsis', 'prop_deaths_sepsis'),
        ('deaths_infection_non_sepsis', 'prop_deaths_infection_non_sepsis'),
        ('deaths_drug_toxicity', 'prop_deaths_drug_toxicity'),
    ]
    
    if all(col in df.columns for col, _ in death_cause_cols):
        denominator_col = 'total_deaths' if 'total_deaths' in df.columns else None
        if denominator_col:
            new_columns = [
                pl.when(pl.col(denominator_col) != 0)
                .then(pl.col(col) / pl.col(denominator_col))
                .otherwise(pl.lit(float('nan')))
                .alias(prop)
                for col, prop in death_cause_cols
            ]
            df = df.with_columns(new_columns)
    
    # Microbiome acquisition metrics (if enabled)
    if enable_microbiome_aggregates:
        on_suffix = '_microbiome_acquisitions_on_drug'
        off_suffix = '_microbiome_acquisitions_off_drug'
        on_columns = [col for col in df.columns if col.endswith(on_suffix)]
        
        for on_col in on_columns:
            slug = on_col[:-len(on_suffix)]
            off_col = f"{slug}{off_suffix}"
            
            if off_col not in df.columns:
                continue
            
            total_col = f"{slug}_microbiome_acquisitions_total"
            share_on_col = f"{slug}_microbiome_acquisitions_share_on_drug"
            share_off_col = f"{slug}_microbiome_acquisitions_share_off_drug"
            
            total_expr = pl.col(on_col) + pl.col(off_col)
            
            new_columns = [
                total_expr.alias(total_col),
                pl.when(total_expr != 0)
                .then(pl.col(on_col) / total_expr)
                .otherwise(pl.lit(float('nan')))
                .alias(share_on_col),
                pl.when(total_expr != 0)
                .then(pl.col(off_col) / total_expr)
                .otherwise(pl.lit(float('nan')))
                .alias(share_off_col),
            ]
            
            if 'total_population' in df.columns:
                new_columns.extend([
                    (pl.col(on_col) / pl.col('total_population') * 1e5).alias(f"{slug}_microbiome_acquisitions_on_drug_per_100k"),
                    (pl.col(off_col) / pl.col('total_population') * 1e5).alias(f"{slug}_microbiome_acquisitions_off_drug_per_100k"),
                    (total_expr / pl.col('total_population') * 1e5).alias(f"{slug}_microbiome_acquisitions_total_per_100k"),
                ])
            
            df = df.with_columns(new_columns)
        
        # Microbiome clearance metrics
        clearance_on_suffix = '_microbiome_clearances_on_drug'
        clearance_off_suffix = '_microbiome_clearances_off_drug'
        clearance_on_columns = [col for col in df.columns if col.endswith(clearance_on_suffix)]
        
        for on_col in clearance_on_columns:
            slug = on_col[:-len(clearance_on_suffix)]
            off_col = f"{slug}{clearance_off_suffix}"
            
            if off_col not in df.columns:
                continue
            
            total_col = f"{slug}_microbiome_clearances_total"
            share_on_col = f"{slug}_microbiome_clearances_share_on_drug"
            share_off_col = f"{slug}_microbiome_clearances_share_off_drug"
            
            total_expr = pl.col(on_col) + pl.col(off_col)
            
            new_columns = [
                total_expr.alias(total_col),
                pl.when(total_expr != 0)
                .then(pl.col(on_col) / total_expr)
                .otherwise(pl.lit(float('nan')))
                .alias(share_on_col),
                pl.when(total_expr != 0)
                .then(pl.col(off_col) / total_expr)
                .otherwise(pl.lit(float('nan')))
                .alias(share_off_col),
            ]
            
            if 'total_population' in df.columns:
                new_columns.extend([
                    (pl.col(on_col) / pl.col('total_population') * 1e5).alias(f"{slug}_microbiome_clearances_on_drug_per_100k"),
                    (pl.col(off_col) / pl.col('total_population') * 1e5).alias(f"{slug}_microbiome_clearances_off_drug_per_100k"),
                    (total_expr / pl.col('total_population') * 1e5).alias(f"{slug}_microbiome_clearances_total_per_100k"),
                ])
            
            df = df.with_columns(new_columns)

    # Match the pandas carrier/non-carrier acquisition derivations.
    carrier_suffix = '_infection_acquisition_events_carrier_at_acquisition'
    non_carrier_suffix = '_infection_acquisition_events_non_carrier_at_acquisition'
    carrier_columns = [col for col in df.columns if col.endswith(carrier_suffix)]

    for carrier_col in carrier_columns:
        slug = carrier_col[:-len(carrier_suffix)]
        non_carrier_col = f"{slug}{non_carrier_suffix}"
        presence_col = f"{slug}_presence_microbiome"
        if non_carrier_col not in df.columns or presence_col not in df.columns:
            continue

        carrier_rolling = pl.col(carrier_col).rolling_sum(window_size=365, min_samples=1)
        non_carrier_rolling = pl.col(non_carrier_col).rolling_sum(
            window_size=365,
            min_samples=1,
        )
        total_rolling = carrier_rolling + non_carrier_rolling
        derived = [
            carrier_rolling.alias(
                f"{slug}_infection_acquisition_events_carrier_rolling_year"
            ),
            non_carrier_rolling.alias(
                f"{slug}_infection_acquisition_events_non_carrier_rolling_year"
            ),
            pl.when(total_rolling > 0)
            .then(carrier_rolling / total_rolling)
            .otherwise(pl.lit(float('nan')))
            .alias(f"{slug}_new_infection_share_from_carriers"),
        ]

        if 'total_population' in df.columns:
            non_carrier_population = (
                pl.col('total_population').cast(pl.Float64)
                - pl.col(presence_col).cast(pl.Float64)
            ).clip(lower_bound=0.0)
            derived.extend([
                pl.when(pl.col(presence_col) > 0)
                .then(carrier_rolling / pl.col(presence_col) * 1e5)
                .otherwise(pl.lit(float('nan')))
                .alias(f"{slug}_infection_acquisition_events_per_100k_carriers"),
                pl.when(non_carrier_population > 0)
                .then(non_carrier_rolling / non_carrier_population * 1e5)
                .otherwise(pl.lit(float('nan')))
                .alias(f"{slug}_infection_acquisition_events_per_100k_non_carriers"),
            ])

        df = df.with_columns(derived)
    
    logger.info("Polars preprocessing completed")
    return df


def polars_to_pandas(df: "pl.DataFrame"):
    """
    Convert Polars DataFrame to pandas DataFrame for downstream compatibility.
    
    Uses optimized conversion with memory-efficient data types.
    
    Args:
        df: Polars DataFrame
        
    Returns:
        pandas DataFrame
    """
    if df is None:
        return None
    
    try:
        import gc
        
        # Use use_pyarrow_extension_array=True for memory-efficient conversion
        # This avoids creating large intermediate numpy arrays
        result = df.to_pandas(use_pyarrow_extension_array=True)
        
        # Force garbage collection after conversion
        gc.collect()
        
        return result
    except Exception as e:
        logger.warning(f"PyArrow extension array conversion failed: {e}, trying standard conversion")
        try:
            # Fallback to standard conversion
            return df.to_pandas()
        except Exception as e2:
            logger.error(f"Failed to convert Polars to pandas: {e2}")
            return None


def is_polars_available() -> bool:
    """Check if Polars is available for use."""
    return POLARS_AVAILABLE
