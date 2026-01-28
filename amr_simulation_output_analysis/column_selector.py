#!/usr/bin/env python3
"""
Column selector for memory-efficient CSV loading.

This module determines which columns are needed based on enabled plot types,
allowing the data loader to skip unnecessary columns and reduce memory usage
by 70-90% for large simulation files.
"""

from typing import List, Set, Optional
import re


# Core columns always needed for any analysis
CORE_COLUMNS = {
    'time_step',
    'policy_option', 
    'run_id',
    'time_in_years',
    'total_population',
    'number_in_hospital',
    'number_severely_immunosuppressed',
    'number_with_sepsis',
    'total_currently_infected',
    'infected_10_days_count',
    'infected_30_days_count',
    'total_with_resistance',
    'currently_taking_drug_count',
    'currently_infected_and_on_drug_count',
    'taking_two_drugs_count',
    'newly_infected_count',
    'newly_infected_with_resistance_count',
    'new_drug_initiations_count',
    'new_drug_initiations_count_infected',
    'newly_infected_past_year',
    'total_deaths',
    'deaths_background',
    'deaths_sepsis',
    'deaths_infection_non_sepsis',
    'deaths_drug_toxicity',
    'deaths_past_year',
    'deaths_background_past_year',
    'deaths_sepsis_past_year',
    'deaths_infection_non_sepsis_past_year',
    'deaths_drug_toxicity_past_year',
    'num_age_0_5',
    'num_age_6_14',
    'num_age_15_49',
    'num_age_50_79',
    'num_age_80plus',
    'num_with_any_bacteria_microbiome',
    'people_on_1_drug',
    'people_on_2_drugs',
    'people_on_3plus_drugs',
    'infected_on_drug_with_previous_failure',
}


# Column patterns needed for grouped plots (not the full resistance matrix)
GROUPED_PLOT_PATTERNS = [
    # Per-bacteria infection counts
    r'.*_currently_infected$',
    r'.*_new_sepsis_cases$',
    r'.*_number_with_sepsis$',
    r'.*_newly_infected_carrier$',
    r'.*_newly_infected_non_carrier$',
    r'.*_newly_infected$',
    r'.*_deaths$',
    r'.*_deaths_past_year$',
    
    # Per-bacteria testing
    r'.*_test_identified$',
    r'.*_test_for_resistance$',
    
    # Per-bacteria infection resolution (all resolution types)
    r'.*_infection_resolution_immune_clearance$',
    r'.*_infection_resolution_drug_assisted_clearance$',
    r'.*_infection_resolution_death_from_sepsis$',
    r'.*_infection_resolution_death_from_background$',
    r'.*_infection_resolution_death_from_infection_non_sepsis$',
    r'.*_infection_resolution_death_from_toxicity$',
    
    # Drug usage columns (aggregated per drug, not per bacteria-drug)
    r'^taking_drug_.*',
    r'^new_initiations_drug_.*',
    
    # Activity R sums (aggregated, not per drug)
    r'.*_activity_r_sum$',
    r'.*_infected_and_on_any_drug$',
    
    # Microbiome (aggregated)
    r'.*_microbiome_presence_count$',
    r'.*_microbiome_acquired_today$',
    r'.*_microbiome_cleared_today$',
    
    # Region columns
    r'^region_.*_population$',
    r'^region_.*_deaths$',
    r'^region_.*_infected$',
    
    # Syndrome columns  
    r'^syndrome_\d+_count$',
    r'^syndrome_\d+_deaths$',
]


# Additional patterns for calibration summary
CALIBRATION_PATTERNS = [
    # Drug class usage
    r'^taking_drug_.*',
    
    # Drug presence columns (needed to detect which drugs are modelled)
    r'.*_currently_on_drug$',                # {drug}_currently_on_drug
    
    # Resistance by bacteria (aggregated any_r, not full matrix)
    r'.*_any_r_gt_0_count$',
    r'.*_majority_r_gt_0_count$',
    
    # Deaths by bacteria
    r'.*_deaths$',
    r'.*_deaths_past_year$',
    
    # === RESISTANCE COLUMNS NEEDED FOR CALIBRATION SUMMARY ===
    # These are bacteria×drug columns but required for resistance benchmarks
    r'.*_sum_any_r_.*',                      # Sum of any_r for each bacteria-drug pair
    r'.*_infected_with_any_r_positive_.*',   # Count of infected with resistance
    r'.*_microbiome_r_positive_.*',          # Microbiome resistance counts
    r'.*_presence_microbiome$',              # Microbiome presence per bacteria
]


# Patterns that are EXCLUDED (the big memory hogs)
# NOTE: Calibration columns (sum_any_r, infected_with_any_r_positive) are now
# explicitly included via CALIBRATION_PATTERNS, so they take precedence.
EXCLUDED_PATTERNS = [
    # Full bacteria-drug resistance matrix - MEAN values (not needed for calibration)
    r'.*_any_r_.*_mean$',           # Skip detailed any_r mean per drug
    r'.*_majority_r_.*_mean$',      # Skip majority_r mean per drug
    r'.*_test_r_.*',                # Skip test_r details
    r'.*_microbiome_r_.*_mean$',    # Skip microbiome_r mean per drug
    
    # MIC distributions (huge - 35 bacteria × 52 drugs × 2 = 3640)
    r'.*_mic_lt_2_.*',
    r'.*_mic_.*_count$',
    
    # Drug scores (35 bacteria × 52 drugs = 1820)
    r'.*_drug_score_.*',
    
    # Detailed resistance source tracking
    r'.*_resistance_source_.*',
    r'.*_new_resistance_.*',
    
    # Per-drug hospital columns (redundant with aggregated)
    r'.*_hospital_.*_drug_.*',
]


def get_required_columns(
    all_columns: List[str],
    include_grouped_plots: bool = True,
    include_calibration: bool = True,
    include_detail_plots: bool = False,
) -> List[str]:
    """
    Determine which columns to load based on enabled analysis types.
    
    Args:
        all_columns: List of all available column names in the CSV
        include_grouped_plots: Include columns for grouped figures 1-10
        include_calibration: Include columns for calibration_summary.txt
        include_detail_plots: Include columns for detailed individual plots
        
    Returns:
        List of column names to load
    """
    required: Set[str] = set(CORE_COLUMNS)
    
    # Compile patterns
    include_patterns = []
    if include_grouped_plots:
        include_patterns.extend(GROUPED_PLOT_PATTERNS)
    if include_calibration:
        include_patterns.extend(CALIBRATION_PATTERNS)
    
    # If detail plots enabled, we need everything (no memory savings)
    if include_detail_plots:
        return list(all_columns)
    
    # Compile exclude patterns
    exclude_regexes = [re.compile(p) for p in EXCLUDED_PATTERNS]
    include_regexes = [re.compile(p) for p in include_patterns]
    
    for col in all_columns:
        # Skip if in core (already added)
        if col in required:
            continue
        
        # Check if matches include pattern FIRST (include takes precedence)
        included = any(rx.match(col) for rx in include_regexes)
        if included:
            required.add(col)
            continue
            
        # Only check exclude if not explicitly included
        excluded = any(rx.match(col) for rx in exclude_regexes)
        if excluded:
            continue
            
        # For columns not matching any pattern, include them by default
        # (they're likely small core columns we missed)
        if included:
            required.add(col)
    
    # Ensure we have columns in the order they appear in CSV
    result = [col for col in all_columns if col in required]
    
    return result


def estimate_memory_savings(total_columns: int, selected_columns: int) -> str:
    """Return a human-readable memory savings estimate."""
    if total_columns == 0:
        return "N/A"
    pct_kept = (selected_columns / total_columns) * 100
    pct_saved = 100 - pct_kept
    return f"Loading {selected_columns:,} of {total_columns:,} columns ({pct_kept:.1f}% kept, ~{pct_saved:.0f}% memory saved)"
