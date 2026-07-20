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
    'new_sepsis_cases',
    'total_currently_infected',
    'infected_10_days_count',
    'infected_21_days_count',
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
    
    # Day-7 drug initiation columns (for grouped figure 7)
    r'.*_day_7_evaluations$',
    r'.*_day_7_drug_used$',
    
    # Drug usage columns (aggregated per drug, not per bacteria-drug)
    r'^taking_drug_.*',
    r'^new_initiations_drug_.*',
    
    # Activity R sums (aggregated, not per drug)
    r'.*_activity_r_sum$',
    r'.*_max_possible_activity_r_sum$',
    r'.*_infected_and_on_any_drug$',
    
    # Microbiome (aggregated)
    r'.*_microbiome_presence_count$',
    r'.*_microbiome_acquired_today$',
    r'.*_microbiome_cleared_today$',
    
    # Region columns (for grouped figure 8)
    r'^region_.*_population$',
    r'^region_.*_deaths$',
    r'^region_.*_infected$',
    r'^north_america_population$',
    r'^south_america_population$',
    r'^africa_population$',
    r'^asia_population$',
    r'^europe_population$',
    r'^oceania_population$',
    
    # Drug failure events (for grouped figure 8)
    r'.*_drug_failure_events_.*',
    
    # Syndrome columns (for grouped figure 8)
    r'^syndrome_\d+_infected$',
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
    r'.*_new_resistance_.*',                 # Breakdown of new resistance sources
    r'.*_asymptomatic_microbiome_hgt_events', # Asymptomatic HGT events
    
    # Additional infection locus columns needed for calibration mapping
    r'.*_newly_infected_(north_america|south_america|africa|asia|europe|oceania)$',
    r'.*_newly_infected_hospital_.*',
    r'.*_newly_infected_any_r_hospital$',
    r'.*_newly_infected_any_r_community$',
    r'.*_newly_infected_carrier$',
    r'.*_newly_infected_non_carrier$',
    r'.*_newly_infected_under_5$',
    r'.*_newly_infected_over_65$',
    r'^new_active_infections_by_bacteria$',
    r'.*_deaths_under_5$',
    r'.*_deaths_over_65$',
    r'.*_deaths_hospital_acquired$',
    r'.*_deaths_community_acquired$',
    r'.*_resistant_infected_hospital_count$',
    r'.*_resistant_infected_community_count$',
    r'.*_presence_microbiome_resistant$',

    # Age-specific infection death rates by region (sepsis + infection_non_sepsis, 60 columns)
    r'^.*_prop_age_(0_5|6_14|15_49|50_79|80plus)_deaths_sepsis$',
    r'^.*_prop_age_(0_5|6_14|15_49|50_79|80plus)_deaths_infection_non_sepsis$',
    # Age proportion denominators for computing per-age-group rates (30 columns)
    r'^.*_prop_age_(0_5|6_14|15_49|50_79|80plus)$',
    # Regional population totals (denominator fallback; also in GROUPED_PLOT_PATTERNS)
    r'^(north_america|south_america|africa|asia|europe|oceania)_population$',
]


# Patterns for specific detail plots (opt-in, not all columns)
# Each detail plot type has its own pattern set with estimated column counts
# LIGHTWEIGHT plots (~234 columns each): Can enable multiple together
# HEAVY plots (1000+ columns): Enable one at a time or disable grouped_plots
DETAIL_PLOT_PATTERNS = {
    # === LIGHTWEIGHT PLOTS (~234 columns each, bacteria×region) ===
    # death_rate_by_bacteria_region - needs deaths_infected per region
    'death_rate_by_bacteria_region': [
        r'.*_deaths_infected_.*',            # {bacteria}_deaths_infected_{region}
    ],
    # population_mortality_by_bacteria_region - same columns
    'population_mortality_by_bacteria_region': [
        r'.*_deaths_infected_.*',
    ],
    # incidence_of_infection - needs regional infection columns
    'incidence_of_infection': [
        r'.*_newly_infected_.*',             # {bacteria}_newly_infected_{region}
    ],
    # drug_failure_rate_by_bacteria_region
    'drug_failure_rate_by_bacteria_region': [
        r'.*_drug_failure_events_.*',        # {bacteria}_drug_failure_events_{region}
    ],
    # death_rate_by_region (aggregated, very small)
    'death_rate_by_region': [
        r'^.*_deaths_.*$',                   # Regional death columns
    ],
    # age_distribution_by_region (uses core columns, no extra needed)
    'age_distribution_by_region': [],
    # death_rate_by_syndrome_region
    'death_rate_by_syndrome_region': [
        r'^syndrome_\d+_deaths_.*',          # Syndrome death columns by region
    ],
    # infection_resolution_by_bacteria (already in grouped patterns)
    'infection_resolution_by_bacteria': [],
    # source_of_new_resistance_by_drug_bacteria
    'source_of_new_resistance_by_drug_bacteria': [
        #r'.*_new_resistance_.*',
    ],
    
    # === MODERATE PLOTS (~500-1000 columns) ===
    # incidence_of_infection_hospital
    'incidence_of_infection_hospital': [
        r'.*_newly_infected_hospital_.*',    # {bacteria}_newly_infected_hospital_{region}
        r'.*_hospital_population$',          # Regional hospital populations
    ],
    # proportion_of_people_infected_with_each_bacteria (uses core columns)
    'proportion_of_people_infected_with_each_bacteria': [],
    # proportion_of_people_taking_each_drug (uses taking_drug columns, already loaded)
    'proportion_of_people_taking_each_drug': [],
    # mean_activity_r_by_bacteria (uses activity_r_sum, already loaded)
    'mean_activity_r_by_bacteria': [],
    # microbiome plots (uses microbiome columns, already loaded for calibration)
    'proportion_of_population_with_microbiome_presence_bacteria': [
        r'.*_presence_microbiome$',
    ],
    'microbiome_acquisition_on_off_drug': [
        r'.*_microbiome_acquisitions_on_drug$',
        r'.*_microbiome_acquisitions_off_drug$',
    ],
    'microbiome_clearance_on_off_drug': [
        r'.*_microbiome_clearances_on_drug$',
        r'.*_microbiome_clearances_off_drug$',
    ],
    'proportion_of_microbiome_presence_with_resistance_by_drug': [
        r'.*_presence_microbiome$',
        r'.*_presence_microbiome_resistant$',
        r'.*_microbiome_r_positive_.*',
    ],
    'microbiome_resistance_microbiome_vs_infection': [
        r'.*_presence_microbiome$',
        r'.*_presence_microbiome_resistant$',
        r'.*_infected_carrier_count$',
        r'.*_infected_non_carrier_count$',
        r'.*_resistant_infected_carrier_count$',
        r'.*_resistant_infected_non_carrier_count$',
    ],
    'carrier_infection_share': [
        r'.*_currently_infected$',
        r'.*_infected_carrier_count$',
        r'.*_infected_non_carrier_count$',
    ],
    'carrier_vs_non_carrier_incidence': [
        r'.*_presence_microbiome$',
        r'.*_newly_infected_carrier$',
        r'.*_newly_infected_non_carrier$',
    ],
    'carriage_duration_distribution': [
        r'.*_carriage_duration_days_.*',
    ],
    
    # === HEAVY PLOTS (1000+ columns, bacteria×drug matrix) ===
    # mean_mic_by_drug_for_each_bacteria (~3640 columns - MIC distributions)
    'mean_mic_by_drug_for_each_bacteria': [
        r'.*_mic_lt_2_.*',
        r'.*_mic_.*_count$',
        r'.*_infected_and_mic_lt2_.*',       # MIC<2 count columns
    ],
    # for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2
    'for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2': [
        r'.*_infected_and_mic_lt2_.*',
    ],
    # drug_score plots (~1820 columns - bacteria×drug scores)
    'drug_score_analysis_by_bacteria': [
        r'.*_drug_score_.*',
    ],
    'drug_score_summary': [
        r'.*_drug_score_.*',
    ],
    # mean_any_r_by_drug plots (~1820 columns)
    'mean_any_r_by_drug_for_each_bacteria': [
        r'.*_any_r_.*_mean$',
    ],
    # resistance_mechanism_by_bacteria
    'resistance_mechanism_by_bacteria': [
        r'.*_resistance_source_.*',
        #r'.*_new_resistance_.*',
    ],
    # source_of_new_resistance_by_drug_bacteria
    'source_of_new_resistance_by_drug_bacteria': [
        r'.*_resistance_source_.*',
        #r'.*_new_resistance_.*',
    ],
    # global_antibiotic_activity: needs _currently_infected and _sum_any_r_ columns
    # (~1820 columns — same footprint as mean_any_r plots)
    'global_antibiotic_activity': [
        r'.*_currently_infected$',
        r'.*_sum_any_r_(?!hospital).*',   # exclude hospital split variants
    ],
}


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
    #r'.*_new_resistance_.*',
    
    # Per-drug hospital columns (redundant with aggregated)
    r'.*_hospital_.*_drug_.*',
]


def get_required_columns(
    all_columns: List[str],
    include_grouped_plots: bool = True,
    include_calibration: bool = True,
    include_detail_plots: bool = False,
    enabled_detail_plots: Optional[List[str]] = None,
) -> List[str]:
    """
    Determine which columns to load based on enabled analysis types.
    
    Args:
        all_columns: List of all available column names in the CSV
        include_grouped_plots: Include columns for grouped figures 1-10
        include_calibration: Include columns for calibration_summary.txt
        include_detail_plots: DEPRECATED - use enabled_detail_plots instead.
                              If True and enabled_detail_plots is None, loads ALL columns.
        enabled_detail_plots: List of specific detail plot names to load columns for.
                              Only columns needed for these plots will be loaded.
                              e.g., ['death_rate_by_bacteria_region', 'incidence_of_infection']
        
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
    
    # Add patterns for specific enabled detail plots
    if enabled_detail_plots:
        for plot_name in enabled_detail_plots:
            if plot_name in DETAIL_PLOT_PATTERNS:
                include_patterns.extend(DETAIL_PLOT_PATTERNS[plot_name])
    
    # DEPRECATED: If include_detail_plots=True but no specific plots given,
    # fall back to loading everything for backward compatibility
    if include_detail_plots and not enabled_detail_plots:
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
