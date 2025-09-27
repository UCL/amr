INFECTION JOURNEY ANALYSIS - FIGURE LOCATIONS
==============================================

The infection journey analysis generates figures and saves them automatically to files.
No popup windows will appear when running the analysis scripts.

📁 MAIN FIGURE LOCATIONS:

1. Journey Analysis Figures:
   Location: output_graphs/journey_analysis/
   
   📊 Population-Level Analysis:
   - infection_states_stacked_first_30_days.png    (Stacked bar chart of daily infection states)
   - 30_day_mortality_analysis.png                  (Mortality analysis with demographic breakdowns)
   
   📈 Individual Patient Journeys (examples):
   - journey_1_individual_16_timeline.png          (70.5yr female, C. difficile, 24 days)
   - journey_4_individual_434_timeline.png         (17.9yr female, Vibrio cholerae)
   - journey_39_individual_1360_timeline.png       (Patient with treatment journey)
   - journey_52_individual_95_timeline.png         (Chronic immunodeficiency case)
   - journey_59_individual_1244_timeline.png       (Drug-assisted clearance)
   - [... and more individual timeline examples]

2. Main Simulation Figures:
   Location: output_graphs/
   
   📊 Population Summary Plots:
   - grouped_figure_1.png through grouped_figure_9.png
   
   📁 Detailed Analysis Subdirectories:
   - age_distribution_by_region/
   - death_rate_by_bacteria_region/
   - drug_usage_ddd_per_1000_per_day/
   - incidence_of_infection/
   - resistance_mechanism_by_bacteria/
   - [... and many more specialized analysis folders]

🔧 USAGE INSTRUCTIONS:

To generate figures:
1. Run comprehensive analysis:    python infection_journey_analysis.py
2. Run interactive viewer:        python interactive_journey_viewer.py
3. Generate more examples:        python generate_journey_examples.py

📋 FIGURE DESCRIPTIONS:

Individual Journey Timelines (4-panel plots):
- Top Left: Primary bacteria level over time
- Top Right: Immunity & toxicity levels
- Bottom Left: Treatment & sepsis status 
- Bottom Right: Hospital status progression

Population Analysis:
- Stacked bars show daily counts of: Not Treated, Treated, Sepsis states
- Mortality analysis shows death rates by bacteria type, sepsis, treatment, age

🎯 KEY FEATURES:

✅ All figures saved automatically as high-resolution PNG files
✅ No popup windows - runs completely in background
✅ Comprehensive day-by-day patient progression data
✅ Population-level statistics and trends
✅ 30-day mortality probability analysis
✅ Interactive filtering and exploration tools available

📧 File Naming Convention:
- journey_[ID]_individual_[ID]_timeline.png  (Individual patient journeys)
- infection_states_stacked_first_[N]_days.png (Population infection states)
- 30_day_mortality_analysis.png (Mortality statistics)

Last updated: 2025-09-27
Analysis covers: 147 unique infection journeys from 103 individuals