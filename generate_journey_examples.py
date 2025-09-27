"""
Generate Multiple Journey Figures

This script creates individual journey timeline figures for multiple patients
to build a collection of example cases. All figures are saved as files.
"""

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for saving files only
from infection_journey_analysis import InfectionJourneyAnalyzer
import random

def generate_multiple_journey_figures(num_examples=10):
    """Generate figures for multiple example journeys."""
    print(f"Generating {num_examples} example journey figures...")
    
    analyzer = InfectionJourneyAnalyzer()
    
    # Get all available journey IDs
    all_journey_ids = list(analyzer.df['journey_id'].unique())
    
    # Select diverse examples
    selected_journeys = []
    
    # Try to get journeys with different characteristics
    for resolution_type in ['ImmuneClearance', 'DrugAssistedClearance', 'Ongoing']:
        matching_journeys = analyzer.journey_summary[
            analyzer.journey_summary['resolution_type'] == resolution_type
        ]['journey_id'].tolist()
        
        if matching_journeys:
            # Get 2-3 examples of each type
            sample_size = min(3, len(matching_journeys))
            selected_journeys.extend(random.sample(matching_journeys, sample_size))
    
    # Fill up to num_examples with random additional journeys
    remaining_journeys = [j for j in all_journey_ids if j not in selected_journeys]
    additional_needed = max(0, num_examples - len(selected_journeys))
    if additional_needed > 0 and remaining_journeys:
        additional_count = min(additional_needed, len(remaining_journeys))
        selected_journeys.extend(random.sample(remaining_journeys, additional_count))
    
    # Generate figures for selected journeys
    print(f"\nGenerating figures for {len(selected_journeys)} journeys...")
    
    for i, journey_id in enumerate(selected_journeys[:num_examples], 1):
        print(f"Processing journey {journey_id} ({i}/{min(num_examples, len(selected_journeys))})")
        try:
            analyzer.view_individual_journey(journey_id=journey_id)
        except Exception as e:
            print(f"Error processing journey {journey_id}: {e}")
    
    print(f"\nFigures saved to: {analyzer.output_dir}")
    print("Individual journey timelines show:")
    print("  - Primary bacteria level over time")
    print("  - Immunity & toxicity levels") 
    print("  - Treatment & sepsis status")
    print("  - Hospital status progression")

if __name__ == "__main__":
    generate_multiple_journey_figures(15)  # Generate 15 example journeys