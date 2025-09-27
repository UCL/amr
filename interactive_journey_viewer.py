"""
Interactive Infection Journey Analysis Tool

Simple interface for exploring individual patient journeys from the AMR simulation.
This script provides easy access to the comprehensive analysis functions.
All figures are saved to files instead of displaying popup windows.
"""

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for saving files only
from infection_journey_analysis import InfectionJourneyAnalyzer
import pandas as pd

def interactive_analysis():
    """Interactive analysis session."""
    print("Loading Infection Journey Analyzer...")
    analyzer = InfectionJourneyAnalyzer()
    
    print(f"\n{'='*60}")
    print("INTERACTIVE INFECTION JOURNEY ANALYSIS")
    print(f"{'='*60}")
    
    # Show quick stats
    total_journeys = analyzer.df['journey_id'].nunique()
    total_individuals = analyzer.df['individual_id'].nunique()
    completed_journeys = len(analyzer.journey_summary[analyzer.journey_summary['resolution_type'] != 'Ongoing'])
    
    print(f"Dataset: {total_journeys} journeys from {total_individuals} individuals")
    print(f"Completed journeys: {completed_journeys}")
    
    while True:
        print(f"\n{'='*40}")
        print("Analysis Options:")
        print("1. View individual journey")
        print("2. Show summary statistics")
        print("3. Plot infection states over time")
        print("4. Calculate 30-day mortality")
        print("5. List journeys by criteria")
        print("6. Show random journey example")
        print("7. Generate detailed text report")
        print("0. Exit")
        
        try:
            choice = input("\nEnter your choice (0-6): ").strip()
            
            if choice == '0':
                print("Goodbye!")
                break
                
            elif choice == '1':
                print("\nAvailable journey IDs:", sorted(analyzer.df['journey_id'].unique())[:20], "...")
                try:
                    journey_id = int(input("Enter journey ID: "))
                    analyzer.view_individual_journey(journey_id=journey_id)
                except (ValueError, KeyError):
                    print("Invalid journey ID. Please try again.")
                    
            elif choice == '2':
                analyzer.generate_summary_report()
                
            elif choice == '3':
                try:
                    max_days = int(input("Enter maximum days to show (default 30): ") or "30")
                    analyzer.plot_infection_states_over_time(max_days=max_days)
                except ValueError:
                    analyzer.plot_infection_states_over_time()
                    
            elif choice == '4':
                mortality_analysis = analyzer.calculate_30_day_mortality()
                if mortality_analysis is not None:
                    print(f"Mortality analysis complete. {len(mortality_analysis)} journeys analyzed.")
                    
            elif choice == '5':
                print("\nFilter options:")
                print("a. By bacteria type")
                print("b. By resolution type") 
                print("c. By region")
                print("d. By age group")
                
                filter_choice = input("Enter filter (a-d): ").strip().lower()
                
                if filter_choice == 'a':
                    bacteria_types = analyzer.journey_summary['primary_bacteria'].unique()
                    print(f"\nAvailable bacteria: {', '.join(bacteria_types[:10])}...")
                    bacteria = input("Enter bacteria name: ").strip()
                    filtered = analyzer.journey_summary[analyzer.journey_summary['primary_bacteria'].str.contains(bacteria, case=False, na=False)]
                    print(f"\nFound {len(filtered)} journeys with {bacteria}:")
                    print(filtered[['journey_id', 'individual_id', 'duration_days', 'resolution_type']].head(10))
                    
                elif filter_choice == 'b':
                    resolution_types = analyzer.journey_summary['resolution_type'].unique()
                    print(f"\nAvailable resolution types: {', '.join(resolution_types)}")
                    resolution = input("Enter resolution type: ").strip()
                    filtered = analyzer.journey_summary[analyzer.journey_summary['resolution_type'] == resolution]
                    print(f"\nFound {len(filtered)} journeys with resolution '{resolution}':")
                    print(filtered[['journey_id', 'individual_id', 'primary_bacteria', 'duration_days']].head(10))
                    
                elif filter_choice == 'c':
                    regions = analyzer.journey_summary['region'].unique()
                    print(f"\nAvailable regions: {', '.join(regions)}")
                    region = input("Enter region: ").strip()
                    filtered = analyzer.journey_summary[analyzer.journey_summary['region'] == region]
                    print(f"\nFound {len(filtered)} journeys from {region}:")
                    print(filtered[['journey_id', 'individual_id', 'primary_bacteria', 'resolution_type']].head(10))
                    
                elif filter_choice == 'd':
                    print("\nAge groups: child (<18), adult (18-65), elderly (>65)")
                    age_group = input("Enter age group: ").strip().lower()
                    age_years = analyzer.journey_summary['age_at_onset'] / 365.25
                    
                    if age_group == 'child':
                        filtered = analyzer.journey_summary[age_years < 18]
                    elif age_group == 'adult':
                        filtered = analyzer.journey_summary[(age_years >= 18) & (age_years <= 65)]
                    elif age_group == 'elderly':
                        filtered = analyzer.journey_summary[age_years > 65]
                    else:
                        print("Invalid age group.")
                        continue
                        
                    print(f"\nFound {len(filtered)} journeys in {age_group} age group:")
                    print(filtered[['journey_id', 'individual_id', 'age_at_onset', 'primary_bacteria']].head(10))
                    
            elif choice == '6':
                # Show random journey
                import random
                random_id = random.choice(list(analyzer.df['journey_id'].unique()))
                print(f"\nShowing random journey example (ID: {random_id}):")
                analyzer.view_individual_journey(journey_id=random_id)
                
            elif choice == '7':
                # Generate detailed text report
                print("\nGenerating comprehensive text report for all journeys...")
                text_report_file = analyzer.generate_text_report()
                print(f"Text report generated: {text_report_file}")
                
            else:
                print("Invalid choice. Please try again.")
                
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            continue
        except Exception as e:
            print(f"Error: {e}")
            continue

if __name__ == "__main__":
    interactive_analysis()