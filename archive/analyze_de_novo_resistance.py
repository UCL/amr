#!/usr/bin/env python3
"""
Quick analysis script to examine de novo resistance detection results
"""

import pandas as pd
import numpy as np

def analyze_de_novo_resistance():
    print("Loading infection journeys data...")
    
    try:
        # Load the CSV file
        df = pd.read_csv('infection_journeys.csv')
        print(f"Total snapshots: {len(df)}")
        print(f"Unique journeys: {df['journey_id'].nunique()}")
        
        # Check for de novo resistance flags
        de_novo_true = df[df['has_de_novo_resistance'] == True]
        de_novo_false = df[df['has_de_novo_resistance'] == False]
        
        print(f"\nDe novo resistance summary:")
        print(f"Journeys flagged with de novo resistance: {len(de_novo_true)}")
        print(f"Journeys without de novo resistance: {len(de_novo_false)}")
        
        if len(de_novo_true) > 0:
            print(f"\nDE NOVO RESISTANCE DETECTED!")
            print(f"Unique journeys with de novo resistance: {de_novo_true['journey_id'].nunique()}")
            
            # Show details of flagged journeys
            flagged_journeys = de_novo_true['journey_id'].unique()
            print(f"\nFlagged journey IDs: {flagged_journeys[:10]}")  # Show first 10
            
            # Analysis by bacteria type
            bacteria_counts = de_novo_true['primary_bacteria'].value_counts()
            print(f"\nDe novo resistance by bacteria:")
            print(bacteria_counts.head())
            
            # Check time periods
            print(f"\nTime steps with de novo resistance:")
            print(f"Min time step: {de_novo_true['time_step'].min()}")
            print(f"Max time step: {de_novo_true['time_step'].max()}")
            
            # Check if they had active treatment
            treatment_cases = de_novo_true[de_novo_true['current_drugs'] != '']
            print(f"\nCases with active drug treatment: {len(treatment_cases)}")
            
            if len(treatment_cases) > 0:
                print("Sample cases with treatment and de novo resistance:")
                sample = treatment_cases[['journey_id', 'time_step', 'primary_bacteria', 'current_drugs', 'day_of_journey']].head()
                print(sample)
        else:
            print("\nNo de novo resistance detected.")
            
            # Diagnostic info
            print(f"\nDiagnostic information:")
            
            # Check time range
            print(f"Time step range: {df['time_step'].min()} to {df['time_step'].max()}")
            
            # Check for drug treatment
            treated = df[df['current_drugs'] != '']
            print(f"Snapshots with drug treatment: {len(treated)}")
            
            if len(treated) > 0:
                print(f"Treatment time steps: {treated['time_step'].min()} to {treated['time_step'].max()}")
                print("Sample treated cases:")
                sample = treated[['journey_id', 'time_step', 'primary_bacteria', 'current_drugs', 'day_of_journey']].head()
                print(sample)
            
            # Check resistance fields
            resistance_fields = ['resistance_any_r', 'resistance_majority_r', 'resistance_activity_r']
            for field in resistance_fields:
                non_empty = df[df[field] != '']
                print(f"Non-empty {field}: {len(non_empty)}")
        
    except FileNotFoundError:
        print("infection_journeys.csv not found. Make sure the simulation has completed.")
    except Exception as e:
        print(f"Error analyzing data: {e}")

if __name__ == "__main__":
    analyze_de_novo_resistance()