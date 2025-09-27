#!/usr/bin/env python3
"""
Script to regenerate only the syndrome death rate plots with improved Y-axis scaling.
"""

import sys
import os
sys.path.append('.')

# Import the main analysis script
from amr_analysis import main

# Backup the original config file
import shutil
shutil.copy('amr_simulation_output_analysis/config.py', 'amr_simulation_output_analysis/config.py.backup')

print("Creating minimal config for syndrome plots only...")

# Read the current config
with open('amr_simulation_output_analysis/config.py', 'r') as f:
    config_content = f.read()

# Replace all True values with False except for death_rate_by_syndrome_region
lines = config_content.split('\n')
modified_lines = []

for line in lines:
    if ': bool = True' in line and 'death_rate_by_syndrome_region' not in line:
        # Replace True with False for all boolean configs except our target
        modified_line = line.replace(': bool = True', ': bool = False')
        modified_lines.append(modified_line)
    else:
        modified_lines.append(line)

# Write the modified config
with open('amr_simulation_output_analysis/config.py', 'w') as f:
    f.write('\n'.join(modified_lines))

print("Running analysis with syndrome plots only...")
try:
    # Run the main analysis
    main()
    print("✓ Syndrome plots generated successfully!")
    print("Check output_graphs/death_rate_by_syndrome_region/ for the updated plots")
finally:
    # Restore the original config
    print("Restoring original config...")
    shutil.move('amr_simulation_output_analysis/config.py.backup', 'amr_simulation_output_analysis/config.py')
    print("Done!")