# this reads in the individuals_Log.csv that is optionally created in simulation.rs 
# which contains variable values at each time step for individuals 0 - 9 for
# de-bugging. this allows a certain bacteria to be selected and a certain 
# individual amongst the 10.   


import csv
import sys
from pathlib import Path


# Set this to the index of the bacteria you want to focus on (e.g., 0 for the first bacteria)
'''
BACTERIA INDEX REFERENCE (from population.rs BACTERIA_LIST):
 0: acinetobacter baumannii
 1: citrobacter spp.
 2: enterobacter spp.
 3: enterococcus faecalis
 4: enterococcus faecium
 5: escherichia coli
 6: klebsiella pneumoniae
 7: morganella spp.
 8: proteus spp.
 9: serratia spp.
10: pseudomonas aeruginosa
11: staphylococcus aureus
12: streptococcus pneumoniae
13: salmonella enterica serovar typhi
14: salmonella enterica serovar paratyphi a
15: invasive non-typhoidal salmonella spp.
16: shigella spp.
17: neisseria gonorrhoeae
18: streptococcus pyogenes
19: streptococcus agalactiae
20: haemophilus influenzae
21: chlamydia trachomatis
22: vibrio cholerae
23: neisseria_meningitidis
24: listeria_monocytogenes
25: clostridioides_difficile
26: campylobacter_jejuni
27: enterobacter_cloacae
28: yersinia_enterocolitica
29: moraxella_catarrhalis
30: treponema pallidum
'''

# Set this to the index of the individual you want to print (1 = first individual after header, 2 = second, etc.)
INDIVIDUAL_INDEX = 1 
BACTERIUM_INDEX = -1  # Set to -1 to show all bacteria, or 0-30 for specific bacteria

OUTPUT_FILENAME = 'individual_output.txt'
# Max length for cell display
MAX_CELL_LEN = 20


# Specify the range of time steps to print (inclusive). Set to None to print all.
# Example: TIME_STEP_START = 10; TIME_STEP_END = 20
TIME_STEP_START = None  # e.g., 10
TIME_STEP_END = None    # e.g., 20

# Only print rows where the person is infected with the selected bacterium
ONLY_WHEN_INFECTED = True  # Set to True to enable

# If True, check for ANY bacterial infection; if False, check only the selected bacterium
ANY_BACTERIAL_INFECTION = True  # Set to True to check for any infection, False for specific bacterium

# Variable output toggles - control which variables to print
PRINT_BASIC_INFO = True        # time_step, individual_id, age, etc.
PRINT_INFECTION_LEVELS = True  # level (bacteria infection levels)
PRINT_DRUG_USAGE = True        # cur_use_drug, cur_level_drug
PRINT_RESISTANCE_DATA = True   # any_r, majority_r, activity_r, test_r, microbiome_r, active_infection_activity_r
PRINT_HEALTH_STATUS = True     # hospital_status, sepsis, immunosuppressed, etc.
PRINT_DEMOGRAPHICS = False     # region, sex, etc.
PRINT_INFECTION_RESOLUTION = True  # infection_resolution_this_timestep
PRINT_OTHER_VARIABLES = True   # All other variables not in the above categories
 
# Usage: python individuals_log_view.py [filename]
def print_aligned_csv(filename, max_rows=30):
    output_lines = []
    with open(filename, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        rows = list(reader)

    if not rows:
        print("No data found.")
        return

    # Identify which columns are array columns (by checking for semicolons in the first data row)
    header = rows[0]
    data_rows = rows[1:]
    # Print only the selected individual's data for all time steps (rows with matching individual index)
    if INDIVIDUAL_INDEX < 1 or INDIVIDUAL_INDEX > len(data_rows):
        print(f"INDIVIDUAL_INDEX {INDIVIDUAL_INDEX} is out of range (1-{len(data_rows)})")
        return

    # Find array columns (those with semicolons in the first data row)
    array_col_indices = []
    for i, cell in enumerate(data_rows[0]):
        if ";" in cell:
            array_col_indices.append(i)

    def extract_bacterium(cell, var_name, current_row=None):
        """Extract bacteria data in a readable format, showing only non-zero values"""
        if ";" not in cell:
            return [f"{var_name}: {cell}"]
        
        # Get bacteria names list
        bacteria_names = [
            "acinetobacter_baumannii", "citrobacter_spp", "enterobacter_spp", 
            "enterococcus_faecalis", "enterococcus_faecium", "escherichia_coli",
            "klebsiella_pneumoniae", "morganella_spp", "proteus_spp", "serratia_spp",
            "pseudomonas_aeruginosa", "staphylococcus_aureus", "streptococcus_pneumoniae",
            "salmonella_enterica_serovar_typhi", "salmonella_enterica_serovar_paratyphi_a",
            "invasive_non_typhoidal_salmonella_spp", "shigella_spp", "neisseria_gonorrhoeae",
            "streptococcus_pyogenes", "streptococcus_agalactiae", "haemophilus_influenzae",
            "chlamydia_trachomatis", "vibrio_cholerae", "neisseria_meningitidis",
            "listeria_monocytogenes", "clostridioides_difficile", "campylobacter_jejuni",
            "enterobacter_cloacae", "yersinia_enterocolitica", "moraxella_catarrhalis",
            "treponema_pallidum"
        ]
        
        parts = cell.split(";")
        results = []
        
        # For infection_resolution_this_timestep, show by bacteria and resolution type
        if var_name.lower() == "infection_resolution_this_timestep":
            # This is a flattened array: bacteria_count * resolution_type_count values
            # Resolution types: ImmuneClearance, DrugAssistedClearance, DeathFromSepsis, DeathFromBackground, DeathFromToxicity
            resolution_types = ["ImmuneClearance", "DrugAssistedClearance", "DeathFromSepsis", "DeathFromBackground", "DeathFromToxicity"]
            num_resolution_types = len(resolution_types)
            
            for bact_idx, bacteria_name in enumerate(bacteria_names):
                for res_idx, res_type in enumerate(resolution_types):
                    flat_idx = bact_idx * num_resolution_types + res_idx
                    if flat_idx < len(parts):
                        try:
                            count = int(parts[flat_idx].strip())
                            if count > 0:  # Only show non-zero resolution counts
                                results.append(f"{bacteria_name} {res_type}: {count}")
                        except ValueError:
                            continue
            
            if not results:
                results.append(f"{var_name}: no resolutions this timestep")
            return results
        
        # For clearance_hazard, only show values for bacteria the person is infected with
        if var_name.lower() == "clearance_hazard" and current_row is not None:
            try:
                level_idx = next(
                    i for i, var in enumerate(header) if var.strip().lower() == "level"
                )
                level_cell = current_row[level_idx]
                if ";" in level_cell:
                    level_parts = level_cell.split(";")
                    # Only show clearance_hazard for bacteria with level > 0
                    for i, value in enumerate(parts):
                        if i < len(bacteria_names) and i < len(level_parts):
                            try:
                                level_val = float(level_parts[i].strip())
                                if level_val > 0:
                                    value_stripped = value.strip()
                                    bacteria_name = bacteria_names[i]
                                    try:
                                        val_float = float(value_stripped)
                                        if val_float < 0.01:
                                            formatted_val = f"{val_float:.6f}"
                                        else:
                                            formatted_val = f"{val_float:.3f}"
                                        results.append(
                                            f"{bacteria_name} {var_name}: {formatted_val}"
                                        )
                                    except ValueError:
                                        if value_stripped:
                                            results.append(
                                                f"{bacteria_name} {var_name}: {value_stripped}"
                                            )
                            except ValueError:
                                continue
                    if not results:
                        results.append(f"{var_name}: no infections")
                    return results
            except StopIteration:
                pass  # No level column found, fall back to normal processing
        
        # Normal processing for all other variables
        # Determine which names to use based on variable type
        if any(drug_var in var_name.lower() for drug_var in ['cur_level_drug', 'cur_use_drug', 'ever_taken_drug']):
            # Use drug names for drug-related variables
            name_list = [
                "sulfanilamide", "penicilling", "ampicillin", "amoxicillin",
                "piperacillin", "ticarcillin", "cephalexin", "cefazolin",
                "cefuroxime", "ceftriaxone", "ceftazidime", "cefepime", "ceftaroline", "meropenem", "imipenem_c",
                "ertapenem", "aztreonam", "erythromycin", "azithromycin", "clarithromycin", "clindamycin",
                "gentamicin", "tobramycin", "amikacin", "ciprofloxacin", "levofloxacin", "moxifloxacin",
                "ofloxacin", "tetracycline", "doxyclycline", "minocycline", "vancomycin", "teicoplanin",
                "linezolid", "tedizolid", "quinu_dalfo", "trim_sulf", "chlorampheni", "nitrofurantoin",
                "retapamulin", "fusidic_a", "metronidazole", "furazolidone",
                "amoxicillin_clavulanate", "piperacillin_tazobactam", "ampicillin_sulbactam", "ticarcillin_clavulanate",
                "ceftazidime_avibactam", "meropenem_vaborbactam", "colistin"
            ]
        else:
            # Use bacteria names for bacteria-related variables
            name_list = bacteria_names
        
        # Show individual items with non-zero values
        for i, value in enumerate(parts):
            if i < len(name_list):
                value_stripped = value.strip().lower()
                item_name = name_list[i]
                
                # Handle boolean values - only show 'true' values
                if value_stripped in ['true', 'false']:
                    if value_stripped == 'true':
                        results.append(f"{item_name} {var_name}: true")
                    # Skip 'false' values - we can infer they're false
                else:
                    # Handle numeric values
                    try:
                        val_float = float(value_stripped)
                        if val_float > 0:  # Only show non-zero values
                            if val_float < 0.01:
                                formatted_val = f"{val_float:.6f}"  # More precision for very small values
                            else:
                                formatted_val = f"{val_float:.3f}"  # 3 decimal places for normal values
                            results.append(f"{item_name} {var_name}: {formatted_val}")
                    except ValueError:
                        # Non-numeric, non-boolean value, show as is if non-empty and not null
                        if value_stripped and value_stripped != 'null':
                            results.append(f"{item_name} {var_name}: {value_stripped}")
        
        # If no results found, show a summary based on variable type
        if not results:
            if any(part.strip().lower() in ['true', 'false'] for part in parts):
                results.append(f"{var_name}: all false")
            else:
                results.append(f"{var_name}: all zero")
            
        return results

    def should_print_variable(var_name):
        """Determine if a variable should be printed based on toggles"""
        var_lower = var_name.lower()
        
        # Basic info
        if var_lower in ['time_step', 'individual_id', 'age']:
            return PRINT_BASIC_INFO
        
        # Infection levels
        if 'level' in var_lower and 'drug' not in var_lower:
            return PRINT_INFECTION_LEVELS
        
        # Drug usage
        if any(x in var_lower for x in ['cur_use_drug', 'cur_level_drug', 'drug']):
            return PRINT_DRUG_USAGE
        
        # Health status
        if any(x in var_lower for x in ['hospital', 'sepsis', 'immunosuppressed', 'mortality', 'death']):
            return PRINT_HEALTH_STATUS
        
        # Demographics
        if any(x in var_lower for x in ['region', 'sex', 'vaccination']):
            return PRINT_DEMOGRAPHICS
        
        # Special case for active_infection_activity_r (single value, not array)
        if var_lower == 'active_infection_activity_r':
            return PRINT_RESISTANCE_DATA
        
        # Resistance data (handled separately but check here too)
        if any(x in var_lower for x in ['any_r', 'majority_r', 'activity_r', 'test_r', 'microbiome_r']):
            return PRINT_RESISTANCE_DATA
        
        # Infection resolution data
        if 'infection_resolution' in var_lower:
            return PRINT_INFECTION_RESOLUTION
        
        # Other variables
        return PRINT_OTHER_VARIABLES

    # Get the individual id for the selected index (assume column 2 is individual id)
    selected_row = data_rows[INDIVIDUAL_INDEX-1]
    individual_id = selected_row[2] if len(selected_row) > 2 else None

    # Print only rows for this individual id
    DRUG_SHORT_NAMES = [
        "sulfanilamide", "penicilling", "ampicillin", "amoxicillin",
        "piperacillin", "ticarcillin", "cephalexin", "cefazolin",
        "cefuroxime", "ceftriaxone", "ceftazidime", "cefepime", "ceftaroline", "meropenem", "imipenem_c",
        "ertapenem", "aztreonam", "erythromycin", "azithromycin", "clarithromycin", "clindamycin",
        "gentamicin", "tobramycin", "amikacin", "ciprofloxacin", "levofloxacin", "moxifloxacin",
        "ofloxacin", "tetracycline", "doxyclycline", "minocycline", "vancomycin", "teicoplanin",
        "linezolid", "tedizolid", "quinu_dalfo", "trim_sulf", "chlorampheni", "nitrofurantoin",
        "retapamulin", "fusidic_a", "metronidazole", "furazolidone",
        "amoxicillin_clavulanate", "piperacillin_tazobactam", "ampicillin_sulbactam", "ticarcillin_clavulanate",
        "ceftazidime_avibactam", "meropenem_vaborbactam", "colistin"
    ]
    resistance_vars = ["any_r", "activity_r", "majority_r", "test_r", "microbiome_r"]
    for row in data_rows:
        # Filter by individual id
        if len(row) > 2 and row[2] == individual_id:
            # Filter by infection status if enabled (use 'level' variable for selected bacterium or any bacterium)
            if ONLY_WHEN_INFECTED:
                try:
                    level_idx = next(i for i, var in enumerate(header) if var.strip().lower() == "level")
                    infection_resolution_idx = next(i for i, var in enumerate(header) if var.strip().lower() == "infection_resolution_this_timestep")
                except StopIteration:
                    level_idx = None
                    infection_resolution_idx = None
                
                if level_idx is not None:
                    cell = row[level_idx]
                    has_current_infection = False
                    has_resolution_this_timestep = False
                    
                    # Check for current infections
                    if ";" in cell:
                        parts = cell.split(";")
                        if ANY_BACTERIAL_INFECTION:
                            # Check if infected with ANY bacteria (any value > 0)
                            has_current_infection = any(float(part.strip()) > 0 for part in parts if part.strip().replace('.','').replace('-','').isdigit())
                        else:
                            # Check only the selected bacterium
                            if BACTERIUM_INDEX < len(parts):
                                val = parts[BACTERIUM_INDEX].strip()
                                try:
                                    val_f = float(val)
                                    has_current_infection = val_f > 0
                                except ValueError:
                                    pass
                    
                    # Check for infection resolution events this timestep
                    if infection_resolution_idx is not None:
                        resolution_cell = row[infection_resolution_idx]
                        if ";" in resolution_cell:
                            resolution_parts = resolution_cell.split(";")
                            # Check if any resolution events occurred (any non-zero values)
                            has_resolution_this_timestep = any(int(part.strip()) > 0 for part in resolution_parts if part.strip().isdigit())
                    
                    # Include row if there's either a current infection OR a resolution event this timestep
                    if not (has_current_infection or has_resolution_this_timestep):
                        continue
            # Filter by time step range if specified (assume time step is in column 0 and is integer)
            if TIME_STEP_START is not None or TIME_STEP_END is not None:
                try:
                    timestep = int(row[0])
                except (ValueError, IndexError):
                    timestep = None
                if TIME_STEP_START is not None and (timestep is None or timestep < TIME_STEP_START):
                    continue
                if TIME_STEP_END is not None and (timestep is None or timestep > TIME_STEP_END):
                    continue
            # Print variables based on toggles
            for i, (var, cell) in enumerate(zip(header, row)):
                if should_print_variable(var):
                    if i in array_col_indices:
                        lines = extract_bacterium(cell, var, row)
                        output_lines.extend(lines)
                    else:
                        val = cell
                        if len(str(val)) > MAX_CELL_LEN:
                            val = str(val)[:MAX_CELL_LEN-3] + "..."
                        line = f"{var}: {val}"
                        output_lines.append(line)
            
            # Print resistance data if enabled
            if PRINT_RESISTANCE_DATA:
                for res_var in resistance_vars:
                    # Find the column for this resistance variable (should be an array column)
                    res_indices = [i for i, var in enumerate(header) if res_var in var and i in array_col_indices]
                    if res_indices:
                        # For each, extract the value for the selected bacterium (semicolon-separated string for all drugs)
                        for idx in res_indices:
                            cell = row[idx]
                            values = cell.split(";")
                            # Only show non-zero resistance values
                            for drug, value in zip(DRUG_SHORT_NAMES, values):
                                value_stripped = value.strip().lower()
                                
                                # Handle boolean values - only show 'true'
                                if value_stripped in ['true', 'false']:
                                    if value_stripped == 'true':
                                        output_lines.append(f"{drug} {res_var}: true")
                                else:
                                    # Handle numeric values
                                    try:
                                        val_float = float(value_stripped)
                                        if val_float > 0:  # Only show non-zero resistance
                                            formatted_val = f"{val_float:.3f}" if val_float >= 0.01 else f"{val_float:.6f}"
                                            output_lines.append(f"{drug} {res_var}: {formatted_val}")
                                    except ValueError:
                                        if value_stripped:  # Non-numeric, non-boolean but non-empty
                                            output_lines.append(f"{drug} {res_var}: {value_stripped}")
            output_lines.append("")

    if OUTPUT_FILENAME:
        with open(OUTPUT_FILENAME, 'w', encoding='utf-8') as outf:
            outf.write("\n".join(output_lines))

if __name__ == "__main__":
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    else:
        # Try to auto-detect a likely log file (now looks for individuals_log*.csv)
        candidates = list(Path('.').glob('individuals_log*.csv'))
        if candidates:
            csv_file = str(candidates[0])
        else:
            print("Usage: python individuals_log_view.py <csv_file>")
            sys.exit(1)
    print_aligned_csv(csv_file)
