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
BACTERIUM_INDEX = 0

# Optional: set this to a filename to save the output (e.g., 'output.txt'). Leave as None to disable.
OUTPUT_FILENAME = 'individual_output.txt'
# Max length for cell display
MAX_CELL_LEN = 20
 
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

    def extract_bacterium(cell):
        # For array columns, return only the value for the selected bacteria
        if ";" in cell:
            parts = cell.split(";")
            if BACTERIUM_INDEX < len(parts):
                return parts[BACTERIUM_INDEX]
            else:
                return ""
        return cell

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
        if len(row) > 2 and row[2] == individual_id:
            # Print all variables as before
            for i, (var, cell) in enumerate(zip(header, row)):
                if i in array_col_indices:
                    val = extract_bacterium(cell)
                else:
                    val = cell
                if len(str(val)) > MAX_CELL_LEN:
                    val = str(val)[:MAX_CELL_LEN-3] + "..."
                line = f"{var}: {val}"
                output_lines.append(line)
            # For each resistance variable, print each drug as its own line
            for res_var in resistance_vars:
                # Find the column for this resistance variable (should be an array column)
                res_indices = [i for i, var in enumerate(header) if res_var in var and i in array_col_indices]
                if res_indices:
                    # For each, extract the value for the selected bacterium (semicolon-separated string for all drugs)
                    for idx in res_indices:
                        cell = row[idx]
                        values = cell.split(";")
                        for drug, value in zip(DRUG_SHORT_NAMES, values):
                            output_lines.append(f"{res_var}_{drug}: {value}")
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
