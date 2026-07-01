You are working in the UCL/amr repository.

Goal:
Propose and, if the relevant parameter file is unambiguous, implement a conservative first-pass adjustment to antibiotic drug-class selection parameters to improve calibration of:

“Figure 3. Calibration: 2025 antibiotic use by drug class”

The latest Figure 3 suggests the model’s simulated 2025 drug-class shares differ from target/benchmark shares. Use the generated calibration data to quantify the mismatch exactly, then update only the relevant drug-use parameter values.

Do not change target values.
Do not change model logic.
Do not change treatment-selection code unless no parameterised control exists.

Visual directional guidance from the latest Figure 3:
Assuming coloured bars are simulated values and hatched bars are targets, the approximate direction appears to be:

- Beta-lactamase combinations (J01CR): close to target, possibly slightly high; leave unchanged or very small decrease.
- Penicillins (J01C): under target; increase.
- Macrolides (J01F): substantially over target; decrease strongly.
- Quinolones/fluoroquinolones (J01M): over target; decrease modestly.
- Cephalosporins 1-2G: over target; decrease.
- Tetracyclines (J01A): close to target; leave unchanged or very small adjustment.
- Cephalosporins 3G: over target; decrease modestly.
- Sulfonamides (J01E): close to target; leave unchanged.
- Nitrofurans (J01XE): under target; increase.
- Carbapenems (J01DH): slightly under target; increase only cautiously, preserving reserve-use logic.

Mandatory first step:
Audit the exact values behind Figure 3 before editing parameters.

Inspect:
- make_paper_tables.py
- parse_calibration.py
- output_graphs/calibration_summary_*.txt
- any simulation_summary_*.csv files used for Figure 3
- parameter/config files controlling antibiotic class choice, empiric treatment, targeted treatment, class weights, drug shares, or calibration targets

Search terms:
- drug share
- antibiotic use
- class share
- J01C
- J01CR
- J01F
- J01M
- J01A
- J01E
- J01XE
- J01DH
- penicillin
- macrolide
- quinolone
- cephalosporin
- nitrofuran
- carbapenem
- empiric
- targeted
- treatment weight
- class_weight
- drug_class_weight
- prescribing

Report the exact current simulated value, target value, absolute difference, and relative ratio for each Figure 3 class.

Parameter adjustment strategy:
Use a damped multiplicative update to avoid overcorrection.

If the model has class selection weights or multipliers, compute candidate updates as:

new_weight = old_weight * clamp((target_share / simulated_share)^0.6, 0.65, 1.35)

Use exact values from the calibration outputs, not the screenshot.

If parameters are probabilities that must sum to 1:
- apply the candidate multipliers;
- renormalise according to the existing config convention;
- preserve any existing constraints or reserve-use caps.

If parameters are independent multipliers:
- apply the multiplier directly;
- do not renormalise unless the code already expects normalised values.

Suggested qualitative bounds:
- Macrolides: substantial decrease, roughly 25-35% if exact ratios support it.
- Cephalosporins 1-2G: decrease roughly 15-25%.
- Quinolones/fluoroquinolones: decrease roughly 10-20%.
- Cephalosporins 3G: decrease roughly 10-20%.
- Penicillins: increase roughly 15-25%.
- Nitrofurans: increase roughly 15-30%.
- Carbapenems: increase cautiously, roughly 5-15%, and only through the existing reserve/severe-infection pathway if such a parameter exists.
- Beta-lactamase combinations: unchanged or decrease no more than 5%.
- Tetracyclines: unchanged unless exact mismatch is meaningful.
- Sulfonamides: unchanged unless exact mismatch is meaningful.

Important:
Do not adjust calibration target values to make the fit look better.
Do not change ATC class mappings.
Do not change drug-class labels.
Do not change treatment initiation rates unless the audit shows the mismatch is caused by total antibiotic exposure rather than class allocation.
Do not change resistance, infection, sepsis, testing, or mortality logic.

If multiple parameter layers exist:
Prefer the most direct class-allocation parameter for 2025 drug-share calibration.

For example:
1. drug class selection weights;
2. empiric class preference weights;
3. syndrome-specific class weights;
4. hospital/community class weights.

Avoid changing broad infection or treatment probabilities unless class-specific controls do not exist.

Special caution:
For carbapenems, do not make them broadly common first-line empiric therapy. If carbapenems are under target, prefer a small increase in existing hospital/severe/Gram-negative/reserve-use pathways rather than a general population-wide increase.

Outputs:
- Update the relevant config/parameter file if the correct parameters are unambiguous.
- Add a short calibration note if the repository has a suitable notes/changelog location.
- Do not create new model-output CSV files.
- Do not add new figures or tables.
- Do not change make_paper_tables.py except if needed only to identify/report the exact values.

Validation:
Run:
- cargo fmt --all --check, if Rust/config formatting is affected and this is applicable.
- cargo test --test runtime_config_validation, if present.
- cargo test --test config_invariants, if present.
- python -m py_compile make_paper_tables.py, only if Python files were touched.

Do not run a long calibration simulation unless the repository has a clearly established quick smoke-test command.
If no new model run is performed, state that Figure 3 improvement requires a new run to confirm.

Final response should include:
- exact Figure 3 simulated and target shares used;
- mismatch direction by class;
- files changed;
- parameter names changed;
- old and new parameter values;
- multiplier applied for each class;
- whether values were renormalised;
- whether carbapenem reserve-use logic was preserved;
- checks run;
- confirmation that target values and model logic were not changed;
- whether a new model run is needed to assess the updated Figure 3 fit.