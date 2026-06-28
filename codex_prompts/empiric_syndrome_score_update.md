You are working in the UCL/amr repository.

Goal:
Revise the empiric active-infection drug scoring templates in src/config.rs so they are less dominated by extreme carbapenem and BL/BLI scores, while still keeping carbapenems and BL/BLI agents meaningfully represented because previous less-extreme templates substantially undershot these classes.

Scope constraints:
- Edit src/config.rs.
- Do not change Rust treatment-selection logic.
- Do not change drug-initiation probabilities.
- Do not change reserve-drug gate logic.
- Do not change drug-selection temperature.
- Do not change syndrome 0 / no-active-modelled-bacterial-infection scores in this task.
- Do not change drug potency, resistance, mortality, sepsis, or calibration scoring logic.
- Do not add new parameter families.
- Use the existing `syndrome_<id>_empiric_drug_<drug>_score` parameter mechanism.
- Update comments in src/config.rs so they no longer describe values as “aggressive calibration” where the values have been moderated.
- Treat the current code as the source of truth and verify all drug names against DRUG_SHORT_NAMES before editing.

Files to inspect:
- @src/config.rs
- @src/rules/mod.rs
- @tests/config_invariants.rs
- @tests/probability_invariants.rs
- @tests/dimension_invariants.rs
- @tests/runtime_config_validation.rs

Context:
The current empiric active-infection syndrome templates contain very high carbapenem and BL/BLI scores. These were introduced to address previous under-use of beta-lactamase inhibitor combinations and carbapenems. However, the latest drug-class calibration now shows:
- Penicillins too low;
- BL/BLI somewhat high;
- fluoroquinolones high;
- macrolides high;
- carbapenems slightly high;
- tetracyclines, nitrofurans, fosfomycin, aminoglycosides, glycopeptides, oxazolidinones/lipoglycopeptides, nitroimidazoles low.

Important:
Do not revert carbapenems to the old very low values. Instead, reduce extreme values to a moderate high range. Carbapenems should remain strong options in severe syndromes, especially bloodstream, intra-abdominal, bone/joint, and other severe/device infections, but should no longer have raw scores in the hundreds-to-thousands that overwhelm other plausible empiric choices.

Also reduce BL/BLI scores in common non-severe/high-incidence syndromes, but do not erase them. BL/BLI agents should remain plausible empiric choices, especially in severe or mixed-flora syndromes.

Specific requested changes:
In the `empiric_syndrome_templates` block, update the following values.

1. Syndrome 1 — UTI / Genitourinary

Change to approximately:

- trim_sulf: 1.0
- amoxicillin_clavulanate: 16.0
- amoxicillin: 3.0
- ciprofloxacin: 1.5
- ampicillin: 1.5
- levofloxacin: 1.5
- nitrofurantoin: 15.0
- fosfomycin: 13.0
- cephalexin: 10.0
- ceftriaxone: 2.5
- cefazolin: 8.0
- cefuroxime: 8.5
- piperacillin_tazobactam: 15.0
- cefepime: 4.0
- ceftazidime: 4.0
- meropenem: 25.0
- imipenem_c: 25.0
- ertapenem: 25.0
- meropenem_vaborbactam: 18.0
- ceftazidime_avibactam: 12.0
- aztreonam_avibactam: 1.0
- cefixime: 7.0
- colistin: 0.4
- vancomycin: 0.3
- linezolid: 0.1

Rationale:
Keep nitrofurantoin/fosfomycin prominent; retain BL/BLI and carbapenem availability but reduce extreme broad-spectrum dominance; reduce fluoroquinolone contribution.

2. Syndrome 2 — Skin / soft tissue

Change to approximately:

- flucloxacillin: 12.0
- amoxicillin_clavulanate: 12.0
- amoxicillin: 6.0
- cephalexin: 14.0
- ampicillin: 5.0
- penicillin_g: 5.0
- cefazolin: 14.0
- clindamycin: 4.0
- trim_sulf: 0.5
- doxycycline: 3.5
- minocycline: 3.0
- linezolid: 10.0
- tedizolid: 9.0
- dalbavancin: 9.0
- vancomycin: 13.0
- quinu_dalfo: 8.0
- rifampicin: 0.5
- ciprofloxacin: 1.5
- piperacillin_tazobactam: 12.0

Rationale:
Restore anti-staphylococcal/Access beta-lactams as plausible SSTI empiric choices; keep MRSA agents for severe contexts; reduce piperacillin-tazobactam dominance.

3. Syndrome 3 — Respiratory

Change to approximately:

- amoxicillin_clavulanate: 18.0
- amoxicillin: 8.0
- penicillin_g: 4.0
- ampicillin: 4.0
- azithromycin: 2.5
- clarithromycin: 2.5
- ceftriaxone: 9.0
- erythromycin: 3.0
- cefuroxime: 8.5
- piperacillin_tazobactam: 16.0
- levofloxacin: 2.0
- moxifloxacin: 2.0
- cefixime: 6.5
- aztreonam_avibactam: 0.01
- cefepime: 7.5
- cephalexin: 7.0
- doxycycline: 5.0
- vancomycin: 8.0
- meropenem: 28.0
- imipenem_c: 28.0
- ofloxacin: 2.0
- linezolid: 7.0
- minocycline: 3.5

Rationale:
Move respiratory prescribing back toward amoxicillin/doxycycline/macrolide/ceftriaxone balance; keep BL/BLI and carbapenems present for severe/hospital respiratory infection; reduce fluoroquinolones and macrolides modestly.

4. Syndrome 4 — Bloodstream / bacteremia

Change to approximately:

- piperacillin_tazobactam: 24.0
- meropenem: 180.0
- imipenem_c: 140.0
- meropenem_vaborbactam: 100.0
- ceftazidime_avibactam: 22.0
- aztreonam_avibactam: 2.0
- cefepime: 14.0
- ceftazidime: 12.0
- ceftriaxone: 10.0
- ampicillin_sulbactam: 18.0
- amoxicillin_clavulanate: 16.0
- ampicillin: 8.0
- amoxicillin: 8.0
- penicillin_g: 4.0
- flucloxacillin: 8.0
- vancomycin: 14.0
- linezolid: 10.0
- tedizolid: 9.0
- dalbavancin: 8.0
- quinu_dalfo: 8.5
- gentamicin: 12.0
- tobramycin: 11.0
- amikacin: 12.0
- colistin: 0.3
- cefazolin: 7.0
- ciprofloxacin: 2.0
- levofloxacin: 2.0
- cephalexin: 3.0
- rifampicin: 0.5

Rationale:
Bloodstream infection should still support strong carbapenem use. Do not revert to low values. But reduce 1200/900/700 to high-but-interpretable values so piperacillin-tazobactam, cefepime, ceftriaxone, vancomycin, and aminoglycosides can still compete.

5. Syndrome 5 — Intra-abdominal

Change to approximately:

- metronidazole: 8.0
- piperacillin_tazobactam: 24.0
- ampicillin_sulbactam: 16.0
- amoxicillin_clavulanate: 16.0
- meropenem: 120.0
- imipenem_c: 100.0
- ertapenem: 80.0
- ceftazidime: 9.0
- cefepime: 9.0
- ceftriaxone: 9.0
- ceftazidime_avibactam: 18.0
- aztreonam_avibactam: 2.0
- meropenem_vaborbactam: 75.0
- ciprofloxacin: 1.5
- levofloxacin: 1.5
- ampicillin: 5.0
- amoxicillin: 5.0
- trim_sulf: 0.1
- gentamicin: 8.0
- amikacin: 8.0
- colistin: 0.5

Rationale:
Keep piperacillin-tazobactam and carbapenems prominent for severe intra-abdominal infection, but reduce extreme 650–900 values. Increase metronidazole because nitroimidazole use is low, while recognising that the model may not fully represent combination therapy.

6. Syndrome 6 — Central nervous system

Change to approximately:

- ceftriaxone: 15.0
- ceftazidime: 12.0
- cefepime: 12.0
- penicillin_g: 11.0
- ampicillin: 10.0
- vancomycin: 13.0
- linezolid: 10.0
- cefixime: 1.0
- meropenem: 35.0
- imipenem_c: 25.0
- chloramphenicol: 2.0
- rifampicin: 1.0
- piperacillin_tazobactam: 6.0

Rationale:
CNS therapy should centre on ceftriaxone/ampicillin/penicillin/vancomycin/meropenem. Piperacillin-tazobactam should not be a major CNS empiric competitor.

7. Syndrome 7 — Gastrointestinal, non-invasive

Change to approximately:

- ciprofloxacin: 1.5
- azithromycin: 3.0
- amoxicillin_clavulanate: 10.0
- amoxicillin: 4.0
- ampicillin: 4.0
- levofloxacin: 1.5
- ampicillin_sulbactam: 8.0
- trim_sulf: 0.5
- doxycycline: 3.5
- minocycline: 2.5
- cefixime: 4.0
- penicillin_g: 2.0
- cephalexin: 4.0
- cefuroxime: 4.0
- furazolidone: 0.2
- metronidazole: 2.5
- rifampicin: 0.5

Rationale:
GI is a high-incidence syndrome, so avoid using it as a BL/BLI calibration lever. Keep antibiotics less dominant and more focused on selected agents; increase metronidazole modestly; reduce fluoroquinolones.

8. Syndrome 8 — Genital / pelvic

Change to approximately:

- penicillin_g: 4.0
- azithromycin: 2.5
- ceftriaxone: 13.0
- cefixime: 10.5
- doxycycline: 6.5
- amoxicillin_clavulanate: 10.0
- amoxicillin: 4.0
- cefuroxime: 8.0
- clindamycin: 4.0
- ampicillin: 5.0
- ampicillin_sulbactam: 8.0
- ciprofloxacin: 1.0
- levofloxacin: 1.5
- cephalexin: 5.0
- trim_sulf: 0.3
- metronidazole: 4.0
- rifampicin: 0.5

Rationale:
Genital/pelvic empiric therapy should be centred more on ceftriaxone, doxycycline, azithromycin, and metronidazole. Reduce BL/BLI and fluoroquinolone influence.

9. Syndrome 9 — Bone / joint / hardware-associated

Change to approximately:

- flucloxacillin: 12.0
- cefazolin: 15.0
- vancomycin: 14.0
- ampicillin: 8.0
- ceftriaxone: 11.0
- cephalexin: 12.0
- penicillin_g: 6.0
- linezolid: 11.0
- tedizolid: 10.0
- dalbavancin: 10.0
- clindamycin: 2.0
- ciprofloxacin: 2.5
- levofloxacin: 2.5
- trim_sulf: 0.5
- meropenem: 80.0
- piperacillin_tazobactam: 15.0
- rifampicin: 6.0

Rationale:
Keep carbapenems available for severe resistant Gram-negative bone/joint infection, but reduce 700 substantially. Restore anti-staphylococcal beta-lactams as central empiric choices.

10. Syndrome 10 — Other severe / device-related catch-all

Change to approximately:

- vancomycin: 14.0
- linezolid: 10.0
- piperacillin_tazobactam: 24.0
- cefepime: 10.0
- ceftriaxone: 7.0
- meropenem: 120.0
- imipenem_c: 100.0
- aztreonam_avibactam: 2.0
- ciprofloxacin: 2.5
- azithromycin: 2.0

Rationale:
Keep broad severe/device coverage and preserve carbapenem use, but remove extreme 800–900 dominance. Keep vancomycin prominent because device/line infections are often Gram-positive.

Comment cleanup:
- Remove or update comments such as “Aggressive calibration: major empiric carbapenem escalation.”
- Replace with wording like:
  “High severe-syndrome carbapenem score retained to avoid prior carbapenem under-use, but moderated to avoid overwhelming other empiric options.”
- Where BL/BLI scores are moderated, use wording like:
  “BL/BLI retained as a common broad empiric option, but reduced from prior calibration-boosted value.”

Do not change:
- syndrome 0 background/no-active-modelled-bacterial-infection scores;
- prophylaxis scores;
- targeted therapy scoring;
- reserve-drug gates;
- selection temperature;
- drug potency matrix;
- treatment-failure logic.

Validation:
After editing, run:
cargo fmt --all --check
cargo test --test config_invariants
cargo test --test probability_invariants
cargo test --test dimension_invariants
cargo test --test runtime_config_validation

If a listed test does not exist, run the closest relevant test and report exactly what was run.

Do not run a long simulation unless requested.

Acceptance criteria:
- Empiric active-infection syndrome templates are updated with the approximate values above.
- Carbapenem scores are reduced from extreme 700–1200 values but remain substantially above older low values.
- BL/BLI scores are reduced in high-incidence common syndromes but remain meaningful.
- Access beta-lactams regain more plausible roles in respiratory, skin, bone/joint, and bloodstream syndromes.
- Fluoroquinolone scores are reduced in common syndromes.
- Nitrofurantoin and fosfomycin remain prominent in UTI.
- Metronidazole is increased in intra-abdominal and pelvic/genital contexts.
- No treatment-selection algorithm or model mechanics are changed.
- Final response should include:
  - file changed;
  - full list of syndrome template values changed;
  - confirmation that syndrome 0 was not changed;
  - confirmation that all drug names were checked against DRUG_SHORT_NAMES;
  - tests run;
  - any drug names from the requested list that were absent or skipped.