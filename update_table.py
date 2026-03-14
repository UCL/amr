import re

with open('MODEL_DESCRIPTION.md', 'r', encoding='utf-8') as f:
    text = f.read()

old_table = """| Clinical Setup Variable | Default | Net Mechanism |
|----------|---------|-------------|
| ntibiotic_initiation_base_log_odds | -5.5 | Baseline likelihood of accidental/non-symptomatic prescribing (~0.4% baseline risk). |
| ntibiotic_initiation_log_odds_symptomatic_infection | +6.0 | Very powerful driver escalating probability from ~0.4% to roughly ~62% purely off symptoms breaching threshold. |
| ntibiotic_initiation_log_odds_sepsis | +6.0 | Further powerful boost representing emergent, life-saving care. |
| ntibiotic_initiation_log_odds_immunodeficiency | +2.08 | (ln 8) Heavy prophylactic bias to initiate regimens even tentatively due to patient vulnerability. |
| ntibiotic_initiation_log_odds_no_indication | -1.05 | (ln 0.35) Protective penalty. dampens probability of inappropriate prescribing if patient doesn't functionally have active infection. |
| ntibiotic_initiation_log_odds_test_identified | +0.92 | (ln 2.5) Lab confirmations prompt targeted initiation. |
| ntibiotic_initiation_log_odds_already_on_drug | +0.18 | (ln 1.2) Modest boost representing combined/layered therapy logic once a patient is already functionally linked to a pharmacy loop. |"""

new_table = """| Clinical Setup Variable | Default (Log-odds ratio/modifier) | Meaning |
|----------|---------|-------------|
| ntibiotic_initiation_base_log_odds | -5.5 | Baseline intercept log-odds; likelihood of non-indicated prescribing (~0.4% baseline risk without modifiers). |
| ntibiotic_initiation_log_odds_symptomatic_infection | +6.0 | Escalates probability from ~0.4% to roughly ~62% when symptoms breach the clinical threshold. |
| ntibiotic_initiation_log_odds_sepsis | +6.0 | Additional log-odds multiplier representing emergent, life-saving care escalation. |
| ntibiotic_initiation_log_odds_immunodeficiency | +2.08 | Reflects proactive prophylactic bias (odds ratio ~8) to initiate regimens for vulnerable patients. |
| ntibiotic_initiation_log_odds_no_indication | -1.05 | Negative log-odds modifier (odds ratio ~0.35) dampening inappropriate prescribing without active infection. |
| ntibiotic_initiation_log_odds_test_identified | +0.92 | Log-odds addition (odds ratio ~2.5) when laboratory conformations prompt targeted initiation. |
| ntibiotic_initiation_log_odds_already_on_drug | +0.18 | Minor modifier (odds ratio ~1.2) representing layered therapy logic when a patient is already treated. |"""

text = text.replace(old_table, new_table)

with open('MODEL_DESCRIPTION.md', 'w', encoding='utf-8') as f:
    f.write(text)

print("Markdown updated!")

