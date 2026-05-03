from __future__ import annotations

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
AUDIT_CSV = REPO_ROOT / "potency_audit_matrix.csv"

RESIDUAL_SPECIAL_PATHOGENS = {
    "bacteroides_fragilis",
    "mdr_mycobacterium_tuberculosis",
}

CLEAR_INTRINSIC_INACTIVITY = {
    ("bacteroides_fragilis", "amikacin"),
    ("bacteroides_fragilis", "gentamicin"),
    ("bacteroides_fragilis", "tobramycin"),
    ("bacteroides_fragilis", "azithromycin"),
    ("bacteroides_fragilis", "clarithromycin"),
    ("bacteroides_fragilis", "erythromycin"),
    ("bacteroides_fragilis", "daptomycin"),
    ("bacteroides_fragilis", "dalbavancin"),
    ("bacteroides_fragilis", "linezolid"),
    ("bacteroides_fragilis", "tedizolid"),
    ("bacteroides_fragilis", "teicoplanin"),
    ("bacteroides_fragilis", "vancomycin"),
    ("bacteroides_fragilis", "quinu_dalfo"),
    ("bacteroides_fragilis", "retapamulin"),
    ("bacteroides_fragilis", "fusidic_a"),
    ("mdr_mycobacterium_tuberculosis", "dalbavancin"),
    ("mdr_mycobacterium_tuberculosis", "daptomycin"),
    ("mdr_mycobacterium_tuberculosis", "fidaxomicin"),
    ("mdr_mycobacterium_tuberculosis", "flucloxacillin"),
    ("mdr_mycobacterium_tuberculosis", "linezolid"),
    ("mdr_mycobacterium_tuberculosis", "quinu_dalfo"),
    ("mdr_mycobacterium_tuberculosis", "retapamulin"),
    ("mdr_mycobacterium_tuberculosis", "tedizolid"),
    ("mdr_mycobacterium_tuberculosis", "teicoplanin"),
    ("mdr_mycobacterium_tuberculosis", "vancomycin"),
}

INTRINSIC_NOTES = {
    ("bacteroides_fragilis", "amikacin"): "Aminoglycosides should not be interpreted as meaningful baseline activity for B. fragilis; the current low residual values are best treated as explicit inactivity placeholders.",
    ("bacteroides_fragilis", "gentamicin"): "Aminoglycosides should not be interpreted as meaningful baseline activity for B. fragilis; the current low residual values are best treated as explicit inactivity placeholders.",
    ("bacteroides_fragilis", "tobramycin"): "Aminoglycosides should not be interpreted as meaningful baseline activity for B. fragilis; the current low residual values are best treated as explicit inactivity placeholders.",
    ("bacteroides_fragilis", "azithromycin"): "Macrolides should not be interpreted as meaningful baseline activity for B. fragilis in this matrix; the current residual values look like placeholders rather than curated activity.",
    ("bacteroides_fragilis", "clarithromycin"): "Macrolides should not be interpreted as meaningful baseline activity for B. fragilis in this matrix; the current residual values look like placeholders rather than curated activity.",
    ("bacteroides_fragilis", "erythromycin"): "Macrolides should not be interpreted as meaningful baseline activity for B. fragilis in this matrix; the current residual values look like placeholders rather than curated activity.",
    ("bacteroides_fragilis", "daptomycin"): "Daptomycin is not a meaningful baseline anti-Bacteroides agent; low residual values should be treated as explicit inactivity placeholders.",
    ("bacteroides_fragilis", "dalbavancin"): "Lipoglycopeptides are not meaningful anti-Bacteroides baseline agents; low residual values should be treated as explicit inactivity placeholders.",
    ("bacteroides_fragilis", "linezolid"): "Oxazolidinones should not be interpreted as meaningful baseline activity for B. fragilis in this matrix; low residual values should be treated as explicit inactivity placeholders.",
    ("bacteroides_fragilis", "tedizolid"): "Oxazolidinones should not be interpreted as meaningful baseline activity for B. fragilis in this matrix; low residual values should be treated as explicit inactivity placeholders.",
    ("bacteroides_fragilis", "teicoplanin"): "Glycopeptides are not meaningful anti-Bacteroides baseline agents; low residual values should be treated as explicit inactivity placeholders.",
    ("bacteroides_fragilis", "vancomycin"): "Glycopeptides are not meaningful anti-Bacteroides baseline agents; low residual values should be treated as explicit inactivity placeholders.",
    ("bacteroides_fragilis", "quinu_dalfo"): "Streptogramins should not be interpreted as meaningful baseline activity for B. fragilis in this matrix; residual values look inherited rather than curated.",
    ("bacteroides_fragilis", "retapamulin"): "Pleuromutilins are not meaningful baseline agents for B. fragilis in this matrix; residual values should not be interpreted as real activity.",
    ("bacteroides_fragilis", "fusidic_a"): "Fusidic acid should not be interpreted as meaningful baseline activity for B. fragilis in this matrix; residual values look inherited rather than curated.",
    ("mdr_mycobacterium_tuberculosis", "dalbavancin"): "Lipoglycopeptides are not meaningful baseline agents for MDR tuberculosis in this matrix; low residual values should be treated as explicit inactivity placeholders.",
    ("mdr_mycobacterium_tuberculosis", "daptomycin"): "Daptomycin is not a meaningful baseline anti-tubercular agent; low residual values should be treated as explicit inactivity placeholders.",
    ("mdr_mycobacterium_tuberculosis", "fidaxomicin"): "Fidaxomicin is not a meaningful anti-tubercular agent in this matrix; low residual values should not be interpreted as real activity.",
    ("mdr_mycobacterium_tuberculosis", "flucloxacillin"): "Anti-staphylococcal penicillin activity should not be carried into MDR tuberculosis; low residual values are best treated as explicit inactivity placeholders.",
    ("mdr_mycobacterium_tuberculosis", "linezolid"): "The current low linezolid row for MDR tuberculosis is suspicious enough to flag as an explicit cleanup/review target rather than interpret as a meaningful settled potency.",
    ("mdr_mycobacterium_tuberculosis", "quinu_dalfo"): "Streptogramins are not meaningful baseline agents for MDR tuberculosis in this matrix; residual values look inherited rather than curated.",
    ("mdr_mycobacterium_tuberculosis", "retapamulin"): "Pleuromutilins are not meaningful baseline agents for MDR tuberculosis in this matrix; residual values should not be interpreted as real activity.",
    ("mdr_mycobacterium_tuberculosis", "tedizolid"): "The current low tedizolid-style row for MDR tuberculosis is suspicious enough to flag as an explicit cleanup/review target rather than interpret as a meaningful settled potency.",
    ("mdr_mycobacterium_tuberculosis", "teicoplanin"): "Glycopeptides are not meaningful baseline agents for MDR tuberculosis in this matrix; residual values indicate incomplete explicit normalization.",
    ("mdr_mycobacterium_tuberculosis", "vancomycin"): "Glycopeptides are not meaningful baseline agents for MDR tuberculosis in this matrix; residual values indicate incomplete explicit normalization.",
}

DRUG_SPECIFIC_REVIEW_NOTES = {
    "bacteroides_fragilis": {
        "metronidazole": "Metronidazole is a core review anchor for B. fragilis and the current very high value is directionally plausible, but it should still be retained as an explicit evidence-reviewed assignment rather than assumed by default.",
        "meropenem": "Carbapenem activity in B. fragilis is a high-value review target and the current very high value should be retained only after explicit evidence review.",
        "amoxicillin_clavulanate": "BL/BLI activity in B. fragilis is directionally plausible here, but the species-specific value should be checked explicitly rather than accepted from broad class logic.",
        "ampicillin_sulbactam": "BL/BLI activity in B. fragilis is directionally plausible here, but the species-specific value should be checked explicitly rather than accepted from broad class logic.",
        "piperacillin_tazobactam": "BL/BLI activity in B. fragilis is directionally plausible here, but the species-specific value should be checked explicitly rather than accepted from broad class logic.",
        "chloramphenicol": "Chloramphenicol retains substantial activity here and should be reviewed explicitly rather than accepted as a broad historical carryover.",
        "moxifloxacin": "Fluoroquinolone activity in B. fragilis should be reviewed explicitly rather than accepted from broad class inheritance.",
        "ciprofloxacin": "Fluoroquinolone activity in B. fragilis should be reviewed explicitly rather than accepted from broad class inheritance.",
        "levofloxacin": "Fluoroquinolone activity in B. fragilis should be reviewed explicitly rather than accepted from broad class inheritance.",
        "doxycycline": "Tetracycline-family activity in B. fragilis should be reviewed explicitly rather than accepted from broad class inheritance.",
        "minocycline": "Tetracycline-family activity in B. fragilis should be reviewed explicitly rather than accepted from broad class inheritance.",
        "tetracycline": "Tetracycline-family activity in B. fragilis should be reviewed explicitly rather than accepted from broad class inheritance.",
        "trim_sulf": "Trim-sulf activity in B. fragilis should be reviewed explicitly rather than accepted from broad inherited values.",
        "cefiderocol": "The current intermediate cefiderocol value in B. fragilis is suspicious and should be reviewed explicitly rather than accepted as a generic Gram-negative inheritance artifact.",
        "cefixime": "The current cefixime value in B. fragilis should be reviewed explicitly rather than accepted from broad cephalosporin inheritance.",
    },
    "mdr_mycobacterium_tuberculosis": {
        "rifampicin": "Rifampicin is a defining drug for tuberculosis, so the current very low MDR-tuberculosis value is directionally plausible but should remain an explicit evidence-reviewed assignment rather than a casual low fallback.",
        "moxifloxacin": "Fluoroquinolone activity in MDR tuberculosis is a high-value review target and should be checked explicitly rather than accepted from broad class inheritance.",
        "levofloxacin": "Fluoroquinolone activity in MDR tuberculosis is a high-value review target and should be checked explicitly rather than accepted from broad class inheritance.",
        "ciprofloxacin": "Fluoroquinolone activity in MDR tuberculosis is a high-value review target and should be checked explicitly rather than accepted from broad class inheritance.",
        "amikacin": "Aminoglycoside activity in MDR tuberculosis should be reviewed explicitly rather than accepted from broad inherited values.",
        "gentamicin": "Aminoglycoside activity in MDR tuberculosis should be reviewed explicitly rather than accepted from broad inherited values.",
        "tobramycin": "Aminoglycoside activity in MDR tuberculosis should be reviewed explicitly rather than accepted from broad inherited values.",
        "linezolid": "Linezolid is a high-value review target in MDR tuberculosis and the current very low row should not be accepted without explicit evidence review.",
        "clarithromycin": "Macrolide activity in MDR tuberculosis should be reviewed explicitly rather than accepted from broad inherited values.",
        "azithromycin": "Macrolide activity in MDR tuberculosis should be reviewed explicitly rather than accepted from broad inherited values.",
        "erythromycin": "Macrolide activity in MDR tuberculosis should be reviewed explicitly rather than accepted from broad inherited values.",
        "doxycycline": "Tetracycline-family activity in MDR tuberculosis should be reviewed explicitly rather than accepted from broad inherited values.",
        "minocycline": "Tetracycline-family activity in MDR tuberculosis should be reviewed explicitly rather than accepted from broad inherited values.",
        "tetracycline": "Tetracycline-family activity in MDR tuberculosis should be reviewed explicitly rather than accepted from broad inherited values.",
        "meropenem": "Carbapenem activity in MDR tuberculosis is a high-value review target and should be checked explicitly rather than accepted from broad class inheritance.",
        "chloramphenicol": "Chloramphenicol activity in MDR tuberculosis should be reviewed explicitly rather than accepted from broad inherited values.",
        "trim_sulf": "Trim-sulf activity in MDR tuberculosis should be reviewed explicitly rather than accepted from broad inherited values.",
        "metronidazole": "Metronidazole activity in MDR tuberculosis should be reviewed explicitly rather than accepted from low residual inheritance.",
        "cefiderocol": "The current low cefiderocol value in MDR tuberculosis is suspicious and should be reviewed explicitly rather than accepted as a generic fallback.",
        "cefixime": "The current cefixime value in MDR tuberculosis should be reviewed explicitly rather than accepted from broad cephalosporin inheritance.",
    },
}


def annotate_row(row: dict[str, str]) -> None:
    bacteria = row["bacteria"]
    drug = row["drug"]
    potency = float(row["potency_no_r"])

    if bacteria not in RESIDUAL_SPECIAL_PATHOGENS:
        return

    key = (bacteria, drug)
    if key in CLEAR_INTRINSIC_INACTIVITY:
        row["review_status"] = "clear intrinsic inactivity"
        row["evidence_notes"] = INTRINSIC_NOTES[key]
        row["decision"] = "later normalize to explicit inactivity after full review"
        return

    organism_notes = DRUG_SPECIFIC_REVIEW_NOTES.get(bacteria, {})
    if drug in organism_notes:
        row["review_status"] = "needs species-specific review"
        row["evidence_notes"] = organism_notes[drug]
        row["decision"] = "defer potency change until evidence review"
        return

    if potency <= 0.1:
        row["review_status"] = "fallback default; needs explicit review"
        row["evidence_notes"] = (
            "Residual special-pathogen first-pass review: this row remains at a very low residual potency and should be treated as an unreviewed fallback until explicitly adjudicated."
        )
        row["decision"] = "defer potency change until evidence review"
        return

    row["review_status"] = "needs species-specific review"
    row["evidence_notes"] = (
        "Residual special-pathogen first-pass review: retained activity should be reviewed explicitly for this organism rather than accepted from broad class inheritance."
    )
    row["decision"] = "defer potency change until evidence review"


def main() -> None:
    with AUDIT_CSV.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        fieldnames = list(rows[0].keys())

    for row in rows:
        annotate_row(row)

    with AUDIT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    reviewed = sum(
        1 for row in rows if row["bacteria"] in RESIDUAL_SPECIAL_PATHOGENS and row["review_status"]
    )
    print(f"Annotated {reviewed} residual special-pathogen rows in {AUDIT_CSV.name}")


if __name__ == "__main__":
    main()