from __future__ import annotations

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
AUDIT_CSV = REPO_ROOT / "potency_audit_matrix.csv"

GRAM_POSITIVES = {
    "clostridioides_difficile",
    "enterococcus_faecalis",
    "enterococcus_faecium",
    "listeria_monocytogenes",
    "staphylococcus_aureus",
    "staphylococcus_epidermidis",
    "streptococcus_agalactiae",
    "streptococcus_pneumoniae",
    "streptococcus_pyogenes",
}

ROW_LEVEL_INTRINSIC_INACTIVITY = {
    ("clostridioides_difficile", "flucloxacillin"),
    ("enterococcus_faecalis", "metronidazole"),
    ("enterococcus_faecium", "metronidazole"),
    ("listeria_monocytogenes", "metronidazole"),
    ("staphylococcus_aureus", "fidaxomicin"),
    ("staphylococcus_aureus", "metronidazole"),
    ("staphylococcus_epidermidis", "fidaxomicin"),
    ("staphylococcus_epidermidis", "metronidazole"),
    ("streptococcus_agalactiae", "fidaxomicin"),
    ("streptococcus_agalactiae", "metronidazole"),
    ("streptococcus_pneumoniae", "fidaxomicin"),
    ("streptococcus_pneumoniae", "metronidazole"),
    ("streptococcus_pyogenes", "fidaxomicin"),
    ("streptococcus_pyogenes", "metronidazole"),
}

ROW_LEVEL_INTRINSIC_NOTES = {
    ("clostridioides_difficile", "flucloxacillin"): "Anti-staphylococcal penicillin activity should not be carried into C. difficile; the current near-zero value looks like an explicit inactivity placeholder.",
    ("enterococcus_faecalis", "metronidazole"): "Metronidazole should not be interpreted as baseline anti-enterococcal activity; low residual values should be treated as explicit inactivity placeholders.",
    ("enterococcus_faecium", "metronidazole"): "Metronidazole should not be interpreted as baseline anti-enterococcal activity; low residual values should be treated as explicit inactivity placeholders.",
    ("listeria_monocytogenes", "metronidazole"): "Metronidazole should not be interpreted as meaningful baseline anti-Listeria activity; low residual values should be treated as explicit inactivity placeholders.",
    ("staphylococcus_aureus", "fidaxomicin"): "Fidaxomicin is not a meaningful systemic antistaphylococcal agent in this matrix; low residual values should be treated as explicit inactivity placeholders.",
    ("staphylococcus_aureus", "metronidazole"): "Metronidazole should not be interpreted as baseline anti-staphylococcal activity; low residual values should be treated as explicit inactivity placeholders.",
    ("staphylococcus_epidermidis", "fidaxomicin"): "Fidaxomicin is not a meaningful systemic anti-staphylococcal agent in this matrix; low residual values should be treated as explicit inactivity placeholders.",
    ("staphylococcus_epidermidis", "metronidazole"): "Metronidazole should not be interpreted as baseline anti-staphylococcal activity; low residual values should be treated as explicit inactivity placeholders.",
    ("streptococcus_agalactiae", "fidaxomicin"): "Fidaxomicin is not a meaningful non-C. difficile streptococcal agent in this matrix; low residual values should be treated as explicit inactivity placeholders.",
    ("streptococcus_agalactiae", "metronidazole"): "Metronidazole should not be interpreted as baseline anti-streptococcal activity; low residual values should be treated as explicit inactivity placeholders.",
    ("streptococcus_pneumoniae", "fidaxomicin"): "Fidaxomicin is not a meaningful non-C. difficile streptococcal agent in this matrix; low residual values should be treated as explicit inactivity placeholders.",
    ("streptococcus_pneumoniae", "metronidazole"): "Metronidazole should not be interpreted as baseline anti-streptococcal activity; low residual values should be treated as explicit inactivity placeholders.",
    ("streptococcus_pyogenes", "fidaxomicin"): "Fidaxomicin is not a meaningful non-C. difficile streptococcal agent in this matrix; low residual values should be treated as explicit inactivity placeholders.",
    ("streptococcus_pyogenes", "metronidazole"): "Metronidazole should not be interpreted as baseline anti-streptococcal activity; low residual values should be treated as explicit inactivity placeholders.",
}

DRUG_SPECIFIC_REVIEW_NOTES = {
    "vancomycin": "Glycopeptide activity across Gram positives is directionally plausible in many rows, but the current map should still be checked species by species rather than accepted wholesale.",
    "teicoplanin": "Glycopeptide activity across Gram positives is directionally plausible in many rows, but the current map should still be checked species by species rather than accepted wholesale.",
    "linezolid": "Oxazolidinone values are broadly high across the Gram-positive block except for suspicious low rows, so the full species map still needs explicit evidence review.",
    "daptomycin": "Uniformly low daptomycin values across this Gram-positive block are a high-priority review target and should not be accepted without explicit evidence checking.",
    "ampicillin": "Penicillin-class baseline activity differs materially across Gram-positive species and should be reviewed explicitly rather than left as inherited values.",
    "amoxicillin": "Penicillin-class baseline activity differs materially across Gram-positive species and should be reviewed explicitly rather than left as inherited values.",
    "penicillin_g": "Natural-penicillin activity should be reviewed species by species across the Gram-positive block rather than accepted from broad inheritance.",
    "flucloxacillin": "Anti-staphylococcal penicillin activity should be reviewed explicitly across Gram positives because meaningful activity is highly species-dependent.",
    "cefazolin": "Cephalosporin activity varies sharply across Gram-positive species and should be reviewed explicitly, especially where low default-like values remain.",
    "ceftriaxone": "Cephalosporin activity varies sharply across Gram-positive species and should be reviewed explicitly, especially where low default-like values remain.",
    "cefuroxime": "Cephalosporin activity in Gram positives should be reviewed species by species rather than inherited by class.",
    "cephalexin": "Cephalosporin activity in Gram positives should be reviewed species by species rather than inherited by class.",
    "cefepime": "Later-generation cephalosporin activity in Gram positives should be reviewed explicitly rather than inferred from class membership.",
    "ceftazidime": "Ceftazidime is not a routine Gram-positive anchor drug and any retained activity should be reviewed explicitly rather than accepted by inheritance.",
    "cefixime": "Oral cephalosporin activity in Gram positives should be reviewed explicitly rather than inferred from class membership.",
    "ceftaroline": "Ceftaroline is one of the more Gram-positive-relevant cephalosporins in the matrix, so its species-level values should be checked deliberately rather than inherited from broad cephalosporin logic.",
    "clindamycin": "Lincosamide activity varies substantially across Gram-positive species and should be reviewed explicitly rather than accepted as a broad class pattern.",
    "erythromycin": "Macrolide activity varies substantially across Gram-positive species and should be reviewed explicitly rather than accepted as a broad class pattern.",
    "azithromycin": "Macrolide activity varies substantially across Gram-positive species and should be reviewed explicitly rather than accepted as a broad class pattern.",
    "clarithromycin": "Macrolide activity varies substantially across Gram-positive species and should be reviewed explicitly rather than accepted as a broad class pattern.",
    "tetracycline": "Tetracycline-family activity should be reviewed species by species across Gram positives rather than accepted as broad class inheritance.",
    "doxycycline": "Tetracycline-family activity should be reviewed species by species across Gram positives rather than accepted as broad class inheritance.",
    "minocycline": "Tetracycline-family activity should be reviewed species by species across Gram positives rather than accepted as broad class inheritance.",
    "trim_sulf": "Trim-sulf activity differs materially across Gram-positive species and should be reviewed explicitly rather than accepted from inherited values.",
    "ciprofloxacin": "Fluoroquinolone activity varies across Gram-positive species and should be reviewed explicitly rather than accepted as a class-wide default.",
    "levofloxacin": "Fluoroquinolone activity varies across Gram-positive species and should be reviewed explicitly rather than accepted as a class-wide default.",
    "moxifloxacin": "Fluoroquinolone activity varies across Gram-positive species and should be reviewed explicitly rather than accepted as a class-wide default.",
    "nitrofurantoin": "Nitrofuran activity is highly species-dependent within Gram positives and should be reviewed explicitly rather than inferred from a few archetypes.",
    "fidaxomicin": "Fidaxomicin should be reviewed explicitly across Gram positives, especially because the current matrix appears to leave many non-C. difficile rows at low defaults while C. difficile itself remains suspiciously low.",
    "metronidazole": "Metronidazole should be reviewed explicitly across Gram positives because the C. difficile row is high while most other rows are low residual placeholders.",
    "colistin": "Polymyxin activity should not be assumed meaningful across Gram positives; any retained values should be reviewed explicitly rather than inherited from generic rules.",
    "aztreonam": "Monobactam activity should not be assumed meaningful across Gram positives; any retained values should be reviewed explicitly rather than inherited from generic Gram-negative logic.",
    "piperacillin": "Broad-spectrum beta-lactam values in Gram positives should be reviewed explicitly rather than inferred from Gram-negative-oriented class logic.",
    "ticarcillin": "Broad-spectrum beta-lactam values in Gram positives should be reviewed explicitly rather than inferred from Gram-negative-oriented class logic.",
    "amoxicillin_clavulanate": "BL/BLI values in Gram positives should be reviewed explicitly rather than inferred from Enterobacterales-oriented class logic.",
    "ampicillin_sulbactam": "BL/BLI values in Gram positives should be reviewed explicitly rather than inferred from Enterobacterales-oriented class logic.",
    "piperacillin_tazobactam": "BL/BLI values in Gram positives should be reviewed explicitly rather than inferred from Enterobacterales-oriented class logic.",
    "ticarcillin_clavulanate": "BL/BLI values in Gram positives should be reviewed explicitly rather than inferred from Enterobacterales-oriented class logic.",
    "gentamicin": "Aminoglycoside activity in Gram positives is context-dependent and should be reviewed explicitly rather than interpreted as straightforward baseline monotherapy potency.",
    "tobramycin": "Aminoglycoside activity in Gram positives is context-dependent and should be reviewed explicitly rather than interpreted as straightforward baseline monotherapy potency.",
    "amikacin": "Aminoglycoside activity in Gram positives is context-dependent and should be reviewed explicitly rather than interpreted as straightforward baseline monotherapy potency.",
    "fosfomycin": "Fosfomycin activity varies across Gram-positive species and should be reviewed explicitly rather than left as broad inherited values.",
    "tigecycline": "Tigecycline activity varies across Gram-positive species and should be reviewed explicitly rather than left as broad inherited values.",
    "rifampicin": "Rifampicin values across Gram positives should be reviewed explicitly because current entries may conflate baseline activity with rapid resistance selection concerns.",
    "chloramphenicol": "Chloramphenicol activity across Gram positives should be reviewed explicitly rather than accepted as a generic broad-spectrum carryover.",
    "quinu_dalfo": "Streptogramin activity across Gram positives should be reviewed explicitly rather than accepted as inherited class logic.",
    "dalbavancin": "Lipoglycopeptide activity across Gram positives should be reviewed explicitly, especially where low default-like rows remain.",
    "retapamulin": "Pleuromutilin activity should be reviewed explicitly across Gram positives rather than left as broad inherited values.",
    "tedizolid": "Oxazolidinone activity across Gram positives should be reviewed explicitly rather than accepted wholesale from class inheritance.",
    "fusidic_a": "Fusidic acid activity differs across Gram-positive species and should be reviewed explicitly rather than accepted from broad inheritance.",
}


def annotate_row(row: dict[str, str]) -> None:
    bacteria = row["bacteria"]
    drug = row["drug"]
    potency = float(row["potency_no_r"])

    if bacteria not in GRAM_POSITIVES:
        return

    row_key = (bacteria, drug)
    if row_key in ROW_LEVEL_INTRINSIC_INACTIVITY:
        row["review_status"] = "clear intrinsic inactivity"
        row["evidence_notes"] = ROW_LEVEL_INTRINSIC_NOTES[row_key]
        row["decision"] = "later normalize to explicit inactivity after full review"
        return

    if drug in DRUG_SPECIFIC_REVIEW_NOTES:
        row["review_status"] = "needs species-specific review"
        row["evidence_notes"] = DRUG_SPECIFIC_REVIEW_NOTES[drug]
        row["decision"] = "defer potency change until evidence review"
        return

    if potency <= 0.1:
        row["review_status"] = "fallback default; needs explicit review"
        row["evidence_notes"] = (
            "Gram-positive first-pass review: this row remains at a very low residual potency and should be treated as an unreviewed fallback until explicitly adjudicated."
        )
        row["decision"] = "defer potency change until evidence review"
        return

    row["review_status"] = "needs species-specific review"
    row["evidence_notes"] = (
        "Gram-positive first-pass review: retained activity should be reviewed explicitly by species rather than accepted from broad class inheritance."
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

    reviewed = sum(1 for row in rows if row["bacteria"] in GRAM_POSITIVES and row["review_status"])
    print(f"Annotated {reviewed} Gram-positive rows in {AUDIT_CSV.name}")


if __name__ == "__main__":
    main()