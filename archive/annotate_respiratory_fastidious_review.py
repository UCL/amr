from __future__ import annotations

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
AUDIT_CSV = REPO_ROOT / "potency_audit_matrix.csv"

RESPIRATORY_FASTIDIOUS = {
    "bordetella_pertussis",
    "haemophilus_influenzae",
    "legionella_pneumophila",
    "moraxella_catarrhalis",
    "neisseria_meningitidis",
}

CLEAR_INTRINSIC_INACTIVITY = {
    "clindamycin",
    "dalbavancin",
    "daptomycin",
    "fidaxomicin",
    "flucloxacillin",
    "fusidic_a",
    "linezolid",
    "quinu_dalfo",
    "retapamulin",
    "tedizolid",
    "teicoplanin",
    "vancomycin",
}

INTRINSIC_NOTES = {
    "clindamycin": "Lincosamides should not be interpreted as baseline activity for this respiratory/fastidious Gram-negative block; any retained values are best treated as placeholders rather than reviewed spectrum.",
    "dalbavancin": "Lipoglycopeptides are Gram-positive-directed and should not carry baseline activity for these fastidious Gram-negative organisms.",
    "daptomycin": "Daptomycin does not have meaningful baseline activity against these respiratory/fastidious Gram-negative organisms; low residual values should be treated as cleanup targets, not real spectrum.",
    "fidaxomicin": "Fidaxomicin is not a meaningful respiratory/fastidious Gram-negative agent; low residual values should not be interpreted as real activity.",
    "flucloxacillin": "Anti-staphylococcal penicillin activity should not be carried into this respiratory/fastidious Gram-negative block; low residual values are best treated as explicit inactivity placeholders.",
    "fusidic_a": "Fusidic acid should be treated as inactive for these respiratory/fastidious Gram-negative organisms; weak residual values look inherited rather than reviewed.",
    "linezolid": "Oxazolidinones are Gram-positive-directed and should not retain baseline activity in this respiratory/fastidious Gram-negative block.",
    "quinu_dalfo": "Streptogramins should be treated as inactive for these respiratory/fastidious Gram-negative organisms; residual values look inherited rather than curated.",
    "retapamulin": "Pleuromutilins are not meaningful agents for this respiratory/fastidious Gram-negative block; residual values should not be interpreted as real activity.",
    "tedizolid": "Oxazolidinones are Gram-positive-directed and should not retain baseline activity in this respiratory/fastidious Gram-negative block.",
    "teicoplanin": "Glycopeptides are not meaningful anti-fastidious-Gram-negative agents; any residual values indicate incomplete explicit normalization.",
    "vancomycin": "Glycopeptides are not meaningful anti-fastidious-Gram-negative agents; any residual values indicate incomplete explicit normalization.",
}

SPECIFIC_REVIEW_NOTES = {
    "amoxicillin": "Penicillin-class activity differs sharply across this respiratory/fastidious Gram-negative block, with plausible activity for some species and very low residual values for others; the species map should be reviewed explicitly rather than inherited by class.",
    "ampicillin": "Penicillin-class activity differs sharply across this respiratory/fastidious Gram-negative block, with plausible activity for some species and very low residual values for others; the species map should be reviewed explicitly rather than inherited by class.",
    "penicillin_g": "Natural-penicillin activity should be reviewed species by species across this respiratory/fastidious Gram-negative block rather than accepted from broad class inheritance.",
    "amoxicillin_clavulanate": "BL/BLI values in this block should be reviewed explicitly because plausible beta-lactam activity in Haemophilus, Moraxella, and meningococcus should not be conflated with the weak residual rows in Bordetella and Legionella.",
    "ampicillin_sulbactam": "BL/BLI values in this block should be reviewed explicitly because plausible beta-lactam activity in Haemophilus, Moraxella, and meningococcus should not be conflated with the weak residual rows in Bordetella and Legionella.",
    "piperacillin_tazobactam": "Broad-spectrum BL/BLI values in this block should be reviewed explicitly rather than inherited from general Gram-negative logic.",
    "ticarcillin_clavulanate": "Broad-spectrum BL/BLI values in this block should be reviewed explicitly rather than inherited from general Gram-negative logic.",
    "cefuroxime": "Cephalosporin activity is species-dependent across this block and should be reviewed explicitly, especially where Legionella and Bordetella retain very low residual values.",
    "ceftriaxone": "Cephalosporin activity is species-dependent across this block and should be reviewed explicitly, especially where Legionella and Bordetella retain very low residual values.",
    "cefixime": "Uniformly high cefixime values in several fastidious organisms, including Bordetella and Legionella, are suspicious and should be reviewed explicitly rather than accepted as class inheritance.",
    "cefazolin": "Early-generation cephalosporin activity in this block should be reviewed explicitly rather than inherited from broader class patterns.",
    "cephalexin": "Early-generation cephalosporin activity in this block should be reviewed explicitly rather than inherited from broader class patterns.",
    "cefepime": "Later-generation cephalosporin activity in this block should be reviewed explicitly rather than accepted as broad class inheritance.",
    "ceftazidime": "Cephalosporin activity in this block should be reviewed explicitly rather than accepted as broad class inheritance.",
    "ceftaroline": "Ceftaroline is not an obvious anchor drug for this fastidious Gram-negative block and any retained activity should be reviewed explicitly.",
    "cefiderocol": "Uniformly high cefiderocol values across this respiratory/fastidious block are suspicious and should be confirmed as intentional evidence-based assignments rather than broad placeholders.",
    "ceftolozane_tazobactam": "Ceftolozane-tazobactam values in this block should be reviewed explicitly rather than inherited from broader Gram-negative rules.",
    "aztreonam": "Monobactam activity should be reviewed species by species across this fastidious block rather than inferred from generic Gram-negative logic.",
    "aztreonam_avibactam": "Novel BL/BLI activity in this fastidious block should be reviewed explicitly rather than accepted from class inheritance.",
    "ceftazidime_avibactam": "Novel BL/BLI activity in this fastidious block should be reviewed explicitly rather than accepted from class inheritance.",
    "meropenem_vaborbactam": "Novel BL/BLI activity in this fastidious block should be reviewed explicitly rather than accepted from class inheritance.",
    "ertapenem": "Carbapenem activity differs across this block and should be reviewed explicitly rather than accepted as broad class inheritance.",
    "imipenem_c": "Carbapenem activity differs across this block and should be reviewed explicitly rather than accepted as broad class inheritance.",
    "meropenem": "Carbapenem activity differs across this block and should be reviewed explicitly rather than accepted as broad class inheritance.",
    "azithromycin": "Macrolides are central to some species in this block and peripheral or weaker in others, so the species map should be reviewed explicitly rather than accepted as a broad class pattern.",
    "clarithromycin": "Macrolides are central to some species in this block and peripheral or weaker in others, so the species map should be reviewed explicitly rather than accepted as a broad class pattern.",
    "erythromycin": "Macrolides are central to some species in this block and peripheral or weaker in others, so the species map should be reviewed explicitly rather than accepted as a broad class pattern.",
    "ciprofloxacin": "Fluoroquinolone activity is broadly high across this block and should be reviewed explicitly rather than accepted as a class-wide default.",
    "levofloxacin": "Fluoroquinolone activity is broadly high across this block and should be reviewed explicitly rather than accepted as a class-wide default.",
    "moxifloxacin": "Fluoroquinolone activity is broadly high across this block and should be reviewed explicitly rather than accepted as a class-wide default.",
    "trim_sulf": "Trim-sulf activity differs materially across this block, with especially weak residual values in Legionella; the species map should be reviewed explicitly rather than inherited.",
    "sulfanilamide": "Sulfonamide activity across this fastidious block should be reviewed explicitly rather than accepted as broad inherited spectrum.",
    "tetracycline": "Tetracycline-family activity is species-dependent across this fastidious block and should be reviewed explicitly.",
    "doxycycline": "Tetracycline-family activity is species-dependent across this fastidious block and should be reviewed explicitly.",
    "minocycline": "Tetracycline-family activity is species-dependent across this fastidious block and should be reviewed explicitly.",
    "rifampicin": "Rifampicin values vary across this block and should be reviewed explicitly rather than accepted from broad inheritance, especially given its distinct role in meningococcal and Legionella contexts.",
    "colistin": "Very low colistin values across this respiratory/fastidious block should be reviewed explicitly rather than left as residual placeholders.",
    "nitrofurantoin": "Nitrofurantoin values across this respiratory/fastidious block are low residuals and should be reviewed explicitly rather than inherited from generic defaults.",
    "metronidazole": "Metronidazole values across this respiratory/fastidious block are very low residuals and should be reviewed explicitly rather than inherited from generic defaults.",
}


def annotate_row(row: dict[str, str]) -> None:
    bacteria = row["bacteria"]
    drug = row["drug"]
    potency = float(row["potency_no_r"])

    if bacteria not in RESPIRATORY_FASTIDIOUS:
        return

    if drug in CLEAR_INTRINSIC_INACTIVITY:
        row["review_status"] = "clear intrinsic inactivity"
        row["evidence_notes"] = INTRINSIC_NOTES[drug]
        row["decision"] = "later normalize to explicit inactivity after full review"
        return

    if drug in SPECIFIC_REVIEW_NOTES:
        row["review_status"] = "needs species-specific review"
        row["evidence_notes"] = SPECIFIC_REVIEW_NOTES[drug]
        row["decision"] = "defer potency change until evidence review"
        return

    if potency <= 0.1:
        row["review_status"] = "fallback default; needs explicit review"
        row["evidence_notes"] = (
            "Respiratory/fastidious Gram-negative first-pass review: this row remains at a very low residual potency and should be treated as an unreviewed fallback until explicitly adjudicated."
        )
        row["decision"] = "defer potency change until evidence review"
        return

    row["review_status"] = "needs species-specific review"
    row["evidence_notes"] = (
        "Respiratory/fastidious Gram-negative first-pass review: retained activity should be reviewed explicitly by species rather than accepted from broad class inheritance."
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
        1 for row in rows if row["bacteria"] in RESPIRATORY_FASTIDIOUS and row["review_status"]
    )
    print(f"Annotated {reviewed} respiratory/fastidious Gram-negative rows in {AUDIT_CSV.name}")


if __name__ == "__main__":
    main()