from __future__ import annotations

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
AUDIT_CSV = REPO_ROOT / "potency_audit_matrix.csv"

STI_ATYPICAL = {
    "chlamydia_trachomatis",
    "helicobacter_pylori",
    "mycoplasma_genitalium",
    "mycoplasma_pneumoniae",
    "neisseria_gonorrhoeae",
    "treponema_pallidum",
}

CLEAR_INTRINSIC_INACTIVITY = {
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
    "dalbavancin": "Lipoglycopeptides are not meaningful baseline agents for this STI/atypical block; any residual values should be treated as explicit inactivity placeholders.",
    "daptomycin": "Daptomycin does not have a meaningful role across this STI/atypical block; low residual values should be treated as cleanup targets, not real spectrum.",
    "fidaxomicin": "Fidaxomicin is not a meaningful STI/atypical agent in this matrix; low residual values should not be interpreted as real activity.",
    "flucloxacillin": "Anti-staphylococcal penicillin activity should not be carried into this STI/atypical block; residual values are best treated as explicit inactivity placeholders.",
    "fusidic_a": "Fusidic acid should not be treated as a baseline agent across this STI/atypical block; weak residual values look inherited rather than reviewed.",
    "linezolid": "Oxazolidinones are not meaningful baseline agents for this STI/atypical block; low residual values should be treated as explicit inactivity placeholders.",
    "quinu_dalfo": "Streptogramins should not be treated as baseline agents across this STI/atypical block; residual values look inherited rather than curated.",
    "retapamulin": "Pleuromutilins are not meaningful baseline agents for this STI/atypical block; residual values should not be interpreted as real activity.",
    "tedizolid": "Oxazolidinones are not meaningful baseline agents for this STI/atypical block; low residual values should be treated as explicit inactivity placeholders.",
    "teicoplanin": "Glycopeptides are not meaningful baseline agents for this STI/atypical block; any residual values indicate incomplete explicit normalization.",
    "vancomycin": "Glycopeptides are not meaningful baseline agents for this STI/atypical block; any residual values indicate incomplete explicit normalization.",
}

CELL_WALL_DRUGS = {
    "amoxicillin",
    "amoxicillin_clavulanate",
    "ampicillin",
    "ampicillin_sulbactam",
    "aztreonam",
    "aztreonam_avibactam",
    "cefazolin",
    "cefixime",
    "cefiderocol",
    "cefepime",
    "ceftaroline",
    "ceftazidime",
    "ceftazidime_avibactam",
    "ceftolozane_tazobactam",
    "ceftriaxone",
    "cefuroxime",
    "cephalexin",
    "ertapenem",
    "imipenem_c",
    "meropenem",
    "meropenem_vaborbactam",
    "penicillin_g",
    "piperacillin",
    "piperacillin_tazobactam",
    "ticarcillin",
    "ticarcillin_clavulanate",
}

CELL_WALL_FREE_SPECIES = {
    "chlamydia_trachomatis",
    "mycoplasma_genitalium",
    "mycoplasma_pneumoniae",
}

DRUG_SPECIFIC_REVIEW_NOTES = {
    "penicillin_g": "Natural-penicillin activity varies sharply across this STI/atypical block, from treponemal high activity to near-null atypical rows, so the species map should be reviewed explicitly rather than inherited by class.",
    "amoxicillin": "Beta-lactam activity is highly species-dependent in this STI/atypical block and should be reviewed explicitly rather than accepted from broad class inheritance.",
    "ampicillin": "Beta-lactam activity is highly species-dependent in this STI/atypical block and should be reviewed explicitly rather than accepted from broad class inheritance.",
    "amoxicillin_clavulanate": "BL/BLI values in this STI/atypical block should be reviewed explicitly rather than accepted from broad class inheritance.",
    "ampicillin_sulbactam": "BL/BLI values in this STI/atypical block should be reviewed explicitly rather than accepted from broad class inheritance.",
    "ceftriaxone": "Cephalosporin activity is species-dependent across this STI/atypical block and should be reviewed explicitly, especially where gonococcal and treponemal rows diverge sharply from atypical near-null rows.",
    "cefixime": "Current cefixime values are especially suspicious in this block, with a partial gonococcal value and near-null atypical rows; this should be reviewed explicitly rather than accepted as class inheritance.",
    "cefiderocol": "Current cefiderocol values in this STI/atypical block look suspiciously inherited, especially in gonococcus and the atypical organisms, and should be reviewed explicitly.",
    "azithromycin": "Macrolides are central to several organisms in this STI/atypical block, but the species map still needs explicit review rather than acceptance as a broad class pattern.",
    "clarithromycin": "Macrolides are central to several organisms in this STI/atypical block, but the species map still needs explicit review rather than acceptance as a broad class pattern.",
    "erythromycin": "Macrolides are central to several organisms in this STI/atypical block, but the species map still needs explicit review rather than acceptance as a broad class pattern.",
    "doxycycline": "Tetracycline-family activity is central to several organisms in this STI/atypical block and should be reviewed explicitly by species rather than accepted from broad class inheritance.",
    "minocycline": "Tetracycline-family activity is central to several organisms in this STI/atypical block and should be reviewed explicitly by species rather than accepted from broad class inheritance.",
    "tetracycline": "Tetracycline-family activity is central to several organisms in this STI/atypical block and should be reviewed explicitly by species rather than accepted from broad class inheritance.",
    "ciprofloxacin": "Fluoroquinolone activity differs materially across this STI/atypical block and should be reviewed explicitly rather than accepted as a class-wide default.",
    "levofloxacin": "Fluoroquinolone activity differs materially across this STI/atypical block and should be reviewed explicitly rather than accepted as a class-wide default.",
    "moxifloxacin": "Fluoroquinolone activity differs materially across this STI/atypical block and should be reviewed explicitly rather than accepted as a class-wide default.",
    "gentamicin": "Aminoglycoside activity in this STI/atypical block is uneven and context-dependent, so the species map should be reviewed explicitly rather than interpreted as straightforward monotherapy potency.",
    "trim_sulf": "Trim-sulf activity differs materially across this STI/atypical block and should be reviewed explicitly rather than accepted from inherited values.",
    "metronidazole": "Metronidazole values in this block are highly species-dependent, with Helicobacter standing apart from near-null rows elsewhere, so the species map should be reviewed explicitly.",
    "rifampicin": "Rifampicin values in this STI/atypical block should be reviewed explicitly rather than accepted from broad inheritance, especially given its distinct roles in atypical respiratory and Helicobacter contexts.",
    "colistin": "Colistin values in this STI/atypical block should be reviewed explicitly rather than left as residual placeholders or inherited Gram-negative defaults.",
    "nitrofurantoin": "Nitrofurantoin values in this STI/atypical block should be reviewed explicitly rather than left as residual placeholders or inherited defaults.",
    "cefuroxime": "Cephalosporin activity in this STI/atypical block should be reviewed explicitly rather than accepted from class inheritance.",
    "cefazolin": "Cephalosporin activity in this STI/atypical block should be reviewed explicitly rather than accepted from class inheritance.",
    "cephalexin": "Cephalosporin activity in this STI/atypical block should be reviewed explicitly rather than accepted from class inheritance.",
    "cefepime": "Cephalosporin activity in this STI/atypical block should be reviewed explicitly rather than accepted from class inheritance.",
    "ceftazidime": "Cephalosporin activity in this STI/atypical block should be reviewed explicitly rather than accepted from class inheritance.",
    "ceftaroline": "Ceftaroline is not an obvious anchor drug for this STI/atypical block and any retained activity should be reviewed explicitly.",
    "ceftolozane_tazobactam": "Ceftolozane-tazobactam values in this STI/atypical block should be reviewed explicitly rather than accepted from broad Gram-negative inheritance.",
    "aztreonam": "Monobactam activity should be reviewed species by species across this STI/atypical block rather than inferred from generic Gram-negative logic.",
    "aztreonam_avibactam": "Novel BL/BLI activity in this STI/atypical block should be reviewed explicitly rather than accepted from class inheritance.",
    "ceftazidime_avibactam": "Novel BL/BLI activity in this STI/atypical block should be reviewed explicitly rather than accepted from class inheritance.",
    "meropenem_vaborbactam": "Novel BL/BLI activity in this STI/atypical block should be reviewed explicitly rather than accepted from class inheritance.",
    "ertapenem": "Carbapenem activity in this STI/atypical block should be reviewed explicitly rather than accepted as broad class inheritance.",
    "imipenem_c": "Carbapenem activity in this STI/atypical block should be reviewed explicitly rather than accepted as broad class inheritance.",
    "meropenem": "Carbapenem activity in this STI/atypical block should be reviewed explicitly rather than accepted as broad class inheritance.",
}


def annotate_row(row: dict[str, str]) -> None:
    bacteria = row["bacteria"]
    drug = row["drug"]
    potency = float(row["potency_no_r"])

    if bacteria not in STI_ATYPICAL:
        return

    if drug in CLEAR_INTRINSIC_INACTIVITY:
        row["review_status"] = "clear intrinsic inactivity"
        row["evidence_notes"] = INTRINSIC_NOTES[drug]
        row["decision"] = "later normalize to explicit inactivity after full review"
        return

    if bacteria in CELL_WALL_FREE_SPECIES and drug in CELL_WALL_DRUGS:
        row["review_status"] = "clear intrinsic inactivity"
        row["evidence_notes"] = (
            "This organism lacks the usual cell-wall target for beta-lactams and related cell-wall-active agents, so the current low residual value should be treated as an explicit inactivity placeholder rather than real activity."
        )
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
            "STI/atypical first-pass review: this row remains at a very low residual potency and should be treated as an unreviewed fallback until explicitly adjudicated."
        )
        row["decision"] = "defer potency change until evidence review"
        return

    row["review_status"] = "needs species-specific review"
    row["evidence_notes"] = (
        "STI/atypical first-pass review: retained activity should be reviewed explicitly by species rather than accepted from broad class inheritance."
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

    reviewed = sum(1 for row in rows if row["bacteria"] in STI_ATYPICAL and row["review_status"])
    print(f"Annotated {reviewed} STI/atypical rows in {AUDIT_CSV.name}")


if __name__ == "__main__":
    main()