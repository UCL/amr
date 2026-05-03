from __future__ import annotations

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
AUDIT_CSV = REPO_ROOT / "potency_audit_matrix.csv"

ENTERIC_PATHOGENS = {
    "campylobacter_jejuni",
    "invasive_non-typhoidal_salmonella_spp.",
    "salmonella_enterica_serovar_paratyphi_a",
    "salmonella_enterica_serovar_typhi",
    "shigella_spp.",
    "vibrio_cholerae",
    "yersinia_enterocolitica",
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
    "clindamycin": "Lincosamides should not be interpreted as baseline activity for these enteric Gram-negative pathogens; any retained values are best treated as placeholders rather than reviewed spectrum.",
    "dalbavancin": "Lipoglycopeptides are Gram-positive-directed and should not carry baseline activity for this enteric pathogen block.",
    "daptomycin": "Daptomycin does not have meaningful baseline activity against these enteric Gram-negative pathogens; low residual values should be treated as cleanup targets, not real spectrum.",
    "fidaxomicin": "Fidaxomicin is not a meaningful enteric Gram-negative systemic agent in this matrix; low residual values should not be interpreted as real activity.",
    "flucloxacillin": "Anti-staphylococcal penicillin activity should not be carried into this enteric pathogen block; residual values are best treated as explicit inactivity placeholders.",
    "fusidic_a": "Fusidic acid should be treated as inactive for these enteric Gram-negative pathogens; weak residual values look inherited rather than reviewed.",
    "linezolid": "Oxazolidinones are Gram-positive-directed and should not retain baseline activity in this enteric pathogen block.",
    "quinu_dalfo": "Streptogramins should be treated as inactive for these enteric Gram-negative pathogens; residual values look inherited rather than curated.",
    "retapamulin": "Pleuromutilins are not meaningful agents for this enteric pathogen block; residual values should not be interpreted as real activity.",
    "tedizolid": "Oxazolidinones are Gram-positive-directed and should not retain baseline activity in this enteric pathogen block.",
    "teicoplanin": "Glycopeptides are not meaningful anti-enteric-Gram-negative agents; any residual values indicate incomplete explicit normalization.",
    "vancomycin": "Glycopeptides are not meaningful anti-enteric-Gram-negative agents; any residual values indicate incomplete explicit normalization.",
}

SPECIFIC_REVIEW_NOTES = {
    "amoxicillin": "Beta-lactam activity differs sharply across this enteric pathogen block, with plausible activity in several Enterobacterales-like organisms but very weak Campylobacter values, so the species map should be reviewed explicitly rather than inherited by class.",
    "ampicillin": "Beta-lactam activity differs sharply across this enteric pathogen block, with plausible activity in several Enterobacterales-like organisms but very weak Campylobacter values, so the species map should be reviewed explicitly rather than inherited by class.",
    "penicillin_g": "Natural-penicillin activity should be reviewed species by species across this enteric pathogen block rather than accepted from broad class inheritance.",
    "amoxicillin_clavulanate": "BL/BLI values in this enteric pathogen block should be reviewed explicitly because plausible activity in several species should not be conflated with the weak Campylobacter row.",
    "ampicillin_sulbactam": "BL/BLI values in this enteric pathogen block should be reviewed explicitly rather than accepted from broad class inheritance.",
    "piperacillin_tazobactam": "Broad-spectrum BL/BLI values in this enteric pathogen block should be reviewed explicitly rather than accepted from general Gram-negative inheritance.",
    "ticarcillin_clavulanate": "Broad-spectrum BL/BLI values in this enteric pathogen block should be reviewed explicitly rather than accepted from general Gram-negative inheritance.",
    "ceftriaxone": "Cephalosporin activity is species-dependent across this enteric pathogen block and should be reviewed explicitly, especially where Campylobacter remains very low while several other species are high.",
    "cefixime": "Uniformly high cefixime values across this enteric pathogen block, including Campylobacter, look suspicious and should be reviewed explicitly rather than accepted as class inheritance.",
    "cefiderocol": "Uniformly high cefiderocol values across this enteric pathogen block, including Campylobacter, look suspicious and should be confirmed as intentional evidence-based assignments rather than broad placeholders.",
    "cefuroxime": "Cephalosporin activity across this enteric pathogen block should be reviewed explicitly rather than accepted from class inheritance.",
    "cefazolin": "Early-generation cephalosporin activity across this enteric pathogen block should be reviewed explicitly rather than accepted from class inheritance.",
    "cephalexin": "Early-generation cephalosporin activity across this enteric pathogen block should be reviewed explicitly rather than accepted from class inheritance.",
    "cefepime": "Later-generation cephalosporin activity across this enteric pathogen block should be reviewed explicitly rather than accepted from class inheritance.",
    "ceftazidime": "Cephalosporin activity across this enteric pathogen block should be reviewed explicitly rather than accepted from class inheritance.",
    "ceftaroline": "Ceftaroline is not an obvious anchor drug for this enteric pathogen block and any retained activity should be reviewed explicitly.",
    "ceftolozane_tazobactam": "Ceftolozane-tazobactam values in this enteric pathogen block should be reviewed explicitly rather than accepted from broad Gram-negative inheritance.",
    "aztreonam": "Monobactam activity should be reviewed species by species across this enteric pathogen block rather than inferred from generic Gram-negative logic.",
    "aztreonam_avibactam": "Novel BL/BLI activity in this enteric pathogen block should be reviewed explicitly rather than accepted from class inheritance.",
    "ceftazidime_avibactam": "Novel BL/BLI activity in this enteric pathogen block should be reviewed explicitly rather than accepted from class inheritance.",
    "meropenem_vaborbactam": "Novel BL/BLI activity in this enteric pathogen block should be reviewed explicitly rather than accepted from class inheritance.",
    "ertapenem": "Carbapenem activity differs across this block and should be reviewed explicitly rather than accepted as broad class inheritance.",
    "imipenem_c": "Carbapenem activity differs across this block and should be reviewed explicitly rather than accepted as broad class inheritance.",
    "meropenem": "Carbapenem activity differs across this block and should be reviewed explicitly rather than accepted as broad class inheritance.",
    "azithromycin": "Macrolide activity differs materially across this enteric pathogen block, with stronger values in Shigella, Vibrio, and Campylobacter than in the salmonellae and Yersinia rows, so the species map should be reviewed explicitly.",
    "clarithromycin": "Macrolide activity differs materially across this enteric pathogen block and should be reviewed explicitly rather than accepted as a broad class pattern.",
    "erythromycin": "Macrolide activity differs materially across this enteric pathogen block and should be reviewed explicitly rather than accepted as a broad class pattern.",
    "ciprofloxacin": "Fluoroquinolone activity is broadly high across this enteric pathogen block and should be reviewed explicitly rather than accepted as a class-wide default.",
    "levofloxacin": "Fluoroquinolone activity is broadly high across this enteric pathogen block and should be reviewed explicitly rather than accepted as a class-wide default.",
    "moxifloxacin": "Fluoroquinolone activity is broadly high across this enteric pathogen block and should be reviewed explicitly rather than accepted as a class-wide default.",
    "trim_sulf": "Trim-sulf activity differs across this enteric pathogen block, especially between Campylobacter and several enteric species, so the species map should be reviewed explicitly rather than inherited.",
    "sulfanilamide": "Sulfonamide activity across this enteric pathogen block should be reviewed explicitly rather than accepted as broad inherited spectrum.",
    "tetracycline": "Tetracycline-family activity is species-dependent across this enteric pathogen block and should be reviewed explicitly.",
    "doxycycline": "Tetracycline-family activity is species-dependent across this enteric pathogen block and should be reviewed explicitly.",
    "minocycline": "Tetracycline-family activity is species-dependent across this enteric pathogen block and should be reviewed explicitly.",
    "chloramphenicol": "Chloramphenicol remains broadly high across this enteric pathogen block and should be reviewed explicitly rather than accepted as a historical broad-spectrum carryover.",
    "rifampicin": "Rifampicin values in this enteric pathogen block should be reviewed explicitly rather than accepted from broad inheritance.",
    "colistin": "Polymyxin activity across this enteric pathogen block should be reviewed explicitly rather than accepted from broad Gram-negative inheritance, especially given the relatively high non-Campylobacter rows.",
    "nitrofurantoin": "Nitrofurantoin values across this enteric pathogen block are low residuals and should be reviewed explicitly rather than left as inherited defaults.",
    "metronidazole": "Metronidazole values across this enteric pathogen block are near-zero residuals and should be reviewed explicitly rather than left as inherited defaults.",
    "fosfomycin": "Fosfomycin activity across this enteric pathogen block should be reviewed explicitly rather than accepted from broad inherited values.",
    "tigecycline": "Tigecycline activity across this enteric pathogen block should be reviewed explicitly rather than accepted from broad inherited values.",
    "gentamicin": "Aminoglycoside activity in this enteric pathogen block should be reviewed explicitly rather than interpreted as straightforward monotherapy baseline potency.",
    "tobramycin": "Aminoglycoside activity in this enteric pathogen block should be reviewed explicitly rather than interpreted as straightforward monotherapy baseline potency.",
    "amikacin": "Aminoglycoside activity in this enteric pathogen block should be reviewed explicitly rather than interpreted as straightforward monotherapy baseline potency.",
}


def annotate_row(row: dict[str, str]) -> None:
    bacteria = row["bacteria"]
    drug = row["drug"]
    potency = float(row["potency_no_r"])

    if bacteria not in ENTERIC_PATHOGENS:
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
            "Enteric pathogen first-pass review: this row remains at a very low residual potency and should be treated as an unreviewed fallback until explicitly adjudicated."
        )
        row["decision"] = "defer potency change until evidence review"
        return

    row["review_status"] = "needs species-specific review"
    row["evidence_notes"] = (
        "Enteric pathogen first-pass review: retained activity should be reviewed explicitly by species rather than accepted from broad class inheritance."
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

    reviewed = sum(1 for row in rows if row["bacteria"] in ENTERIC_PATHOGENS and row["review_status"])
    print(f"Annotated {reviewed} enteric pathogen rows in {AUDIT_CSV.name}")


if __name__ == "__main__":
    main()