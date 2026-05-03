from __future__ import annotations

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
AUDIT_CSV = REPO_ROOT / "potency_audit_matrix.csv"

NONFERMENTERS = {
    "acinetobacter_baumannii",
    "burkholderia_cepacia_complex",
    "pseudomonas_aeruginosa",
    "stenotrophomonas_maltophilia",
}

CLEAR_INTRINSIC_INACTIVITY = {
    "azithromycin",
    "clarithromycin",
    "clindamycin",
    "dalbavancin",
    "daptomycin",
    "fidaxomicin",
    "flucloxacillin",
    "fusidic_a",
    "linezolid",
    "metronidazole",
    "nitrofurantoin",
    "penicillin_g",
    "quinu_dalfo",
    "retapamulin",
    "tedizolid",
    "teicoplanin",
    "vancomycin",
}

INTRINSIC_NOTES = {
    "azithromycin": "Macrolides should not be interpreted as baseline anti-non-fermenter activity in this matrix; residual values look like inherited placeholders rather than intentional potency assignments.",
    "clarithromycin": "Macrolides should not be interpreted as baseline anti-non-fermenter activity in this matrix; residual values look like inherited placeholders rather than intentional potency assignments.",
    "clindamycin": "Lincosamides should be treated as inactive against these non-fermenters; residual low values are placeholder-like rather than curated.",
    "dalbavancin": "Lipoglycopeptides are Gram-positive-directed and should not carry baseline activity for these non-fermenters.",
    "daptomycin": "Daptomycin does not have meaningful baseline activity against these non-fermenters; any weak entries should be treated as cleanup targets, not real spectrum.",
    "fidaxomicin": "Fidaxomicin is not a meaningful anti-non-fermenter agent; low residual values should not be interpreted as real activity.",
    "flucloxacillin": "Anti-staphylococcal penicillin activity should not be carried into this non-fermenter block; low residual values are best treated as explicit inactivity placeholders.",
    "fusidic_a": "Fusidic acid should be treated as inactive for these non-fermenters; weak residual values look inherited rather than reviewed.",
    "linezolid": "Oxazolidinones are Gram-positive-directed and should not retain baseline non-fermenter activity.",
    "metronidazole": "Metronidazole is not a meaningful anti-non-fermenter agent; residual low values should eventually be normalized to explicit inactivity.",
    "nitrofurantoin": "Nitrofurantoin is a urinary drug with no meaningful baseline role against these non-fermenters; low values should be treated as explicit inactivity placeholders.",
    "penicillin_g": "Penicillin G should not retain baseline activity in this non-fermenter block; weak residual values are fallback-like rather than evidence-based.",
    "quinu_dalfo": "Streptogramins should be treated as inactive for these non-fermenters; low residual values look inherited rather than curated.",
    "retapamulin": "Pleuromutilins are not meaningful anti-non-fermenter agents in this matrix; residual values should not be interpreted as real activity.",
    "tedizolid": "Oxazolidinones are Gram-positive-directed and should not retain baseline non-fermenter activity.",
    "teicoplanin": "Glycopeptides are not anti-non-fermenter agents; any residual values indicate incomplete explicit normalization.",
    "vancomycin": "Glycopeptides are not anti-non-fermenter agents; any residual values indicate incomplete explicit normalization.",
}

SPECIFIC_REVIEW_NOTES = {
    "colistin": "Polymyxin activity is strongly species-dependent across non-fermenters; the current map needs explicit review rather than class-wide acceptance.",
    "tobramycin": "Aminoglycoside baseline activity varies substantially across non-fermenters and should be reviewed species by species.",
    "amikacin": "Aminoglycoside baseline activity varies substantially across non-fermenters and should be reviewed species by species.",
    "gentamicin": "Aminoglycoside baseline activity varies substantially across non-fermenters and should be reviewed species by species.",
    "ciprofloxacin": "Fluoroquinolone activity is clearly species-shaped in this block and should be reviewed explicitly rather than left as broad class inheritance.",
    "levofloxacin": "Fluoroquinolone activity is clearly species-shaped in this block and should be reviewed explicitly rather than left as broad class inheritance.",
    "moxifloxacin": "Moxifloxacin is not a standard anti-non-fermenter workhorse, so any retained activity should be checked species by species.",
    "ofloxacin": "Fluoroquinolone activity is clearly species-shaped in this block and should be reviewed explicitly rather than left as broad class inheritance.",
    "trim_sulf": "Trim-sulf activity is highly species-dependent across non-fermenters and should be reviewed explicitly, especially for Stenotrophomonas and Burkholderia.",
    "sulfanilamide": "Sulfonamide activity in this block looks placeholder-like in places and should be reviewed species by species rather than accepted by inheritance.",
    "chloramphenicol": "Chloramphenicol values across non-fermenters should be reviewed explicitly; broad retained activity here may reflect inheritance rather than an evidence-based map.",
    "rifampicin": "Rifampicin entries in non-fermenters should be treated cautiously and reviewed explicitly rather than inferred from class-level patterns.",
    "fosfomycin": "Fosfomycin activity in non-fermenters is species-dependent and often limited; the current map should be reviewed explicitly.",
    "tigecycline": "Tigecycline activity should be reviewed explicitly across non-fermenters; very low retained values may reflect defaults while true species differences are being missed.",
    "minocycline": "Minocycline activity is species-dependent across non-fermenters and deserves explicit review, especially for Acinetobacter and Stenotrophomonas.",
    "doxycycline": "Tetracycline-family activity is species-dependent across non-fermenters and should be reviewed explicitly.",
    "tetracycline": "Tetracycline-family activity is species-dependent across non-fermenters and should be reviewed explicitly.",
    "piperacillin": "Anti-pseudomonal beta-lactam activity should be reviewed species by species across this non-fermenter block rather than accepted as a generic class pattern.",
    "ticarcillin": "Carboxypenicillin activity should be reviewed species by species across this non-fermenter block rather than accepted as a generic class pattern.",
    "amoxicillin": "Weak aminopenicillin activity in this block is unlikely to be broadly meaningful and should be reviewed explicitly rather than left as inherited values.",
    "ampicillin": "Weak aminopenicillin activity in this block is unlikely to be broadly meaningful and should be reviewed explicitly rather than left as inherited values.",
    "amoxicillin_clavulanate": "BL/BLI activity is strongly species-dependent in non-fermenters and should be reviewed explicitly rather than inferred from Enterobacterales logic.",
    "ampicillin_sulbactam": "BL/BLI activity is strongly species-dependent in non-fermenters and should be reviewed explicitly rather than inferred from Enterobacterales logic.",
    "piperacillin_tazobactam": "Anti-pseudomonal BL/BLI activity varies by species in this block and should be reviewed explicitly rather than accepted as a broad class assignment.",
    "ticarcillin_clavulanate": "BL/BLI activity is strongly species-dependent in non-fermenters and should be reviewed explicitly rather than inferred from Enterobacterales logic.",
    "cefazolin": "Early-generation cephalosporin activity in non-fermenters is unlikely to be broadly meaningful; any retained values should be reviewed explicitly.",
    "cefuroxime": "Cephalosporin activity is species-dependent in this block and should be reviewed explicitly rather than inherited by class.",
    "cephalexin": "Early-generation cephalosporin activity in non-fermenters is unlikely to be broadly meaningful; any retained values should be reviewed explicitly.",
    "ceftriaxone": "Cephalosporin activity is species-dependent in this block and should be reviewed explicitly rather than inherited by class.",
    "ceftazidime": "Ceftazidime is a key differentiator across non-fermenters and should be reviewed species by species rather than accepted as a broad default.",
    "cefepime": "Cefepime is a key differentiator across non-fermenters and should be reviewed species by species rather than accepted as a broad default.",
    "ceftaroline": "Ceftaroline is not a standard non-fermenter drug and any retained activity should be reviewed explicitly rather than inferred from cephalosporin class membership.",
    "cefixime": "Oral cephalosporin activity in non-fermenters is unlikely to be broadly meaningful; any retained values should be reviewed explicitly.",
    "ceftolozane_tazobactam": "This is a high-value review target in non-fermenters because the current low values for key species may reflect unreviewed defaults rather than intentional spectrum choices.",
    "cefiderocol": "This is a high-value review target in non-fermenters because the current low values for key species may reflect unreviewed defaults rather than intentional spectrum choices.",
    "aztreonam": "Monobactam activity differs materially across non-fermenters and should be reviewed explicitly rather than inferred from generic Gram-negative logic.",
    "aztreonam_avibactam": "Novel BL/BLI activity in non-fermenters should be reviewed species by species; current values should not be treated as settled without evidence review.",
    "ceftazidime_avibactam": "Novel BL/BLI activity in non-fermenters should be reviewed species by species; current values should not be treated as settled without evidence review.",
    "meropenem_vaborbactam": "Novel BL/BLI activity in non-fermenters should be reviewed species by species; current values should not be treated as settled without evidence review.",
    "ertapenem": "Carbapenem activity is species-dependent in this block and should be reviewed explicitly; ertapenem in particular should not be inherited casually across non-fermenters.",
    "imipenem_c": "Carbapenem activity is species-dependent in this block and should be reviewed explicitly rather than accepted as a broad class pattern.",
    "meropenem": "Carbapenem activity is species-dependent in this block and should be reviewed explicitly rather than accepted as a broad class pattern.",
}


def annotate_row(row: dict[str, str]) -> None:
    bacteria = row["bacteria"]
    drug = row["drug"]
    potency = float(row["potency_no_r"])

    if bacteria not in NONFERMENTERS:
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
            "Non-fermenter first-pass review: this row remains at a very low residual potency and should be treated as an unreviewed fallback until explicitly adjudicated."
        )
        row["decision"] = "defer potency change until evidence review"
        return

    row["review_status"] = "needs species-specific review"
    row["evidence_notes"] = (
        "Non-fermenter first-pass review: retained activity should be reviewed explicitly by species rather than accepted from broad Gram-negative inheritance."
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
        1 for row in rows if row["bacteria"] in NONFERMENTERS and row["review_status"]
    )
    print(f"Annotated {reviewed} non-fermenter rows in {AUDIT_CSV.name}")


if __name__ == "__main__":
    main()