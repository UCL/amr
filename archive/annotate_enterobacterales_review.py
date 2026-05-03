from __future__ import annotations

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
AUDIT_CSV = REPO_ROOT / "potency_audit_matrix.csv"

ENTEROBACTERALES = {
    "citrobacter_spp.",
    "enterobacter_cloacae",
    "enterobacter_spp.",
    "escherichia_coli",
    "klebsiella_pneumoniae",
    "morganella_spp.",
    "p_stuartii",
    "proteus_spp.",
    "serratia_spp.",
}

CLEAR_INTRINSIC_INACTIVITY = {
    "azithromycin",
    "clarithromycin",
    "flucloxacillin",
    "erythromycin",
    "clindamycin",
    "dalbavancin",
    "daptomycin",
    "fidaxomicin",
    "fusidic_a",
    "linezolid",
    "metronidazole",
    "quinu_dalfo",
    "retapamulin",
    "tedizolid",
    "teicoplanin",
    "vancomycin",
}

INTRINSIC_NOTES = {
    "azithromycin": "Macrolides are not meaningful baseline therapy for Enterobacterales; current 0.05/0.1 residues look like fallback placeholders rather than intentional activity.",
    "clarithromycin": "Macrolides are not meaningful baseline therapy for Enterobacterales; current 0.05/0.1 residues look like fallback placeholders rather than intentional activity.",
    "flucloxacillin": "Anti-staphylococcal penicillin activity should not be carried into Enterobacterales; current near-zero residual values are best treated as explicit inactivity placeholders.",
    "erythromycin": "Macrolides are not meaningful baseline therapy for Enterobacterales; current 0.05/0.1 residues look like fallback placeholders rather than intentional activity.",
    "clindamycin": "Lincosamides should be treated as intrinsically inactive for Enterobacterales; mixed 0/0.05/0.1 values indicate incomplete cleanup of fallback entries.",
    "dalbavancin": "Lipoglycopeptides are Gram-positive-directed and should not carry residual Enterobacterales activity.",
    "daptomycin": "Daptomycin does not have meaningful Enterobacterales activity; uniform 0.1 values are especially suggestive of untouched defaults.",
    "fidaxomicin": "Fidaxomicin is not a meaningful Enterobacterales-active drug; uniform 0.1 values look like inherited defaults rather than evidence-based assignments.",
    "fusidic_a": "Fusidic acid should be treated as inactive for Enterobacterales; residual 0.05 values are weak-placeholder style entries.",
    "linezolid": "Oxazolidinones are Gram-positive-directed and should not retain weak default Enterobacterales activity.",
    "metronidazole": "Metronidazole is not active against Enterobacterales; residual 0.05 values should eventually be normalized to explicit inactivity.",
    "quinu_dalfo": "Streptogramins should be treated as inactive for Enterobacterales; current weak residual values look inherited rather than curated.",
    "retapamulin": "Pleuromutilins are not meaningful Enterobacterales agents in this matrix; residual 0.05 values should not be interpreted as real activity.",
    "tedizolid": "Oxazolidinones are Gram-positive-directed and should not retain weak default Enterobacterales activity.",
    "teicoplanin": "Glycopeptides are not Enterobacterales-active; mixed 0/0.05/0.1 values indicate incomplete explicit normalization.",
    "vancomycin": "Glycopeptides are not Enterobacterales-active; mixed 0/0.05/0.1 values indicate incomplete explicit normalization.",
}

ENTERO_REVIEW_BUCKETS = {
    "colistin": "review species-specific intrinsic resistance map",
    "nitrofurantoin": "review urinary-spectrum species map",
    "fosfomycin": "review species-specific activity",
    "tigecycline": "review species-specific activity",
    "tetracycline": "review tetracycline family by species",
    "doxycycline": "review tetracycline family by species",
    "minocycline": "review tetracycline family by species",
    "amoxicillin": "review beta-lactam baseline by species",
    "ampicillin": "review beta-lactam baseline by species",
    "piperacillin": "review beta-lactam baseline by species",
    "ticarcillin": "review beta-lactam baseline by species",
    "amoxicillin_clavulanate": "review BL/BLI baseline by species",
    "ampicillin_sulbactam": "review BL/BLI baseline by species",
    "piperacillin_tazobactam": "review BL/BLI baseline by species",
    "ticarcillin_clavulanate": "review BL/BLI baseline by species",
    "cefazolin": "review cephalosporin baseline by species",
    "cefuroxime": "review cephalosporin baseline by species",
    "cephalexin": "review cephalosporin baseline by species",
    "ceftriaxone": "review cephalosporin baseline by species",
    "ceftazidime": "review cephalosporin baseline by species",
    "cefepime": "review cephalosporin baseline by species",
    "ceftaroline": "review cephalosporin baseline by species",
    "cefixime": "review cephalosporin baseline by species",
    "ceftolozane_tazobactam": "review newer cephalosporin/BLI baseline by species",
    "cefiderocol": "review siderophore cephalosporin baseline by species",
    "gentamicin": "review aminoglycoside baseline by species",
    "tobramycin": "review aminoglycoside baseline by species",
    "amikacin": "review aminoglycoside baseline by species",
    "ciprofloxacin": "review fluoroquinolone baseline by species",
    "levofloxacin": "review fluoroquinolone baseline by species",
    "moxifloxacin": "review fluoroquinolone baseline by species",
    "ofloxacin": "review fluoroquinolone baseline by species",
    "trim_sulf": "review folate antagonist baseline by species",
    "sulfanilamide": "review folate antagonist baseline by species",
    "chloramphenicol": "review chloramphenicol baseline by species",
    "rifampicin": "review rifamycin baseline by species",
    "aztreonam": "review monobactam baseline by species",
    "aztreonam_avibactam": "review novel BL/BLI baseline by species",
    "ceftazidime_avibactam": "review novel BL/BLI baseline by species",
    "meropenem_vaborbactam": "review novel BL/BLI baseline by species",
    "ertapenem": "review carbapenem baseline by species",
    "imipenem_c": "review carbapenem baseline by species",
    "meropenem": "review carbapenem baseline by species",
}

SPECIFIC_REVIEW_NOTES = {
    "colistin": "Current Enterobacterales map likely under-represents intrinsic polymyxin resistance beyond p_stuartii; Proteus, Morganella, Serratia, and Providencia-pattern organisms need explicit adjudication.",
    "nitrofurantoin": "Current pattern gives substantial activity to several non-E. coli Enterobacterales; urinary-spectrum species differences need explicit review, especially Proteus/Morganella/Serratia/Providencia-pattern organisms.",
    "fosfomycin": "Uniform 0.1 across Enterobacterales is almost certainly a fallback state, not a reviewed organism-specific evidence map.",
    "tigecycline": "Uniform 0.1 across Enterobacterales is implausible for a drug retained specifically for Gram-negative MDR use and should be reviewed from scratch.",
    "tetracycline": "Current Enterobacterales tetracycline values remain broadly high across most species after removal of the explicit calibration block; family should be re-reviewed from microbiology only.",
    "doxycycline": "Current Enterobacterales tetracycline-family values remain broadly high across most species; needs species-specific evidence review.",
    "minocycline": "Current Enterobacterales tetracycline-family values remain broadly high across most species; needs species-specific evidence review.",
    "amoxicillin": "Aminopenicillin values split sharply by species, but the current map still mixes clear species logic with defaults/placeholders and needs evidence normalization.",
    "ampicillin": "Aminopenicillin values split sharply by species, but the current map still mixes clear species logic with defaults/placeholders and needs evidence normalization.",
    "piperacillin": "Anti-pseudomonal penicillin values are broadly high across Enterobacterales with a weaker p_stuartii value; should be checked against species-level intrinsic/acquired susceptibility expectations.",
    "ticarcillin": "Carboxypenicillin values are broadly high across Enterobacterales with a weaker p_stuartii value; should be checked against species-level intrinsic/acquired susceptibility expectations.",
    "amoxicillin_clavulanate": "BL/BLI values appear hand-shaped by species; review should distinguish plausible inhibitor rescue from overstated baseline activity in inducible AmpC groups.",
    "ampicillin_sulbactam": "BL/BLI values appear hand-shaped by species; review should distinguish plausible inhibitor rescue from overstated baseline activity in inducible AmpC groups.",
    "piperacillin_tazobactam": "Broadly very high values may be directionally plausible for some species, but the entire species map should be checked rather than accepted as calibrated truth.",
    "ticarcillin_clavulanate": "Broadly high BL/BLI values need review by species, particularly for organisms with known inducible AmpC or intrinsic resistance patterns.",
    "cefazolin": "First-generation cephalosporin pattern looks partially species-aware, but residual 0.1 defaults remain for several organisms and should be normalized deliberately.",
    "cefuroxime": "Second-generation cephalosporin values are heterogeneous and need explicit species-level evidence review.",
    "cephalexin": "First-generation cephalosporin pattern looks partially species-aware, but residual 0.1 defaults remain for several organisms and should be normalized deliberately.",
    "ceftriaxone": "Third-generation cephalosporin values vary widely by species and include notably lower entries for some AmpC-associated organisms; needs deliberate evidence review.",
    "ceftazidime": "Third-generation cephalosporin values are broadly high but still species-shaped; review should verify this is microbiology-first rather than inherited tuning.",
    "cefepime": "Fourth-generation cephalosporin values are uniformly high across the block; likely directionally plausible but still need explicit species-level confirmation.",
    "ceftaroline": "Anti-MRSA cephalosporin values look suspiciously generous for Enterobacterales and should be treated as a high-priority review target.",
    "cefixime": "Oral third-generation cephalosporin values are uniformly high across the block; this may be directionally plausible in susceptible isolates but should still be confirmed species by species.",
    "ceftolozane_tazobactam": "Uniform 0.8 across the Enterobacterales block looks more like a class-wide inheritance rule than a reviewed species map and should be checked explicitly.",
    "cefiderocol": "Uniform 0.8 across the Enterobacterales block is plausible in direction but should be confirmed as an intentional evidence-based assignment rather than a broad placeholder.",
    "gentamicin": "Aminoglycoside values are broadly high with a weaker p_stuartii entry; direction may be plausible, but species-level baseline activity still needs explicit review.",
    "tobramycin": "Aminoglycoside values are broadly high with a weaker p_stuartii entry; direction may be plausible, but species-level baseline activity still needs explicit review.",
    "amikacin": "Aminoglycoside values are broadly high and especially generous across the block; verify that this is intended species-level baseline activity rather than broad inheritance.",
    "ciprofloxacin": "Fluoroquinolone values are very high across most Enterobacterales with only a weaker p_stuartii entry; this needs an evidence-first species review rather than acceptance by default.",
    "levofloxacin": "Fluoroquinolone values are very high across most Enterobacterales with only a weaker p_stuartii entry; this needs an evidence-first species review rather than acceptance by default.",
    "moxifloxacin": "Moxifloxacin is given fairly broad Enterobacterales activity here; review should confirm that this is microbiologically intended rather than carried over from class inheritance.",
    "ofloxacin": "Fluoroquinolone values are broadly high across the block with a weaker p_stuartii entry; review should confirm species-level realism.",
    "trim_sulf": "Trim-sulf values are broadly high across the block; this may be directionally plausible for susceptible isolates, but species-level baseline activity should still be checked explicitly.",
    "sulfanilamide": "Uniform 0.5 across most Enterobacterales is suspicious for an obsolete agent and should be reviewed as a likely inherited or placeholder-style assignment.",
    "chloramphenicol": "Broadly high chloramphenicol values across Enterobacterales deserve review; the current map may be overly generous for a drug with variable modern susceptibility patterns.",
    "rifampicin": "Rifampicin values look surprisingly generous across Enterobacterales given its limited standalone role and rapid resistance selection; treat as a review target.",
    "aztreonam": "Monobactam values are broadly high across the block with a weaker p_stuartii entry; direction may be plausible but should still be confirmed species by species.",
    "aztreonam_avibactam": "Uniformly maximal values across the Enterobacterales block are directionally plausible for spectrum, but should still be confirmed as an intentional evidence-based assignment.",
    "ceftazidime_avibactam": "Very high novel BL/BLI values are directionally plausible, but the map should still be checked against intended spectrum and intrinsic non-target organisms.",
    "meropenem_vaborbactam": "Uniformly maximal values across the block are directionally plausible but still deserve explicit confirmation during review.",
    "ertapenem": "Carbapenem values are broadly very high and likely directionally plausible, but species-level distinctions should be verified explicitly.",
    "imipenem_c": "Carbapenem values are broadly very high and likely directionally plausible, but species-level distinctions should be verified explicitly.",
    "meropenem": "Carbapenem values are broadly very high and likely directionally plausible, but species-level distinctions should be verified explicitly.",
}

FALLBACK_REVIEW_NOTES = {
    "penicillin_g": "Penicillin G remains at default-like weak activity across Enterobacterales; this should be normalized deliberately rather than left as a soft fallback residue.",
    "furazolidone": "Furazolidone remains at default-like weak activity across Enterobacterales and should be reviewed explicitly rather than inherited from fallback values.",
}


def annotate_row(row: dict[str, str]) -> None:
    bacteria = row["bacteria"]
    drug = row["drug"]
    potency = float(row["potency_no_r"])

    if bacteria not in ENTEROBACTERALES:
        return

    if drug in CLEAR_INTRINSIC_INACTIVITY:
        row["review_status"] = "clear intrinsic inactivity"
        row["evidence_notes"] = INTRINSIC_NOTES[drug]
        row["decision"] = "later normalize to explicit inactivity after full review"
        return

    if drug in ENTERO_REVIEW_BUCKETS:
        row["review_status"] = "needs species-specific review"
        row["evidence_notes"] = SPECIFIC_REVIEW_NOTES[drug]
        row["decision"] = "defer potency change until evidence review"
        return

    if drug in FALLBACK_REVIEW_NOTES and potency <= 0.1 and not row["review_status"]:
        row["review_status"] = "fallback default; needs explicit review"
        row["evidence_notes"] = FALLBACK_REVIEW_NOTES.get(
            drug,
            "Enterobacterales first-pass review: row is still at model default 0.1 rather than an obviously adjudicated organism-drug value.",
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
        1
        for row in rows
        if row["bacteria"] in ENTEROBACTERALES and row["review_status"]
    )
    print(f"Annotated {reviewed} Enterobacterales rows in {AUDIT_CSV.name}")


if __name__ == "__main__":
    main()