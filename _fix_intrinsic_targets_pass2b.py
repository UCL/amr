"""
Pass 2B: Dot resistance targets in resistance_prevalence_values.csv where the
drug has potency < 0.15 (never selected in simulation) AND the resistance
figure represents intrinsic/not-applicable resistance rather than acquired.
These are deferred from pass 1 because they required clinical judgment on
whether potency should be raised (sub-case A, handled in pass 2A) vs target
should be blanked (sub-case B, handled here).
"""
import csv

PATH = "data/resistance_prevalence_values.csv"

PASS2B = {
    # Acinetobacter: aztreonam inactive (monobactam poor outer-membrane penetration
    # in Acinetobacter); tetracycline not used clinically for Acinetobacter.
    "Acinetobacter baumannii": {
        "drugs": {"aztreonam_avibactam", "tetracycline"},
        "note": (
            "aztreonam_avibactam/tetracycline dotted (pass 2): "
            "aztreonam inactive vs Acinetobacter (outer membrane); "
            "tetracycline not used for Acinetobacter treatment (both potency 0.10)."
        ),
    },
    # B. fragilis: penicillin G intrinsically inactive (constitutive BL), potency 0.10.
    "Bacteroides fragilis": {
        "drugs": {"penicillin_g"},
        "note": "penicillin_g dotted (pass 2): constitutive beta-lactamase (potency 0.10).",
    },
    # Citrobacter: intrinsic chromosomal AmpC makes aminopenicillins and 1st-gen
    # cephalosporins inactive; ceftaroline also poorly active (AmpC hydrolysis).
    "Citrobacter spp.": {
        "drugs": {"amoxicillin", "ampicillin", "cefazolin", "cephalexin", "ceftaroline"},
        "note": (
            "amoxicillin/ampicillin/cefazolin/cephalexin/ceftaroline dotted (pass 2): "
            "intrinsic AmpC beta-lactamase resistance; not used clinically (all potency 0.10)."
        ),
    },
    # C. diff: macrolides/clindamycin/daptomycin/erythromycin/tigecycline are tested
    # in surveillance but not used for C. diff treatment in simulation (potency 0.10).
    "Clostridioides difficile": {
        "drugs": {"azithromycin", "clarithromycin", "clindamycin", "daptomycin", "erythromycin", "tigecycline"},
        "note": (
            "6 drugs dotted (pass 2): not used for C. difficile treatment "
            "(azithromycin, clarithromycin, clindamycin, daptomycin, erythromycin, tigecycline; "
            "all potency 0.10 - surveillance data not applicable to treatment resistance)."
        ),
    },
    # Enterobacter: intrinsic AmpC (same as Citrobacter); ceftaroline AmpC-hydrolysed.
    "Enterobacter spp.": {
        "drugs": {"amoxicillin", "ampicillin", "cefazolin", "cephalexin", "ceftaroline"},
        "note": (
            "amoxicillin/ampicillin/cefazolin/cephalexin/ceftaroline dotted (pass 2): "
            "intrinsic AmpC; not used clinically (all potency 0.10)."
        ),
    },
    "Enterobacter cloacae": {
        "drugs": {"ceftaroline"},
        "note": "ceftaroline dotted (pass 2): intrinsic AmpC hydrolysis in E. cloacae (potency 0.10).",
    },
    # Enterococcus: fidaxomicin is a C. diff-only drug.
    "Enterococcus faecalis": {
        "drugs": {"fidaxomicin"},
        "note": "fidaxomicin dotted (pass 2): C. diff-specific drug; not used for Enterococcus (potency 0.10).",
    },
    "Enterococcus faecium": {
        "drugs": {"fidaxomicin"},
        "note": "fidaxomicin dotted (pass 2): C. diff-specific drug; not used for Enterococcus (potency 0.10).",
    },
    # H. pylori: tigecycline is not used for H. pylori eradication.
    "Helicobacter pylori": {
        "drugs": {"tigecycline"},
        "note": "tigecycline dotted (pass 2): not used for H. pylori treatment (potency 0.10).",
    },
    # Klebsiella: intrinsic SHV-1 chromosomal BL makes ampicillin/amoxicillin inactive.
    "Klebsiella pneumoniae": {
        "drugs": {"amoxicillin", "ampicillin"},
        "note": (
            "amoxicillin/ampicillin dotted (pass 2): "
            "intrinsic SHV-1 chromosomal beta-lactamase; not used for Klebsiella (potency 0.10)."
        ),
    },
    # Listeria: ticarcillin/clav not used; tigecycline not standard for Listeria.
    "Listeria monocytogenes": {
        "drugs": {"ticarcillin_clavulanate", "tigecycline"},
        "note": (
            "ticarcillin_clavulanate/tigecycline dotted (pass 2): "
            "not used for Listeria treatment (potency 0.10)."
        ),
    },
    # Moraxella: ceftaroline not standard for Moraxella (amox/azithromycin/FQ used).
    "Moraxella catarrhalis": {
        "drugs": {"ceftaroline"},
        "note": "ceftaroline dotted (pass 2): not standard for Moraxella treatment (potency 0.10).",
    },
    # Morganella: intrinsic AmpC (ceftaroline); intrinsic tetracycline family efflux
    # (Proteus-type tet pumps); tigecycline intrinsic efflux (no EUCAST breakpoint).
    "Morganella spp.": {
        "drugs": {"ceftaroline", "tetracycline", "doxycycline", "minocycline", "tigecycline"},
        "note": (
            "5 drugs dotted (pass 2): intrinsic AmpC (ceftaroline) and "
            "intrinsic tetracycline-family efflux including tigecycline "
            "(no EUCAST clinical breakpoints for Morganella; all potency 0.10)."
        ),
    },
    # Neisseria: ceftaroline/tigecycline not standard gonorrhea/meningitis treatment.
    "Neisseria gonorrhoeae": {
        "drugs": {"ceftaroline", "tigecycline"},
        "note": (
            "ceftaroline/tigecycline dotted (pass 2): "
            "not standard gonorrhea treatment (potency 0.10)."
        ),
    },
    "Neisseria meningitidis": {
        "drugs": {"ceftaroline", "tigecycline"},
        "note": (
            "ceftaroline/tigecycline dotted (pass 2): "
            "not standard meningitis treatment (potency 0.10)."
        ),
    },
    # Proteus: intrinsic tetracycline/tigecycline efflux (no EUCAST breakpoint for Proteus);
    # ceftaroline partially AmpC-mediated. Doxycycline intrinsic efflux.
    "Proteus spp.": {
        "drugs": {"ceftaroline", "tetracycline", "doxycycline", "tigecycline"},
        "note": (
            "4 drugs dotted (pass 2): intrinsic tetracycline-family efflux "
            "and ceftaroline poor activity (all potency 0.10; no EUCAST breakpoints)."
        ),
    },
    # Salmonella typhi/paratyphi/iNTS: ceftaroline not used for enteric fever/iNTS.
    "Salmonella enterica serovar typhi": {
        "drugs": {"ceftaroline"},
        "note": "ceftaroline dotted (pass 2): not used for typhoid treatment (potency 0.10).",
    },
    "Salmonella enterica serovar paratyphi a": {
        "drugs": {"ceftaroline"},
        "note": "ceftaroline dotted (pass 2): not used for paratyphoid treatment (potency 0.10).",
    },
    "Invasive non-typhoidal Salmonella spp.": {
        "drugs": {"ceftaroline"},
        "note": "ceftaroline dotted (pass 2): not standard for iNTS treatment (potency 0.10).",
    },
    # Serratia: intrinsic AmpC (ceftaroline); intrinsic tetracycline-family efflux;
    # tigecycline no EUCAST clinical breakpoint for Serratia.
    "Serratia spp.": {
        "drugs": {"ceftaroline", "tetracycline", "doxycycline", "tigecycline"},
        "note": (
            "4 drugs dotted (pass 2): intrinsic AmpC (ceftaroline) and "
            "intrinsic tetracycline-family efflux including tigecycline "
            "(no EUCAST breakpoints for Serratia; all potency 0.10)."
        ),
    },
    # Shigella: ceftaroline not used for shigellosis.
    "Shigella spp.": {
        "drugs": {"ceftaroline"},
        "note": "ceftaroline dotted (pass 2): not standard shigellosis treatment (potency 0.10).",
    },
    # Staph/Strep: fidaxomicin is a C. diff-only drug.
    "Staphylococcus aureus": {
        "drugs": {"fidaxomicin"},
        "note": "fidaxomicin dotted (pass 2): C. diff-specific drug (potency 0.10).",
    },
    "Streptococcus agalactiae": {
        "drugs": {"fidaxomicin"},
        "note": "fidaxomicin dotted (pass 2): C. diff-specific drug (potency 0.10).",
    },
    "Streptococcus pneumoniae": {
        "drugs": {"fidaxomicin"},
        "note": "fidaxomicin dotted (pass 2): C. diff-specific drug (potency 0.10).",
    },
    "Streptococcus pyogenes": {
        "drugs": {"fidaxomicin"},
        "note": "fidaxomicin dotted (pass 2): C. diff-specific drug (potency 0.10).",
    },
    # Bordetella: tigecycline not standard whooping cough treatment.
    "Bordetella pertussis": {
        "drugs": {"tigecycline"},
        "note": "tigecycline dotted (pass 2): not used for Bordetella treatment (potency 0.10).",
    },
    # Treponema: cefixime (0% target) and tigecycline not used for syphilis.
    "Treponema pallidum": {
        "drugs": {"cefixime", "tigecycline"},
        "note": "cefixime/tigecycline dotted (pass 2): not used for syphilis treatment (potency 0.10).",
    },
    # Vibrio: ceftaroline not used for cholera.
    "Vibrio cholerae": {
        "drugs": {"ceftaroline"},
        "note": "ceftaroline dotted (pass 2): not used for cholera treatment (potency 0.10).",
    },
    # Yersinia: ceftaroline not standard for Yersinia enterocolitica.
    "Yersinia enterocolitica": {
        "drugs": {"ceftaroline"},
        "note": "ceftaroline dotted (pass 2): not standard for Yersinia treatment (potency 0.10).",
    },
}

with open(PATH, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames[:]
    rows = list(reader)

notes_col = next((fn for fn in fieldnames if fn.lower() == "notes"), None)

total_changes = 0
rows_changed = 0
for row in rows:
    bacteria = row.get("Bacteria", "").strip()
    config = PASS2B.get(bacteria)
    if config is None:
        continue
    dotted = []
    for drug in config["drugs"]:
        if drug in row:
            old = row[drug].strip()
            if old not in (".", ""):
                row[drug] = "."
                dotted.append(drug)
                total_changes += 1
    if dotted:
        rows_changed += 1
        if notes_col:
            existing = row.get(notes_col, "").strip()
            new_note = config["note"]
            row[notes_col] = (existing + " " + new_note).strip() if existing else new_note
        print(f"  {bacteria}: dotted {len(dotted)}: {', '.join(sorted(dotted))}")

print(f"\nTotal: {total_changes} cells dotted across {rows_changed} bacteria rows.")

with open(PATH, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, quoting=csv.QUOTE_MINIMAL)
    writer.writeheader()
    writer.writerows(rows)
print("Written to", PATH)
