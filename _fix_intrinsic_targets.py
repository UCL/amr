"""
Pass 1: Dot out resistance targets that represent constitutive/intrinsic resistance
for drugs that are never used in the simulation (potency < 0.15).
Scope: only pairs with unambiguous biological justification (potency <= 0.05,
or Stenotrophomonas where all listed drugs are clearly not used clinically).
Deferred to pass 2: tigecycline, ceftaroline, fidaxomicin at potency 0.10
(need potency review first).
"""
import csv

PATH = "data/resistance_prevalence_values.csv"

# bacteria display name (exact as in CSV) -> drugs to dot + row note to append
PASS1 = {
    # Steno: intrinsic L1/L2 metallobeta-lactamases + SmeABC/SmeDEF efflux pumps
    # render all 30 drugs inactive; none are used clinically for Steno treatment.
    "Stenotrophomonas maltophilia": {
        "drugs": {
            "amikacin", "amoxicillin", "amoxicillin_clavulanate", "ampicillin",
            "ampicillin_sulbactam", "aztreonam", "aztreonam_avibactam", "cefazolin",
            "cefepime", "cefixime", "ceftaroline", "ceftazidime_avibactam",
            "ceftriaxone", "cefuroxime", "cephalexin", "clindamycin", "colistin",
            "ertapenem", "gentamicin", "imipenem_c", "meropenem",
            "meropenem_vaborbactam", "nitrofurantoin", "penicillin_g", "piperacillin",
            "piperacillin_tazobactam", "ticarcillin", "ticarcillin_clavulanate",
            "tigecycline", "tobramycin",
        },
        "note": (
            "30 drugs dotted (pass 1): intrinsic L1/L2 beta-lactamases + efflux "
            "pumps - surveillance figures are intrinsic not acquired "
            "(all potency < 0.15; none selected for treatment in simulation)."
        ),
    },
    # B. fragilis: constitutive chromosomal cephalosporinase (CepA) + cfxA;
    # aztreonam/BLI combos also intrinsically inactive against strict anaerobes.
    "Bacteroides fragilis": {
        "drugs": {
            "amoxicillin", "aztreonam_avibactam", "cefazolin",
            "cephalexin", "piperacillin", "ticarcillin",
        },
        "note": (
            "6 drugs dotted (pass 1): constitutive CepA/cfxA beta-lactamase "
            "+ aztreonam inactive vs anaerobes - intrinsic not acquired "
            "(potency 0.000-0.050)."
        ),
    },
    # Aztreonam has zero intrinsic activity vs Gram-positive cocci (outer membrane absent
    # but PBP binding differs; monobactam spectrum is Gram-negative only).
    "Streptococcus pneumoniae": {
        "drugs": {"aztreonam", "aztreonam_avibactam"},
        "note": (
            "aztreonam/aztreonam_avibactam dotted (pass 1): "
            "no activity against Gram-positive cocci (potency 0.000/0.010)."
        ),
    },
    "Streptococcus pyogenes": {
        "drugs": {"aztreonam", "aztreonam_avibactam"},
        "note": (
            "aztreonam/aztreonam_avibactam dotted (pass 1): "
            "no activity against Gram-positive cocci (potency 0.000/0.010)."
        ),
    },
    "Streptococcus agalactiae": {
        "drugs": {"aztreonam", "aztreonam_avibactam"},
        "note": (
            "aztreonam/aztreonam_avibactam dotted (pass 1): "
            "no activity against Gram-positive cocci (potency 0.000/0.010)."
        ),
    },
    # Nitrofurantoin: Proteus/Morganella/Klebsiella/Serratia are intrinsically
    # resistant via alkaline pH of urine + loss of nitroreductase NfsAB;
    # clinical guidelines contra-indicate nitrofurantoin for these organisms.
    "Klebsiella pneumoniae": {
        "drugs": {"nitrofurantoin"},
        "note": "nitrofurantoin dotted (pass 1): intrinsic resistance (potency 0.050 - not used for Klebsiella treatment).",
    },
    "Morganella spp.": {
        "drugs": {"nitrofurantoin"},
        "note": "nitrofurantoin dotted (pass 1): intrinsic resistance (potency 0.050).",
    },
    "Proteus spp.": {
        "drugs": {"nitrofurantoin"},
        "note": "nitrofurantoin dotted (pass 1): intrinsic resistance (potency 0.050).",
    },
    "Serratia spp.": {
        "drugs": {"nitrofurantoin"},
        "note": "nitrofurantoin dotted (pass 1): intrinsic resistance (potency 0.050).",
    },
    # Listeria: aztreonam intrinsically inactive; meropenem_vaborbactam not used.
    "Listeria monocytogenes": {
        "drugs": {"aztreonam_avibactam", "meropenem_vaborbactam"},
        "note": (
            "aztreonam_avibactam/meropenem_vaborbactam dotted (pass 1): "
            "not used for Listeria treatment (potency 0.010/0.050)."
        ),
    },
    # Yersinia: chromosomal AmpC (blaA gene) confers intrinsic ampicillin resistance.
    "Yersinia enterocolitica": {
        "drugs": {"ampicillin"},
        "note": "ampicillin dotted (pass 1): chromosomal AmpC intrinsic resistance (potency 0.020).",
    },
    # Legionella: beta-lactam/BLI combinations have poor activity; not used.
    "Legionella pneumophila": {
        "drugs": {"amoxicillin_clavulanate"},
        "note": (
            "amoxicillin_clavulanate dotted (pass 1): "
            "beta-lactam/BLI not active against Legionella (potency 0.050)."
        ),
    },
    # H. pylori: cefiderocol/ceftolozane not used for H. pylori treatment.
    "Helicobacter pylori": {
        "drugs": {"cefiderocol", "ceftolozane_tazobactam"},
        "note": (
            "cefiderocol/ceftolozane_tazobactam dotted (pass 1): "
            "not used for H. pylori treatment (potency 0.050)."
        ),
    },
    # Shigella: clindamycin has zero intrinsic activity against Gram-negatives
    # (outer membrane barrier; spectrum is Gram-positive/anaerobe only).
    "Shigella spp.": {
        "drugs": {"clindamycin"},
        "note": "clindamycin dotted (pass 1): no activity against Gram-negatives (potency 0.000).",
    },
}

# Read
with open(PATH, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames[:]
    rows = list(reader)

notes_col = next((fn for fn in fieldnames if fn.lower() == "notes"), None)

total_changes = 0
rows_changed = 0
for row in rows:
    bacteria = row.get("Bacteria", "").strip()
    config = PASS1.get(bacteria)
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
        print(f"  {bacteria}: dotted {len(dotted)} drug(s): {', '.join(sorted(dotted))}")

print(f"\nTotal: {total_changes} cells dotted across {rows_changed} bacteria rows.")

# Write back (utf-8-sig preserves BOM)
with open(PATH, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, quoting=csv.QUOTE_MINIMAL)
    writer.writeheader()
    writer.writerows(rows)

print("Written back to", PATH)
