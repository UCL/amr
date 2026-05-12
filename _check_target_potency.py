"""Check for resistance targets where drug potency is below the selection threshold."""
import csv
import re

THRESHOLD = 0.15

def slug(s):
    s = s.lower().strip()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    return s.strip("_")

# Load audit matrix
potency = {}
sim_bacteria = set()
with open("potency_audit_matrix.csv", newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        b = slug(row["bacteria"])
        d = row["drug"].strip().lower()
        sim_bacteria.add(b)
        try:
            potency[(b, d)] = float(row["potency_no_r"])
        except (ValueError, KeyError):
            pass

# Load targets (utf-8-sig strips BOM from first column name)
hits = []
with open("data/resistance_prevalence_values.csv", newline="", encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        bacteria = row.get("Bacteria", "").strip()
        b_slug = slug(bacteria)
        if b_slug not in sim_bacteria:
            continue
        for drug, val in row.items():
            if drug in ("Bacteria", "notes", "Notes", "note", "Note", ""):
                continue
            val = val.strip()
            if val in (".", ""):
                continue
            try:
                target = float(val)
            except ValueError:
                continue
            p = potency.get((b_slug, drug.strip().lower()))
            if p is not None and p < THRESHOLD:
                hits.append((bacteria, drug, target, p))

hits.sort(key=lambda x: (x[0], x[1]))
print(f"Pairs with numeric target AND potency < {THRESHOLD} (modelled bacteria only): {len(hits)}")
print()
header = f"{'Bacteria':<48} {'Drug':<30} {'Target':>7} {'Potency':>8}"
print(header)
print("-" * len(header))
for bacteria, drug, target, p in hits:
    flag = "  <-- likely intrinsic" if p <= 0.05 else ""
    print(f"{bacteria:<48} {drug:<30} {target:>7.0%} {p:>8.3f}{flag}")
