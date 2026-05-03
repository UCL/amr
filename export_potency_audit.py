from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = REPO_ROOT / "potency_audit_matrix.csv"
POPULATION_RS = REPO_ROOT / "src" / "simulation" / "population.rs"
def run_appendix_dump() -> str:
    result = subprocess.run(
        ["cargo", "run", "--quiet", "--bin", "dump_parameter_appendix"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def extract_b4_table(markdown: str) -> list[dict[str, str]]:
    lines = markdown.splitlines()
    in_section = False
    table_lines: list[str] = []

    for line in lines:
        if line.startswith("### B.4 ") and "Potency Matrix" in line:
            in_section = True
            continue
        if in_section and line.startswith("### "):
            break
        if in_section and line.startswith("|"):
            table_lines.append(line)

    if len(table_lines) < 3:
        raise RuntimeError("Could not find B.4 potency matrix table in appendix output")

    headers = [cell.strip() for cell in table_lines[0].strip("|").split("|")]
    rows: list[dict[str, str]] = []
    for line in table_lines[2:]:
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        if len(cells) != len(headers):
            continue
        rows.append(dict(zip(headers, cells)))
    return rows


def parse_rust_string_array(file_text: str, const_name: str) -> list[str]:
    pattern = re.compile(
        rf"pub const {re.escape(const_name)}:[^=]*=\s*&?\[(.*?)\];",
        re.DOTALL,
    )
    match = pattern.search(file_text)
    if not match:
        raise RuntimeError(f"Could not find Rust array for {const_name}")
    body = match.group(1)
    return re.findall(r'"([^"]+)"', body)


def load_population_axes() -> tuple[list[str], list[str]]:
    file_text = POPULATION_RS.read_text(encoding="utf-8")
    bacteria = parse_rust_string_array(file_text, "BACTERIA_LIST")
    drugs = parse_rust_string_array(file_text, "DRUG_SHORT_NAMES")
    return bacteria, drugs


def classify_bucket(potency: float) -> str:
    if potency == 0.0:
        return "inactive"
    if potency > 1.0:
        return "above-scale"
    if potency <= 0.05:
        return "very-poor-or-none"
    if potency < 0.25:
        return "poor"
    if potency < 0.50:
        return "moderate"
    if potency < 1.0:
        return "good"
    return "excellent"


def review_flags(potency: float, init_multiplier: float) -> str:
    flags: list[str] = []
    if potency > 1.0:
        flags.append("above_scale")
    if potency == 0.1:
        flags.append("default_potency")
    if potency == 0.0:
        flags.append("intrinsic_inactivity_or_zero")
    if init_multiplier != 1.0:
        flags.append("nondefault_init_multiplier")
    return ";".join(flags)


def build_audit_rows(
    rows: list[dict[str, str]], bacteria_list: list[str], drug_list: list[str]
) -> list[dict[str, str]]:
    appendix_lookup = {
        (row["Bacteria"], row["Drug"]): row
        for row in rows
    }
    audit_rows: list[dict[str, str]] = []
    for bacteria in bacteria_list:
        for drug in drug_list:
            row = appendix_lookup.get((bacteria, drug))
            if row is None:
                potency = 0.0
                init_multiplier = 1.0
                source = "omitted_from_appendix_zero_potency_default_init"
            else:
                potency = float(row["Potency (no R)"])
                init_multiplier = float(row["Init multiplier"])
                source = "appendix_dump"

            audit_rows.append(
                {
                    "bacteria": bacteria,
                    "drug": drug,
                    "potency_no_r": f"{potency:.6f}",
                    "init_multiplier": f"{init_multiplier:.6f}",
                    "potency_bucket": classify_bucket(potency),
                    "review_flags": review_flags(potency, init_multiplier),
                    "source": source,
                    "review_status": "",
                    "evidence_notes": "",
                    "decision": "",
                }
            )
    return audit_rows


def write_csv(rows: list[dict[str, str]], output_path: Path) -> None:
    if not rows:
        raise RuntimeError("No potency rows were extracted")
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    markdown = run_appendix_dump()
    rows = extract_b4_table(markdown)
    bacteria_list, drug_list = load_population_axes()
    audit_rows = build_audit_rows(rows, bacteria_list, drug_list)
    write_csv(audit_rows, DEFAULT_OUTPUT)

    total_rows = len(audit_rows)
    above_scale = sum(1 for row in audit_rows if "above_scale" in row["review_flags"])
    defaults = sum(1 for row in audit_rows if "default_potency" in row["review_flags"])
    nondefault_init = sum(
        1 for row in audit_rows if "nondefault_init_multiplier" in row["review_flags"]
    )

    print(f"Wrote {total_rows} potency rows to {DEFAULT_OUTPUT.name}")
    print(f"Above-scale rows: {above_scale}")
    print(f"Default 0.1 potency rows: {defaults}")
    print(f"Rows with non-default initiation multipliers: {nondefault_init}")


if __name__ == "__main__":
    main()