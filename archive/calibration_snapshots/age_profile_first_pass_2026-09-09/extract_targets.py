"""Reproduce complete-continent age targets from UN WPP2024, 1 July 2023.

Usage: python extract_targets.py /path/to/WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx

Uses only the standard library. Writes age_targets.csv,
un_source_age_counts.csv and un_source_manifest.json beside this script.
North America is Northern America + Central America + Caribbean; all six
continents are checked against the workbook's World row at every single age.
Source location codes are read from the selected names, never guessed.
"""

import csv
import hashlib
import json
import re
import sys
import zipfile
from decimal import Decimal
from pathlib import Path
from xml.etree import ElementTree as ET


NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
SOURCE_URL = (
    "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/"
    "EXCEL_FILES/2_Population/WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx"
)
EXPECTED_SHA256 = "78109efb998fee10bd0c501844e04cafd0922836239e00aa90d9f35fc828a94a"
EXPECTED_BYTES = 220318066
YEAR = "2023"
COMPONENTS = {
    "north_america": ("Northern America", "Central America", "Caribbean"),
    "south_america": ("South America",),
    "europe": ("Europe",),
    "asia": ("Asia",),
    "africa": ("Africa",),
    "oceania": ("Oceania",),
}
BANDS = (
    ("0_5", 0, 5), ("6_14", 6, 14), ("15_49", 15, 49),
    ("50_79", 50, 79), ("80plus", 80, 100),
)


def cell_value(cell, strings):
    value = cell.find(NS + "v")
    if value is None:
        return "".join(cell.itertext())
    if cell.attrib.get("t") == "s":
        return strings[int(value.text)]
    return value.text


def extract(source):
    wanted = {name.casefold() for names in COMPONENTS.values() for name in names}
    wanted.add("world")
    selected = {}
    with zipfile.ZipFile(source) as workbook:
        strings = []
        with workbook.open("xl/sharedStrings.xml") as stream:
            for _, element in ET.iterparse(stream, events=("end",)):
                if element.tag == NS + "si":
                    strings.append("".join(element.itertext()))
                    element.clear()
        relations = ET.fromstring(workbook.read("xl/_rels/workbook.xml.rels"))
        targets = {element.attrib["Id"]: element.attrib["Target"] for element in relations}
        sheets = ET.fromstring(workbook.read("xl/workbook.xml")).find(NS + "sheets")
        estimates = next(sheet for sheet in sheets if sheet.attrib["name"] == "Estimates")
        relation = estimates.attrib[
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
        ]
        target = targets[relation]
        sheet_path = target.lstrip("/") if target.startswith("/") else "xl/" + target
        header = None
        with workbook.open(sheet_path) as stream:
            for _, element in ET.iterparse(stream, events=("end",)):
                if element.tag != NS + "row":
                    continue
                cells = {
                    re.sub(r"\d+$", "", cell.attrib["r"]): cell_value(cell, strings)
                    for cell in element
                }
                if cells.get("E") == "Location code" and cells.get("K") == "Year":
                    header = cells
                name = cells.get("C", "").strip()
                key = name.casefold()
                if key in wanted and cells.get("K") == YEAR:
                    if header is None or cells.get("B") != "Estimates":
                        raise ValueError(f"Unexpected source layout or variant for {name}")
                    if key in selected:
                        raise ValueError(f"Duplicate source region: {name}")
                    ages = {
                        int(age.rstrip("+")): Decimal(cells[column])
                        for column, age in header.items()
                        if age.rstrip("+").isdigit()
                    }
                    if set(ages) != set(range(101)):
                        raise ValueError(f"Incomplete single-age coverage for {name}")
                    if not all(value.is_finite() and value >= 0 for value in ages.values()):
                        raise ValueError(f"Invalid population counts for {name}")
                    selected[key] = {
                        "name": name, "code": cells["E"],
                        "source_row": element.attrib["r"], "ages": ages,
                    }
                element.clear()
                if set(selected) == wanted:
                    break
    missing = wanted.difference(selected)
    if missing:
        raise ValueError(f"Missing source regions: {sorted(missing)}")
    if len({row["code"] for row in selected.values()}) != len(selected):
        raise ValueError("Selected source regions have duplicate location codes")
    return selected


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    source = Path(sys.argv[1])
    with source.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if digest != EXPECTED_SHA256 or source.stat().st_size != EXPECTED_BYTES:
        raise ValueError("Workbook identity does not match the recorded WPP2024 source")
    selected = extract(source)
    regional = {
        region: {
            age: sum(selected[name.casefold()]["ages"][age] for name in names)
            for age in range(101)
        }
        for region, names in COMPONENTS.items()
    }
    world = selected["world"]["ages"]
    differences = {
        age: sum(ages[age] for ages in regional.values()) - world[age]
        for age in range(101)
    }
    # Allow at most five persons per age for independently rounded source
    # aggregates, recording the actual discrepancy even when it is zero.
    tolerance = Decimal("0.005")  # Source units are thousands of persons.
    if any(abs(value) > tolerance for value in differences.values()):
        raise ValueError("Six-continent age counts do not partition the World source row")
    world_total = sum(world.values())
    rows = []
    for region, names in COMPONENTS.items():
        ages = regional[region]
        total = sum(ages.values())
        if total <= 0:
            raise ValueError(f"Non-positive regional population for {region}")
        components = [selected[name.casefold()] for name in names]
        for band, lower, upper in BANDS:
            count = sum(ages[age] for age in range(lower, upper + 1))
            rows.append({
                "region": region, "age_group": band,
                "age_share": format(count / total, ".15f"),
                "age_population": str(count * 1000),
                "regional_total_population": str(total * 1000),
                "world_population_share": format(total / world_total, ".15f"),
                "year": YEAR, "sex": "both", "source_revision": "WPP2024",
                "source_sheet": "Estimates", "source_date": "2023-07-01",
                "source_regions": "|".join(item["name"] for item in components),
                "un_location_codes": "|".join(item["code"] for item in components),
                "source_rows": "|".join(item["source_row"] for item in components),
                "source_url": SOURCE_URL,
            })
        share_total = sum(Decimal(row["age_share"]) for row in rows if row["region"] == region)
        if abs(share_total - 1) > Decimal("0.000000000000003"):
            raise ValueError(f"Age shares do not sum to one for {region}")
        print(region, "population:", total * 1000, "age shares:",
              [row["age_share"] for row in rows if row["region"] == region])

    raw_rows = []
    for key in sorted(selected):
        item = selected[key]
        for age, count in sorted(item["ages"].items()):
            raw_rows.append({
                "source_region": item["name"], "un_location_code": item["code"],
                "year": YEAR, "age": "100+" if age == 100 else age,
                "population_thousands": str(count), "source_row": item["source_row"],
            })
    manifest = {
        "source_url": SOURCE_URL, "source_sha256": digest,
        "source_bytes": source.stat().st_size, "retrieved_date": "2026-09-09",
        "source_sheet": "Estimates", "source_date": "2023-07-01", "year": int(YEAR),
        "source_units": "thousands of persons", "output_population_units": "persons",
        "age_units": "completed years; 100+ open interval",
        "regional_components": {region: list(names) for region, names in COMPONENTS.items()},
        "selected_source_regions": [
            {"name": item["name"], "code": item["code"], "source_row": item["source_row"]}
            for key, item in sorted(selected.items())
        ],
        "source_age_rows": len(raw_rows), "target_rows": len(rows),
        "world_population": str(world_total * 1000),
        "six_continent_population": str(sum(sum(ages.values()) for ages in regional.values()) * 1000),
        "max_single_age_partition_difference_persons": str(max(abs(value) for value in differences.values()) * 1000),
        "total_partition_difference_persons": str(sum(differences.values()) * 1000),
        "single_age_partition_tolerance_persons": str(tolerance * 1000),
    }
    destination = Path(__file__).resolve().parent
    write_csv(destination / "age_targets.csv", rows)
    write_csv(destination / "un_source_age_counts.csv", raw_rows)
    (destination / "un_source_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print("World partition verified; wrote", len(rows), "targets and", len(raw_rows), "source age counts.")


if __name__ == "__main__":
    main()
