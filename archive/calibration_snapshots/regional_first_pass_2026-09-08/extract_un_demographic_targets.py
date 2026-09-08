"""Reproduce six regional age targets from the official UN WPP2024 workbook.

Usage: python extract_un_demographic_targets.py /path/to/WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx
Only Python's standard library is needed. Outputs are written beside this script.
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
SOURCE_URL = "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/EXCEL_FILES/2_Population/WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx"
REGIONS = {
    "905": ("north_america", "Northern America"),
    "931": ("south_america", "South America"),
    "908": ("europe", "Europe"),
    "935": ("asia", "Asia"),
    "903": ("africa", "Africa"),
    "909": ("oceania", "Oceania"),
}
BANDS = [("0_5", 0, 5), ("6_14", 6, 14), ("15_49", 15, 49),
         ("50_79", 50, 79), ("80plus", 80, 100)]
YEAR = "2023"


def cell_value(cell, strings):
    value = cell.find(NS + "v")
    if value is None:
        return "".join(cell.itertext())
    if cell.attrib.get("t") == "s":
        return strings[int(value.text)]
    return value.text


def main():
    source = Path(sys.argv[1])
    destination = Path(__file__).resolve().parent
    strings = []
    selected = {}
    with zipfile.ZipFile(source) as workbook:
        for _, element in ET.iterparse(workbook.open("xl/sharedStrings.xml"), events=("end",)):
            if element.tag == NS + "si":
                strings.append("".join(element.itertext()))
                element.clear()
        relations = ET.fromstring(workbook.read("xl/_rels/workbook.xml.rels"))
        targets = {element.attrib["Id"]: element.attrib["Target"] for element in relations}
        sheets = ET.fromstring(workbook.read("xl/workbook.xml")).find(NS + "sheets")
        estimates = next(sheet for sheet in sheets if sheet.attrib["name"] == "Estimates")
        relation = estimates.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
        sheet_path = "xl/" + targets[relation]
        header = None
        for _, element in ET.iterparse(workbook.open(sheet_path), events=("end",)):
            if element.tag != NS + "row":
                continue
            cells = {re.sub(r"\d+$", "", cell.attrib["r"]): cell_value(cell, strings)
                     for cell in element}
            if cells.get("E") == "Location code" and cells.get("K") == "Year":
                header = cells
            code = cells.get("E")
            if code in REGIONS and cells.get("K") == YEAR:
                assert header is not None
                assert cells["C"] == REGIONS[code][1], (code, cells["C"])
                assert cells["B"] == "Estimates"
                assert code not in selected
                ages = {int(age.rstrip("+")): Decimal(cells[column])
                        for column, age in header.items()
                        if age.rstrip("+").isdigit()}
                assert set(ages) == set(range(101))
                assert all(number >= 0 for number in ages.values())
                selected[code] = {"ages": ages, "source_row": element.attrib["r"]}
            element.clear()
            if len(selected) == len(REGIONS):
                break
    assert set(selected) == set(REGIONS)
    rows = []
    raw_rows = []
    for code, (region, source_region) in REGIONS.items():
        data = selected[code]
        ages = data["ages"]
        total = sum(ages.values())
        for age, count in ages.items():
            raw_rows.append({"region": region, "source_region": source_region,
                "un_location_code": code, "year": YEAR,
                "age": "100+" if age == 100 else age,
                "population_thousands": str(count), "source_row": data["source_row"]})
        for band, lower, upper in BANDS:
            count = sum(ages[age] for age in range(lower, upper + 1))
            rows.append({"region": region, "age_group": band,
                "age_share": format(count / total, ".15f"),
                "age_population": str(count * 1000),
                "regional_total_population": str(total * 1000),
                "year": YEAR, "sex": "both", "un_location_code": code,
                "source_region": source_region, "source_revision": "WPP2024",
                "source_sheet": "Estimates", "source_row": data["source_row"],
                "source_url": SOURCE_URL})
        assert abs(sum(Decimal(row["age_share"]) for row in rows if row["region"] == region) - 1) < Decimal("0.000000000000003")
        print(region, "total persons:", total * 1000,
              "shares:", [row["age_share"] for row in rows if row["region"] == region])
    for filename, contents in [("demographic_targets_un_wpp2024.csv", rows),
                               ("demographic_source_single_ages_un_wpp2024.csv", raw_rows)]:
        with (destination / filename).open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=list(contents[0]))
            writer.writeheader()
            writer.writerows(contents)
    manifest = {"source_url": SOURCE_URL, "source_sha256": hashlib.file_digest(source.open("rb"), "sha256").hexdigest(),
                "source_bytes": source.stat().st_size, "retrieved_date": "2026-09-08",
                "source_sheet": "Estimates", "year": int(YEAR), "source_units": "thousands of persons",
                "age_units": "completed years; 100+ open interval", "age_rows": len(raw_rows), "target_rows": len(rows)}
    (destination / "demographic_source_un_wpp2024.manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
