import csv
import json
import unittest
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

from amr_simulation_output_analysis.build_resistance_targets_v1 import (
    PREVALENCE_COMPONENT,
    SEVERITY_COMPONENT,
    TARGET_COLUMNS,
    TARGET_SET_VERSION,
    build_resistance_targets_v1,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
TARGET_PATH = DATA / "resistance_targets_v1.csv"
SOURCE_PATH = DATA / "resistance_target_sources_v1.csv"


def _read_csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_wide_values(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        drugs = [name for name in reader.fieldnames[1:] if name != "notes"]
        values = {}
        for row in reader:
            for drug in drugs:
                values[(row["Bacteria"], drug)] = "" if row[drug] == "." else row[drug]
    return drugs, values


def _row_key(row):
    return row["component"], row["bacteria"], row["drug"]


class ResistanceTargetSchemaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rows = _read_csv(TARGET_PATH)
        cls.sources = _read_csv(SOURCE_PATH)

    def test_generated_files_are_reproducible(self) -> None:
        with TemporaryDirectory() as temp_dir:
            target_output = Path(temp_dir) / TARGET_PATH.name
            source_output = Path(temp_dir) / SOURCE_PATH.name
            build_resistance_targets_v1(
                ROOT,
                target_output=target_output,
                source_output=source_output,
            )

            self.assertEqual(target_output.read_bytes(), TARGET_PATH.read_bytes())
            self.assertEqual(source_output.read_bytes(), SOURCE_PATH.read_bytes())

    def test_long_form_values_exactly_match_wide_matrices(self) -> None:
        prevalence_drugs, prevalence = _read_wide_values(
            DATA / "resistance_prevalence_values.csv"
        )
        severity_drugs, severity = _read_wide_values(
            DATA / "resistance_average_resistant_values.csv"
        )
        self.assertEqual(prevalence_drugs, severity_drugs)

        long_values = {_row_key(row): row["value"] for row in self.rows}
        self.assertEqual(len(long_values), len(self.rows))
        self.assertEqual(len(self.rows), 42 * 61 * 2)

        for (bacterium, drug), value in prevalence.items():
            self.assertEqual(
                long_values[(PREVALENCE_COMPONENT, bacterium, drug)], value
            )
        for (bacterium, drug), value in severity.items():
            self.assertEqual(
                long_values[(SEVERITY_COMPONENT, bacterium, drug)], value
            )

    def test_schema_and_provenance_contract_is_complete(self) -> None:
        with TARGET_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
            self.assertEqual(csv.DictReader(handle).fieldnames, TARGET_COLUMNS)

        schema = json.loads(
            (DATA / "resistance_targets_v1.schema.json").read_text(encoding="utf-8")
        )
        self.assertEqual(schema["required"], TARGET_COLUMNS)
        self.assertEqual(set(schema["properties"]), set(TARGET_COLUMNS))

        source_ids = [row["source_id"] for row in self.sources]
        self.assertEqual(len(source_ids), 43)
        self.assertEqual(len(source_ids), len(set(source_ids)))
        known_sources = set(source_ids)

        for row in self.rows:
            self.assertEqual(row["target_set_version"], TARGET_SET_VERSION)
            if row["source_id"]:
                self.assertIn(row["source_id"], known_sources)
            self.assertEqual(row["evidence_weight"], "")
            if row["include_in_score"] == "true":
                self.assertNotEqual(row["value"], "")
                self.assertEqual(row["score_exclusion_reason"], "")
                self.assertEqual(row["score_row_weight"], "1.0")
            else:
                self.assertNotEqual(row["score_exclusion_reason"], "")
                self.assertEqual(row["score_row_weight"], "0.0")

    def test_cell_statuses_make_legacy_missingness_explicit(self) -> None:
        counts = Counter((row["component"], row["cell_status"]) for row in self.rows)
        self.assertEqual(counts[(PREVALENCE_COMPONENT, "active_target")], 1294)
        self.assertEqual(
            counts[(PREVALENCE_COMPONENT, "legacy_unclassified_missing")], 1268
        )
        self.assertEqual(counts[(SEVERITY_COMPONENT, "active_target")], 1292)
        self.assertEqual(
            counts[(SEVERITY_COMPONENT, "inactive_legacy_prevalence_gate")], 118
        )
        self.assertEqual(
            counts[(SEVERITY_COMPONENT, "legacy_unclassified_missing")], 1152
        )

    def test_include_flags_capture_target_side_score_eligibility(self) -> None:
        included = Counter(
            row["component"]
            for row in self.rows
            if row["include_in_score"] == "true"
        )
        self.assertEqual(included[PREVALENCE_COMPONENT], 1241)
        self.assertEqual(included[SEVERITY_COMPONENT], 1240)

        allowed_reasons = {
            "legacy_prevalence_target_missing",
            "severity_target_missing",
            "hard_exclusion_rifampicin",
            "hard_exclusion_mdr_tb",
            "hard_exclusion_listeria",
        }
        for row in self.rows:
            reasons = {
                reason
                for reason in row["score_exclusion_reason"].split(";")
                if reason
            }
            self.assertTrue(reasons.issubset(allowed_reasons))


if __name__ == "__main__":
    unittest.main()
